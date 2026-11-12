"""Has anything new been published that could join this benchmark?

One deterministic check against a public, keyless API (Europe PMC), compared
with a watermark committed in reference/source_watermarks.json.

No LLM, no API key, no extraction. It cannot add or invent a row -- it only says
"go look". Ingestion stays a human decision, and the row that eventually lands is
written by scripts/append_paper.py, never by a model.

On the query, and why it is tight. This field is small: de-novo designed proteins
with a MEASURED affinity for a small molecule are a few papers a year, not a few
a week. The obvious query ("de novo" AND "small molecule") is dominated by the
opposite problem -- designing a small-molecule drug against a natural protein --
and ran at roughly 2 relevant hits in 23. Requiring the DESIGNED thing to be the
protein and the BOUND thing to be the ligand gives 5 hits over three years, of
which 2 are real. Adding biosensor/sensor terms for recall took it to 9 hits and
brought back SPR instrument papers and imprinted polymers, so it is left out.

Recall is therefore the weak side, deliberately. Treat this as a reading list,
not a verdict, and keep watching the handful of groups that publish here
directly -- a quarterly skim of their output will find things this misses.

  python scripts/check_sources.py            # report
  python scripts/check_sources.py --update   # also write the new watermark
"""

from __future__ import annotations

import argparse
import json
import os
import re
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
WATERMARK = ROOT / "reference" / "source_watermarks.json"
API = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"

DESIGN = ('(TITLE_ABS:"de novo design" OR TITLE_ABS:"de novo designed" '
          'OR TITLE_ABS:"computationally designed" OR TITLE_ABS:"designed miniprotein" '
          'OR TITLE_ABS:"designed protein" OR TITLE_ABS:"designed binding protein")')
BINDS = ('(TITLE_ABS:"small-molecule binding" OR TITLE_ABS:"small molecule binding" '
         'OR TITLE_ABS:"ligand-binding protein" OR TITLE_ABS:"binding pocket" '
         'OR TITLE_ABS:"ligand binder" OR TITLE_ABS:"small-molecule binder")')
AFFINITY = ('(ABSTRACT:"dissociation constant" OR ABSTRACT:"binding affinity" '
            'OR ABSTRACT:"Kd" OR ABSTRACT:"KD")')

# Europe PMC indexes weeks late, so a window shorter than this misses papers that
# were published before the last run but only became searchable after it.
LOOKBACK_DAYS = 120

AFFINITY_NUM = re.compile(
    r"\b\d+(?:\.\d+)?\s*(?:[pnµumM]M|nanomolar|micromolar|picomolar|millimolar)\b", re.I)
ASSAY = re.compile(
    r"\b(SPR|BLI|ITC|biolayer|isothermal titration|fluorescence polari|MST|"
    r"microscale thermophoresis|surface plasmon)\b", re.I)
LIGAND_HINT = re.compile(
    r"\b(ligand|small[- ]molecule|cofactor|metabolite|steroid|heme|haem)\b", re.I)


TAGS = re.compile(r"<[^>]+>")


def clean(text: str) -> str:
    """Europe PMC returns titles with markup (`<i>De novo</i>`), which would land
    verbatim in a GitHub issue."""
    return TAGS.sub("", text or "").replace("  ", " ").strip()


def _search(query: str, page_size: int = 50) -> list[dict]:
    url = API + "?" + urllib.parse.urlencode(
        {"query": query, "format": "json", "pageSize": page_size, "resultType": "core"})
    with urllib.request.urlopen(url, timeout=60) as fh:
        return json.load(fh)["resultList"]["result"]


def rank(paper: dict) -> int:
    """0-3 by how likely the abstract is to yield an extractable row.

    A paper with a number AND an assay is one somebody can pull a (protein,
    ligand, KD) triple out of; one with neither is a reading-list entry at best.
    """
    text = f"{paper.get('title', '')} {paper.get('abstractText', '') or ''}"
    return ((2 if AFFINITY_NUM.search(text) else 0)
            + (1 if ASSAY.search(text) else 0)) if LIGAND_HINT.search(text) else 0


def main(update: bool) -> None:
    since = os.environ.get("SINCE") or _since_date()
    query = f"{DESIGN} AND {BINDS} AND {AFFINITY} AND FIRST_PDATE:[{since} TO 2099-12-31]"
    try:
        papers = _search(query)
    except Exception as exc:                       # a flaky API must not fail the cron
        print(f"Europe PMC unreachable: {exc}")
        _emit(changed=False, summary=f"Europe PMC unreachable: {exc}")
        return

    seen = json.loads(WATERMARK.read_text()) if WATERMARK.is_file() else {}
    known = set(seen.get("europepmc_ids", []))
    fresh = [p for p in papers if p.get("id") and p["id"] not in known]
    fresh.sort(key=rank, reverse=True)

    print(f"Europe PMC: {len(papers)} in window since {since}, {len(fresh)} not seen before")
    lines = []
    for p in fresh:
        r = rank(p)
        doi = p.get("doi") or p.get("id")
        line = f"[{r}/3] {p.get('firstPublicationDate', '?')[:7]}  {clean(p.get('title', ''))[:90]}  {doi}"
        print("  " + line)
        lines.append(line)

    if update:
        WATERMARK.parent.mkdir(parents=True, exist_ok=True)
        WATERMARK.write_text(json.dumps(
            {"europepmc_ids": sorted(known | {p["id"] for p in papers if p.get("id")}),
             "last_window_start": since}, indent=2) + "\n")
        print(f"wrote {WATERMARK.relative_to(ROOT)}")

    summary = ("\n".join(lines) if lines
               else "No paper in the window that was not already seen.")
    # the best new candidate, for the ingest job; only if it looks in scope
    best = next((p for p in fresh if rank(p) >= 2), None)
    top = (best.get("doi") or "") if best else ""
    if top:
        print(f"top candidate for ingest: {top}")
    _emit(changed=bool(fresh), summary=summary, top_doi=top)


def _since_date() -> str:
    import datetime as dt

    return (dt.date.today() - dt.timedelta(days=LOOKBACK_DAYS)).isoformat()


def _emit(changed: bool, summary: str, top_doi: str = "") -> None:
    """top_doi is the single best new candidate, for the ingest job to chain to.
    Empty unless something scored >= 2, because a 0 is usually small-molecule
    DRUG design -- the opposite direction -- and is not worth a model's tokens."""
    out = os.environ.get("GITHUB_OUTPUT")
    if not out:
        return
    with open(out, "a") as fh:
        fh.write(f"changed={'true' if changed else 'false'}\n")
        fh.write(f"top_doi={top_doi}\n")
        fh.write("summary<<EOF\n" + summary + "\nEOF\n")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--update", action="store_true")
    main(ap.parse_args().update)
