"""Extract one paper's rows with the Anthropic API, for the unattended path.

The interactive path is the `add-smallmol-paper` skill: subagents that can browse,
open a PDF and chase an SI. This is the same two-step contract -- extract, then
adversarially verify -- driven from CI, where there is no browser.

WHAT THIS DELIBERATELY CANNOT DO, and the reason the loop still ends at a human:

  * It reads only what Europe PMC will serve as text: the open-access full text,
    or failing that the abstract. Most of this field's sequences live in a PDF
    supplementary table behind a publisher paywall, so for a closed paper this
    will correctly extract NOTHING rather than guess. "Zero rows" is the normal
    and correct outcome for maybe half of candidates; it is not a failure.
  * It never writes `reference/`. It writes a candidate CSV, `append_paper.py`
    ingests it, and the workflow opens a PULL REQUEST. No model output reaches
    main without a person merging it, which is the property that makes this set
    worth citing. qc_dataset.py runs on that PR.
  * It refuses to invent. The prompts below are the agent definitions in
    `.claude/agents/`, read verbatim at runtime rather than copied, so the rules
    that bind the interactive path bind this one and cannot drift from it.

    python scripts/ingest_paper.py --doi 10.1038/s41586-025-xxxxx --out cand.csv
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import os
import re
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
AGENTS = ROOT / ".claude" / "agents"
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest"
COLS = ["protein", "ligand_name", "sequence", "smiles", "kd_nM",
        "is_binder", "source", "source_url", "notes"]
# Full text runs to ~100k characters; the cap keeps one paper to a predictable
# call rather than an unbounded bill on a 300-page supplement.
MAX_CHARS = 400_000
MODEL = "claude-sonnet-5-5"


def _get(url: str, timeout: int = 60) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "denovo-smallmol-benchmark"})
    with urllib.request.urlopen(req, timeout=timeout) as fh:   # noqa: S310 - fixed host
        return fh.read()


def agent_prompt(name: str) -> str:
    """The agent definition's body, minus its YAML frontmatter.

    Read at runtime on purpose: the extraction rules exist once, in the file the
    interactive path also uses. A copy here would drift the first time one was
    edited, and the drift would be invisible -- both paths would still run.
    """
    text = (AGENTS / f"{name}.md").read_text()
    if text.startswith("---"):
        text = text.split("---", 2)[2]
    return text.strip()


def fetch_text(doi: str) -> tuple[str, str, str]:
    """-> (text, source_url, how). `how` names what was actually readable, and
    is reported so a thin extraction is attributable to a thin source."""
    q = urllib.parse.quote(f'DOI:"{doi}"')
    hits = json.loads(_get(f"{EPMC}/search?query={q}&resultType=core&format=json"))
    results = hits.get("resultList", {}).get("result") or []
    if not results:
        return "", f"https://doi.org/{doi}", "not in Europe PMC"
    r = results[0]
    url = f"https://doi.org/{doi}"
    pmcid, oa = r.get("pmcid"), r.get("isOpenAccess") == "Y"
    if pmcid and oa:
        try:
            xml = _get(f"{EPMC}/{pmcid}/fullTextXML").decode("utf-8", "replace")
            # crude de-tagging: the model reads prose, and a real XML parse buys
            # nothing when the next step is a language model
            text = re.sub(r"<[^>]+>", " ", xml)
            text = re.sub(r"\s+", " ", text).strip()
            if len(text) > 500:
                return text[:MAX_CHARS], url, f"open-access full text ({pmcid})"
        except Exception as exc:                      # noqa: BLE001 - fall back to abstract
            print(f"  full text unavailable ({type(exc).__name__}), falling back")
    abstract = (r.get("abstractText") or "").strip()
    how = "ABSTRACT ONLY - closed access" if abstract else "no text available"
    return abstract, url, how


def ask(client, system: str, user: str) -> str:
    msg = client.messages.create(
        model=MODEL, max_tokens=8000,
        system=system, messages=[{"role": "user", "content": user}],
    )
    return "".join(b.text for b in msg.content if b.type == "text")


def parse_csv(text: str) -> list[dict]:
    """Pull the CSV out of a reply that may be fenced or prefaced with prose."""
    if "```" in text:
        blocks = re.findall(r"```(?:csv)?\s*\n(.*?)```", text, re.S)
        text = max(blocks, key=len) if blocks else text
    lines = [ln for ln in text.splitlines() if ln.strip()]
    start = next((i for i, ln in enumerate(lines) if ln.lower().startswith("protein,")), None)
    if start is None:
        return []
    rows = list(csv.DictReader(io.StringIO("\n".join(lines[start:]))))
    return [r for r in rows if (r.get("protein") or "").strip()]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--doi", required=True)
    ap.add_argument("--out", required=True, help="candidate CSV; NOT reference/")
    ap.add_argument("--model", default=MODEL)
    a = ap.parse_args()

    if not os.environ.get("ANTHROPIC_API_KEY"):
        print("ANTHROPIC_API_KEY is not set; nothing attempted.")
        return 1
    import anthropic

    client = anthropic.Anthropic()
    globals()["MODEL"] = a.model

    text, url, how = fetch_text(a.doi)
    print(f"source: {how}  ({len(text)} chars)")
    if not text:
        print("no readable text; writing an empty candidate file.")
        Path(a.out).write_text(",".join(COLS) + "\n")
        return 0

    task = (f"Paper DOI: {a.doi}\nCanonical source_url for every row: {url}\n"
            f"How much of it is readable here: {how}\n\n"
            "Return ONLY the CSV, with the header, and no prose around it. If the "
            "text below does not contain what a row needs -- most often the "
            "designed protein's own sequence, which usually lives in a PDF "
            "supplement this cannot see -- return the header alone. A missing row "
            f"is recoverable; a wrong one is not.\n\n---\n{text}")
    raw = parse_csv(ask(client, agent_prompt("smallmol-extractor"), task))
    print(f"extractor: {len(raw)} candidate row(s)")
    if not raw:
        Path(a.out).write_text(",".join(COLS) + "\n")
        return 0

    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=COLS, extrasaction="ignore")
    w.writeheader()
    w.writerows(raw)
    check = (f"Paper DOI: {a.doi}\nHow much of it was readable: {how}\n\n"
             "Below is the extraction, then the source text it came from. Drop "
             "every row the text does not support. Return ONLY the CSV.\n\n"
             f"EXTRACTION:\n{buf.getvalue()}\n\n---\nSOURCE:\n{text}")
    good = parse_csv(ask(client, agent_prompt("smallmol-verifier"), check))
    print(f"verifier: {len(good)} row(s) survived ({len(raw) - len(good)} dropped)")

    with Path(a.out).open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, extrasaction="ignore")
        w.writeheader()
        w.writerows(good)
    print(f"wrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
