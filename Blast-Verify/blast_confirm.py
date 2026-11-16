#!/usr/bin/env python3
"""
Per-read BLAST confirmation rate, for one file or many.

Each read carries a claimed species in its FASTQ header, formatted
"<read_id> <taxid> <species name>". For every read we take its single best
BLAST hit (highest bitscore) and ask: of the reads claiming species X, what
fraction have a best hit that actually IS X (or a strain within it)?

Single file:
    python blast_confirm.py --fastq corona_REY.fastq --blast corona_REY.blastn.nt.tsv \
        --only-taxid 28295 --match lineage

Many files + many viruses (reads the same targets.tsv the BLAST step used):
    python3 blast_confirm.py --config targets.xlsx --match lineage --min-qcov 95

Denominator = every targeted read (read straight from the FASTQ), so reads
with no BLAST hit are counted as unconfirmed rather than dropped.
"""

import argparse
import gzip
import os
import re
import sys

import pandas as pd

BLAST_COLS = ["qseqid", "sseqid", "pident", "staxid", "staxids", "sskingdom",
              "sscinames", "stitle", "qcovs", "evalue", "bitscore"]

PER_READ_COLS = ["read", "claimed", "best_hit", "best_taxid",
                 "pident", "qcov", "status"]


# ---------------------------------------------------------------------------
# Pull the claimed species + taxid from a header line (no leading '@').
# Format: "<read_id> <taxid> <species name...>"
# ---------------------------------------------------------------------------
def parse_claim(header: str):
    parts = header.split()
    if len(parts) >= 3 and parts[1].isdigit():
        return " ".join(parts[2:]).strip(), parts[1]   # species, taxid
    if len(parts) >= 2:
        return " ".join(parts[1:]).strip(), None
    return parts[0], None


def read_id(header: str) -> str:
    return header.split(None, 1)[0]


def norm(name: str) -> str:
    return re.sub(r"\s+", " ", name.strip().lower())


def names_match(claimed: str, hit_names: str) -> bool:
    c = norm(claimed)
    for h in re.split(r"[;,]", str(hit_names)):
        h = norm(h)
        if h and (c in h or h in c):
            return True
    return False


_LINEAGE = None        # lazy-loaded parent map from nodes.dmp
TAXDUMP_HINT = None     # optional path to nodes.dmp or a taxdump dir (--taxdump)


def _find_nodes_dmp():
    """Locate nodes.dmp: use --taxdump if given, else a cached copy, else
    download NCBI's taxdump once (~60 MB) into the cache. Stdlib only."""
    import tarfile
    import urllib.request

    cands = []
    if TAXDUMP_HINT:
        cands.append(os.path.join(TAXDUMP_HINT, "nodes.dmp")
                     if os.path.isdir(TAXDUMP_HINT) else TAXDUMP_HINT)
    cache = os.path.join(os.path.expanduser("~"), ".cache", "blast_confirm")
    cands.append(os.path.join(cache, "nodes.dmp"))
    for c in cands:
        if c and os.path.exists(c):
            return c

    os.makedirs(cache, exist_ok=True)
    url = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
    tgz = os.path.join(cache, "taxdump.tar.gz")
    sys.stderr.write(f"Downloading NCBI taxdump (~60 MB, one time) into {cache} ...\n")
    try:
        urllib.request.urlretrieve(url, tgz)
        with tarfile.open(tgz) as tf:
            tf.extract("nodes.dmp", cache)
    except Exception as e:
        raise SystemExit(
            f"Could not fetch the NCBI taxonomy ({e}).\n"
            "  Fixes: run where there's internet (login node/laptop), OR point\n"
            "  --taxdump at an existing nodes.dmp, OR use --match name.")
    return os.path.join(cache, "nodes.dmp")


class _Lineage:
    """child_taxid -> parent_taxid map from nodes.dmp, with lineage walk."""
    def __init__(self, nodes_path):
        self.parent = {}
        with open(nodes_path) as fh:
            for line in fh:
                p = line.split("\t|\t")
                if len(p) >= 2 and p[0].isdigit() and p[1].isdigit():
                    self.parent[int(p[0])] = int(p[1])

    def contains(self, ancestor, taxid):
        t, anc, seen = int(taxid), int(ancestor), set()
        while t in self.parent and t not in seen:
            if t == anc:
                return True
            seen.add(t)
            nxt = self.parent[t]
            if nxt == t:            # reached root
                break
            t = nxt
        return t == anc


def lineage_match(claim_taxid, hit_taxids) -> bool:
    """True if the claimed species' taxid is at or above the hit's taxid in the
    NCBI tree -- i.e. the hit is that species or a strain within it."""
    global _LINEAGE
    if _LINEAGE is None:
        _LINEAGE = _Lineage(_find_nodes_dmp())
    if claim_taxid in (None, "") or not str(claim_taxid).isdigit():
        return False
    claim = int(claim_taxid)
    for t in re.split(r"[;,]", str(hit_taxids)):
        t = t.strip()
        if t.isdigit() and _LINEAGE.contains(claim, int(t)):
            return True
    return False


def _open(path):
    """Open a FASTQ as text whether it's plain or gzipped."""
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def load_reads(fastq):
    claim_species, claim_taxid = {}, {}
    with _open(fastq) as fh:
        for i, line in enumerate(fh):
            if i % 4 != 0:
                continue
            header = line[1:].rstrip("\n").rstrip("\r")   # drop leading @, CRs
            rid = read_id(header)
            sp, tx = parse_claim(header)
            claim_species[rid], claim_taxid[rid] = sp, tx
    return claim_species, claim_taxid


def analyze_file(fastq, blast, only_taxid, only_species,
                 match, min_pident, min_qcov, sample_contains=("REY",)):
    claim_species, claim_taxid = load_reads(fastq)
    if not claim_species:
        raise SystemExit(f"No reads parsed from {fastq} -- check parse_claim.")

    # Keep reads whose header contains ANY of the sample tokens, each as a
    # delimited token so 'REY' won't match inside 'grey'. Empty list = keep all.
    toks = ([sample_contains] if isinstance(sample_contains, str)
            else list(sample_contains or []))
    toks = [t for t in toks if t]
    if toks:
        alt = "|".join(re.escape(t) for t in toks)
        pat = re.compile(r"(?:^|[^A-Za-z0-9])(?:" + alt + r")(?:[^A-Za-z0-9]|$)",
                         re.I)
        for dead in [r for r in list(claim_species) if not pat.search(r)]:
            claim_species.pop(dead); claim_taxid.pop(dead, None)
        if not claim_species:
            sys.stderr.write(f"note: no reads in {fastq} match sample "
                             f"{toks}\n")

    # Restrict to the viruses of interest for this file.
    if only_taxid or only_species:
        keep_tax = set(str(t) for t in (only_taxid or []))
        keep_sp = [s.lower() for s in (only_species or [])]

        def wanted(rid):
            if keep_tax and str(claim_taxid.get(rid)) in keep_tax:
                return True
            if keep_sp and any(s in claim_species[rid].lower() for s in keep_sp):
                return True
            return False

        for dead in [r for r in list(claim_species) if not wanted(r)]:
            claim_species.pop(dead); claim_taxid.pop(dead, None)

    best = pd.DataFrame()
    if os.path.exists(blast) and os.path.getsize(blast) > 0:
        # New TSVs carry a column header; older ones don't. Skip it if present.
        with open(blast) as fh:
            has_header = fh.readline().split("\t", 1)[0] == "qseqid"
        df = pd.read_csv(blast, sep="\t", names=BLAST_COLS, dtype=str,
                         skiprows=1 if has_header else 0)
        for c in ("pident", "qcovs", "bitscore"):
            df[c] = pd.to_numeric(df[c], errors="coerce")
        best = (df.sort_values("bitscore", ascending=False)
                  .drop_duplicates("qseqid", keep="first")
                  .set_index("qseqid"))

    rows = []
    for rid, claimed in claim_species.items():
        rec = {"read": rid, "claimed": claimed, "best_hit": None,
               "best_taxid": None, "pident": None, "qcov": None,
               "status": "no_hit"}
        if not best.empty and rid in best.index:
            hit = best.loc[rid]
            rec.update(best_hit=hit["sscinames"], best_taxid=hit["staxid"],
                       pident=hit["pident"], qcov=hit["qcovs"])
            passes = (float(hit["pident"] or 0) >= min_pident and
                      float(hit["qcovs"] or 0) >= min_qcov)
            if match == "taxid":
                ok = passes and claim_taxid.get(rid) not in (None, "") and \
                     str(claim_taxid[rid]) in re.split(r"[;,]", str(hit["staxids"]))
            elif match == "lineage":
                ok = passes and lineage_match(claim_taxid.get(rid), hit["staxids"])
            else:
                ok = passes and names_match(claimed, hit["sscinames"])
            rec["status"] = "confirmed" if ok else "mismatch"
        rows.append(rec)

    return pd.DataFrame(rows, columns=PER_READ_COLS)


def summarize(per_read):
    out = []
    for (fname, sp), grp in per_read.groupby(["file", "claimed"]):
        n = len(grp)
        conf = int((grp.status == "confirmed").sum())
        wrong = grp[grp.status == "mismatch"]["best_hit"].value_counts().head(3)
        out.append({
            "file": fname,
            "claimed_species": sp,
            "n_reads": n,
            "n_confirmed": conf,
            "pct_confirmed": round(100 * conf / n, 1) if n else 0.0,
            "n_mismatch": int((grp.status == "mismatch").sum()),
            "n_no_hit": int((grp.status == "no_hit").sum()),
            "top_mismatched_hits": "; ".join(f"{k} ({v})" for k, v in wrong.items()),
        })
    return pd.DataFrame(out).sort_values(["file", "n_reads"],
                                         ascending=[True, False])


def _read_xlsx_rows(path):
    """Read the first worksheet of an .xlsx using only the standard library
    (no pandas/openpyxl/xlrd). Yields a list of cell strings per row."""
    import zipfile
    import xml.etree.ElementTree as ET

    def loc(tag):
        return tag.rsplit("}", 1)[-1]

    def colnum(ref):
        m = re.match(r"[A-Z]+", ref or "")
        if not m:
            return None
        n = 0
        for ch in m.group(0):
            n = n * 26 + (ord(ch) - 64)
        return n

    z = zipfile.ZipFile(path)
    names = z.namelist()
    shared = []
    if "xl/sharedStrings.xml" in names:
        root = ET.fromstring(z.read("xl/sharedStrings.xml"))
        for si in root:
            shared.append("".join(t.text or "" for t in si.iter()
                                  if loc(t.tag) == "t"))
    sheets = sorted(n for n in names
                    if re.match(r"xl/worksheets/sheet\d+\.xml$", n))
    if not sheets:
        raise SystemExit(f"No worksheet found in {path}")
    root = ET.fromstring(z.read(sheets[0]))
    for row in root.iter():
        if loc(row.tag) != "row":
            continue
        cells, maxc = {}, 0
        for c in row:
            if loc(c.tag) != "c":
                continue
            ci = colnum(c.get("r")) or (len(cells) + 1)
            t = c.get("t")
            raw = None
            for child in c:
                lt = loc(child.tag)
                if lt == "v":
                    raw = child.text
                elif lt == "is":
                    raw = "".join(tt.text or "" for tt in child.iter()
                                 if loc(tt.tag) == "t")
            val = ""
            if raw is not None:
                if t == "s":
                    try:
                        val = shared[int(raw)]
                    except (ValueError, IndexError):
                        val = ""
                else:
                    val = raw
            cells[ci] = val
            maxc = max(maxc, ci)
        yield [str(cells.get(i, "")).strip() for i in range(1, maxc + 1)]


def _config_records(path):
    """Yield one token/cell list per row from tsv/txt/csv/xlsx."""
    low = path.lower()
    if low.endswith((".xlsx", ".xlsm")):
        yield from _read_xlsx_rows(path)
    elif low.endswith(".csv"):
        import csv
        with open(path, newline="") as fh:
            for row in csv.reader(fh):
                yield [c.strip() for c in row]
    else:                                         # tsv / txt: whitespace tokens
        with open(path) as fh:
            for line in fh:
                yield line.rstrip("\n").rstrip("\r").split()


def read_config(path):
    """Read targets from tsv/txt/csv/xlsx. Each row: fastq, then taxids
    (in one comma/space cell or spread across cells), optional trailing
    blast-tsv path. '#' comments and a header row are skipped."""
    jobs = []
    for cells in _config_records(path):
        cells = [c for c in cells if c != ""]
        if not cells or cells[0].startswith("#"):
            continue
        if cells[0].lower() in ("fastq", "file", "filename", "path"):
            continue                              # header row
        fq = cells[0]
        rest = cells[1:]
        blast = None
        if rest and rest[-1].endswith(".tsv"):
            blast = rest[-1]
            rest = rest[:-1]
        taxids = [t for cell in rest for t in re.split(r"[,\s]+", cell) if t]
        jobs.append((fq, taxids, blast))
    return jobs


def _run_folder(fastq):
    """The input's parent-directory name (e.g. MU-CASPER-2026-07-14-a), or ''
    for a bare filename. This is the per-run subfolder the shell writes into."""
    run = os.path.basename(os.path.dirname(fastq))
    return "" if run in (".", "") else run


def default_blast_path(fastq):
    """BLAST TSV path the shell wrote: <run-folder>/<stem>.blastn.nt.tsv, so
    same-named files in different run folders don't collide."""
    b = os.path.basename(fastq)
    for suf in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if b.endswith(suf):
            b = b[:-len(suf)]
            break
    name = b + ".blastn.nt.tsv"
    run = _run_folder(fastq)
    return os.path.join(run, name) if run else name


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", help="targets file (tsv/txt/csv/xlsx): "
                    "<fastq> <taxids> [blast_tsv] per row")
    ap.add_argument("--emit-tsv", metavar="OUT", default=None,
                    help="convert --config to a plain TSV and exit (nothing else)")
    ap.add_argument("--fastq", help="single-file mode: the FASTQ that was BLASTed")
    ap.add_argument("--blast", help="single-file mode: the BLAST TSV")
    ap.add_argument("--match", choices=["name", "taxid", "lineage"], default="name",
                    help="name: species-name substring (default). taxid: exact "
                         "taxid. lineage: hit is the claimed species or a strain "
                         "within it (needs ete3; most correct for coronaviruses)")
    ap.add_argument("--taxdump", default=None,
                    help="path to nodes.dmp (or a dir containing it) for --match "
                         "lineage; if omitted, NCBI's taxdump is downloaded+cached once")
    ap.add_argument("--min-pident", type=float, default=0.0,
                    help="best hit must reach this %% identity to count as confirming")
    ap.add_argument("--min-qcov", type=float, default=0.0,
                    help="best hit must reach this %% query coverage to count")
    ap.add_argument("--only-taxid", nargs="*", default=None,
                    help="single-file mode: restrict to these claimed taxids")
    ap.add_argument("--only-species", nargs="*", default=None,
                    help="single-file mode: restrict to these claimed species")
    ap.add_argument("--blast-dir", default=None,
                    help="root dir holding the run subfolders of *.blastn.nt.tsv "
                         "(match the shell's OUTDIR; default: current dir)")
    ap.add_argument("--sample-contains", nargs="+", default=None,
                    metavar="TOKEN",
                    help="only analyze reads whose header contains ANY of these "
                         "tokens, e.g. --sample-contains Exp3 Exp4. If omitted, "
                         "you're prompted for token(s); blank = all samples.")
    ap.add_argument("--all-samples", action="store_true",
                    help="analyze every sample (skip the prompt, no filtering)")
    ap.add_argument("--out", default="confirmation_summary.csv")
    ap.add_argument("--per-read-out", default="confirmation_per_read.csv")
    args = ap.parse_args()

    # Convert any supported config to a plain TSV and stop (for the shell step).
    if args.emit_tsv:
        if not args.config:
            ap.error("--emit-tsv requires --config")
        jobs = read_config(args.config)
        with open(args.emit_tsv, "w") as out:
            for fq, taxids, blast in jobs:
                cols = [fq, ",".join(taxids)] + ([blast] if blast else [])
                out.write("\t".join(cols) + "\n")
        print(f"Wrote {args.emit_tsv} ({len(jobs)} rows) from {args.config}")
        return

    # Decide the sample filter:
    #   --all-samples            -> no filter (every sample)
    #   --sample-contains A B    -> keep reads matching A or B
    #   neither given            -> ask; a blank answer means all samples
    if args.all_samples:
        sample = []
    elif args.sample_contains is not None:
        sample = args.sample_contains
    elif sys.stdin.isatty():
        resp = input("Sample token(s) to keep, space-separated "
                     "(press Enter for ALL samples): ").strip()
        sample = resp.split() if resp else []
    else:
        sample = []          # non-interactive and nothing specified -> all
    if not sample:
        sys.stderr.write("No sample filter -- analyzing ALL samples.\n")
    global TAXDUMP_HINT
    TAXDUMP_HINT = args.taxdump
    frames = []
    if args.config:
        for fq, taxids, blast in read_config(args.config):
            b = blast or default_blast_path(fq)   # already includes run subfolder
            if args.blast_dir:
                b = os.path.join(args.blast_dir, b)   # prepend root, keep subfolder
            pr = analyze_file(fq, b, taxids, None,
                              args.match, args.min_pident, args.min_qcov, sample)
            run = _run_folder(fq)
            pr["file"] = (f"{run}/{os.path.basename(fq)}" if run
                          else os.path.basename(fq))
            frames.append(pr)
    else:
        if not (args.fastq and args.blast):
            ap.error("provide --config, or both --fastq and --blast")
        pr = analyze_file(args.fastq, args.blast, args.only_taxid, args.only_species,
                          args.match, args.min_pident, args.min_qcov, sample)
        pr["file"] = os.path.basename(args.fastq)
        frames.append(pr)

    per_read = pd.concat(frames, ignore_index=True)
    per_read.to_csv(args.per_read_out, index=False)

    summary = summarize(per_read)
    summary.to_csv(args.out, index=False)
    with pd.option_context("display.max_colwidth", 40, "display.width", 200):
        print(summary.to_string(index=False))
    print(f"\nWrote {args.out} and {args.per_read_out}")


if __name__ == "__main__":
    main()
