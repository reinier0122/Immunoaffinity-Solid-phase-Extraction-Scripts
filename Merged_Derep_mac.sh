#!/bin/bash
# =============================================================================
# Paired-end Illumina amplicon preprocessing pipeline (Alpha-Gal)
#
# Per sample:  merge R1/R2 (bbmerge) -> dereplicate (derep.py)
#              -> trim 5' primer (cutadapt) -> trim 3' primer (cutadapt)
#              -> final FASTA in Final_Sequences/
#
# Accepts BOTH uncompressed (.fastq) and gzipped (.fastq.gz) input in the SAME
# run. bbmerge reads .gz natively, so no manual decompression is required.
#
# Safe to re-run: any sample whose final .derep.cut.fa already exists is skipped.
# =============================================================================

# ------------------------------ Configuration --------------------------------
# Paths
BASE="/Users/reiniersuarez/Library/CloudStorage/Dropbox/Alpha-Gal_Project/Python_Scripts_aGal"
DATA="/Volumes/Rey's ExHDD/Alpha-Gal_Data/Illumina_Sequences/2026_Data/11_November_2026"

# Analysis parameters — the knobs worth reviewing per run / amplicon.
# NOTE: these values are UNCHANGED from your original script; they're just
#       surfaced here so they're easy to find and tweak. No science was altered.
DEREP_MIN=10                          # derep.py: min times a seq must be seen to survive
BBMERGE_MEM=4000m                     # JVM heap for bbmerge (-Xmx / -Xms)
BBMERGE_THREADS=2                     # bbmerge threads
PRIMER_5="acactctttccctacacgacgctcttccgatct"     # forward primer, trimmed at 5' forward orientation (anchored with ^)
PRIMER_3="gtgactggagttcagacgtgtgctcttccgatct"    # reverse primer, trimmed at 3' forward orientation (anchored with $)
CUT_ERR="0.3"                         # cutadapt error tolerance (30%) — LOOSE, see note
# ^ 0.3 forgives primer degeneracy but will also trim on weak matches. If you
#   see over-trimming, lower it (e.g. 0.1–0.15). Left as-is intentionally.
# -----------------------------------------------------------------------------

# Counters
processed=0
skipped=0
current=0

# --- Dependency check --------------------------------------------------------
# Replaces the old `export PATH="...Python/3.x/bin"` line, which had a literal
# "3.x" in it and so added a directory that never exists (a silent no-op). We
# instead verify the tools we actually call are reachable, and fail early with
# a clear message. Install cutadapt via pip/conda and java via Adoptium/Homebrew
# however you already do — this just checks they're on PATH.
for tool in java python3 cutadapt; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "ERROR: required tool '$tool' not found on PATH"
        exit 1
    fi
done

# Work in the data directory so everything below uses relative paths
cd "$DATA" || { echo "ERROR: Could not cd to $DATA"; exit 1; }

# --- Gather input ------------------------------------------------------------
# nullglob makes a non-matching pattern expand to nothing (not the literal
# string), so the list is clean whether you have .fastq, .fastq.gz, or a mix.
shopt -s nullglob
r1_files=(*_R1_001.fastq *_R1_001.fastq.gz)

total=${#r1_files[@]}
echo "Found $total samples to process"
echo ""

if [[ "$total" -eq 0 ]]; then
    echo "No *_R1_001.fastq or *_R1_001.fastq.gz files here — nothing to do."
    echo "(Check you're pointed at the right run folder: $DATA)"
    exit 0
fi

# Create output directories up front
mkdir -p Final_Sequences Intermediate_Files Merged_Files Derep_Files CutInfo_Files
echo "Created all output directories"
echo ""

# --- Per-sample loop ---------------------------------------------------------
for file in "${r1_files[@]}"; do
    current=$((current + 1))

    # Strip the correct suffix depending on whether this file is gzipped, and
    # remember the extension so the R2 name and raw-file moves match it.
    if [[ "$file" == *.fastq.gz ]]; then
        Sampid="${file%_R1_001.fastq.gz}"
        ext=".fastq.gz"
    else
        Sampid="${file%_R1_001.fastq}"
        ext=".fastq"
    fi
    R2="${Sampid}_R2_001${ext}"

    # Resume support: skip samples whose FINAL output already exists.
    # (Keys only on the final file. A sample that died mid-run left its raw R1
    #  in place and no final output, so it re-runs from the top and overwrites
    #  its own stale intermediates by name — which is what you want.)
    if [[ -f "./Final_Sequences/${Sampid}.derep.cut.fa" ]]; then
        skipped=$((skipped + 1))
        echo "[$current/$total] Skipping $Sampid — already completed"
        continue
    fi

    # Guard against an orphaned R1 with no matching R2.
    if [[ ! -f "$R2" ]]; then
        echo "[$current/$total] WARNING: no R2 for $Sampid (expected $R2) — skipping"
        continue
    fi

    echo "================================================"
    echo "[$current/$total] Processing file: $file"
    echo "Sample ID: $Sampid"

    # [Step 1/4] Merge R1/R2 read pairs. bbmerge reads .gz input natively.
    echo "[Step 1/4] Running bbmerge..."
    start_time=$(date +%s)
    java -ea -Xmx${BBMERGE_MEM} -Xms${BBMERGE_MEM} \
        -cp "$BASE/bbmap/current/" jgi.BBMerge \
        threads=${BBMERGE_THREADS} \
        in1="$file" \
        in2="$R2" \
        out="$Sampid.merge.fq" \
        outu="$Sampid.unmerge.fq" > "$Sampid.mergestats.txt" 2>&1 
    echo "bbmerge completed in $(( $(date +%s) - start_time ))s"

    # [Step 2/4] Dereplicate the merged reads (your custom script).
    echo "[Step 2/4] Running dereplication..."
    start_time=$(date +%s)
    python3 "$BASE/derep.py" \
        "$Sampid.merge.fq" \
        "$Sampid.derep.fa" \
        "$DEREP_MIN" > "${Sampid}_derepinfo.txt" 2>&1
    echo "Dereplication completed in $(( $(date +%s) - start_time ))s"

    # [Step 3/4] Trim the forward primer from the 5' end (^ = anchored to start).
    echo "[Step 3/4] Running first cutadapt..."
    start_time=$(date +%s)
    cutadapt -g "^${PRIMER_5}" -e "$CUT_ERR" \
        -o "$Sampid.derep.cut1.fa" \
        "$Sampid.derep.fa" >> "$Sampid.CutInfo.txt" 2>&1
    echo "First cutadapt completed in $(( $(date +%s) - start_time ))s"

    # [Step 4/4] Trim the reverse primer from the 3' end ($ = anchored to end).
    echo "[Step 4/4] Running second cutadapt..."
    start_time=$(date +%s)
    cutadapt -a "${PRIMER_3}\$" -e "$CUT_ERR" \
        -o "$Sampid.derep.cut.fa" \
        "$Sampid.derep.cut1.fa" >> "$Sampid.CutInfo.txt" 2>&1
    echo "Second cutadapt completed in $(( $(date +%s) - start_time ))s"

    # --- File the outputs into their subfolders ------------------------------
    if [[ -f "${Sampid}.derep.cut.fa" ]]; then
        mv "${Sampid}.derep.cut.fa" ./Final_Sequences/
        echo "Moved ${Sampid}.derep.cut.fa to Final_Sequences/"
    else
        echo "WARNING: ${Sampid}.derep.cut.fa not found — check ${Sampid}.mergestats.txt and CutInfo"
    fi

    # Raw reads (carry whatever extension they came in with)
    mv "$file" ./Intermediate_Files/ 2>/dev/null
    mv "$R2"   ./Intermediate_Files/ 2>/dev/null

    # Merged reads + stats
    mv "${Sampid}.merge.fq"       ./Merged_Files/  2>/dev/null
    mv "${Sampid}.mergestats.txt" ./Merged_Files/  2>/dev/null

    # Derep sequences + info
    mv "${Sampid}.derep.fa"      ./Derep_Files/    2>/dev/null
    mv "${Sampid}_derepinfo.txt" ./Derep_Files/    2>/dev/null

    # Intermediate cutadapt file + trim log
    mv "${Sampid}.derep.cut1.fa" ./CutInfo_Files/  2>/dev/null
    mv "${Sampid}.CutInfo.txt"   ./CutInfo_Files/  2>/dev/null

    processed=$((processed + 1))
    echo "✓ Sample $Sampid processed successfully"
    echo "Progress: $processed processed, $skipped skipped out of $current seen"
done

echo ""
echo "================================================"
echo "Pipeline completed at $(date)"
echo "Total samples processed: $processed"
echo "Total samples skipped:   $skipped"
echo "Total samples seen:      $current out of $total"
echo ""
echo "Final sequences in:       Final_Sequences/"
echo "Raw FASTQ in:             Intermediate_Files/"
echo "Merged reads + stats in:  Merged_Files/"
echo "Derep sequences + info:   Derep_Files/"
echo "Cutadapt files + info:    CutInfo_Files/"