#!/bin/bash
#SBATCH --partition=marcjohnson-lab
#SBATCH --account=marcjohnson-lab
#SBATCH --job-name=Blast
#SBATCH --time=5-16:00:00
#SBATCH --export=all
#SBATCH -n 64 # max 128
#SBATCH --mem 245GB # max / node 490GB
#SBATCH --out=Blast-%j.out
#SBATCH --mail-type=begin,end,fail,requeue
#SBATCH --mail-user=rswz6@umsystem.edu

export PATH=$PATH:/cluster/VAST/marcjohnson-lab/jkome/BioinformaticTools/ncbi-blast-2.15.0+/bin
export BLASTDB=/cluster/VAST/marcjohnson-lab/Databases/BLAST/nt

# Python is needed only if CONFIG is .xlsx/.csv (for the auto-conversion below).
# Harmless if it's already a .tsv. Needs pandas+openpyxl for .xlsx:
#   pip install --user pandas openpyxl   (run once on the login node)
module load python/3.13.3 2>/dev/null || true

# ------------------------------------------------------------------
# Edit targets.tsv, not this script. Each row = one FASTQ + the taxids
# to validate in it. Inputs may be .fastq or .fastq.gz (also .fq/.fq.gz).
CONFIG="${CONFIG:-targets.tsv}"
MAXTARGETS=50
# Keep reads whose header (sample name) contains ANY of these tokens, space-separated. 
#This is how samples are tagged (e.g. "REY", "exp3"). Examples:
#   SAMPLE_KEEP="REY"            SAMPLE_KEEP="exp3 exp4"      SAMPLE_KEEP=""(all)
SAMPLE_KEEP="Exp3 Exp4 REY Exp2"
# Root for outputs. A per-run subfolder (named after the input's parent dir,
# e.g. MU-CASPER-2026-07-14-a) is created under this, so files with the same
# name in different run folders don't overwrite each other.
OUTDIR="${OUTDIR:-.}"
outfmt="6 qseqid sseqid pident staxid staxids sskingdom sscinames stitle qcovs evalue bitscore"
# Tab-separated header for the output TSV: the field names from outfmt, minus
# the leading format code "6". Kept in sync automatically if you edit outfmt.
hdr=$(printf '%s' "$outfmt" | cut -d' ' -f2- | tr ' ' '\t')
# ------------------------------------------------------------------

[ -e "$CONFIG" ] || { echo "Config $CONFIG not found."; exit 1; }

# The BLAST loop reads a plain TSV. If CONFIG is .xlsx/.csv, convert it with
# xlsx_to_tsv.py, which uses only the Python standard library (no pandas).
case "${CONFIG,,}" in
    *.tsv|*.txt)
        CFG="$CONFIG" ;;
    *)
        CFG="targets.autogen.tsv"
        echo "Converting $CONFIG -> $CFG"
        python3 "$(dirname "$0")/xlsx_to_tsv.py" "$CONFIG" "$CFG" \
            || { echo "Conversion failed. Convert on a login node, then rerun with the .tsv."; exit 1; }
        ;;
esac

# stream a FASTQ to stdout whether it's plain or gzipped
reader() { case "$1" in *.gz) gzip -dc "$1";; *) cat "$1";; esac; }
# basename with fastq/fq[.gz] stripped
strip_fastq() { local b; b=$(basename "$1"); b=${b%.gz}; b=${b%.fastq}; b=${b%.fq}; printf '%s' "$b"; }

while IFS= read -r line || [ -n "$line" ]; do
    line=${line%$'\r'}                       # strip Windows CR
    line="${line#"${line%%[![:space:]]*}"}"  # left-trim
    [[ -z "$line" || "$line" == \#* ]] && continue

    read -r fq taxids <<< "$line"            # fq = first token, taxids = rest
    taxids=${taxids//,/ }                     # commas -> spaces for awk
    [ -n "$taxids" ] || { echo "No taxids for $fq, skipping"; continue; }
    [ -e "$fq" ] || { echo "File $fq not found, skipping"; continue; }

    base=$(strip_fastq "$fq")
    # Output goes in OUTDIR/<run-folder>/ , where <run-folder> is the input's
    # parent directory name. Bare filenames (no dir) fall back to OUTDIR itself.
    run=$(basename "$(dirname "$fq")")
    if [ "$run" = "." ] || [ -z "$run" ]; then dest="$OUTDIR"; else dest="$OUTDIR/$run"; fi
    mkdir -p "$dest"

    sub="$dest/${base}.target.fastq"
    fasta="$dest/${base}.fasta"
    full="$dest/${base}.blastn.nt.tsv"

    # Keep a read only if (a) its taxid (field 2) is a target for this file AND
    # (b) its header contains ANY SAMPLE_KEEP token as a delimited word.
    reader "$fq" | awk -v keep="$taxids" -v samp="$SAMPLE_KEEP" '
        BEGIN{ n=split(keep,a," "); for(i=1;i<=n;i++) K[a[i]]=1
               m=split(tolower(samp), s, /[ ,]+/); alt=""
               for(i=1;i<=m;i++) if(s[i]!="") alt=(alt==""? s[i] : alt "|" s[i])
               re=(alt=="") ? "" : ("(^|[^a-z0-9])(" alt ")([^a-z0-9]|$)") }
        NR%4==1{ split($0,f," ")
                 wt=(f[2] in K)
                 ws=(re=="") ? 1 : (tolower(f[1]) ~ re)
                 want=(wt && ws) }
        want' > "$sub"

    kept=$(( $(wc -l < "$sub") / 4 ))
    echo "=== $run/$base: kept $kept reads for taxids [$taxids] -> $full ==="
    [ "$kept" -eq 0 ] && { echo "  none matched, skipping"; continue; }

    # FASTQ -> FASTA (blastn needs FASTA; keep full header for the claim).
    awk 'NR%4==1{sub(/^@/,">"); print} NR%4==2{print}' "$sub" > "$fasta"

    blastn -query "$fasta" -db nt -task megablast -num_threads 60 -evalue 1 \
        -max_target_seqs "$MAXTARGETS" \
        -outfmt "$outfmt" -out "${full}.body"
    # Prepend the column-name header, then the BLAST rows.
    { printf '%s\n' "$hdr"; cat "${full}.body"; } > "$full"
    rm -f "${full}.body"
done < "$CFG"

echo "=== Done. Now run: python blast_confirm.py --config $CONFIG --match lineage ==="