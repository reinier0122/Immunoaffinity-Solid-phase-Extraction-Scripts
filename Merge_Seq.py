Sequence-Counts/Merge_Seq.py

#!/usr/bin/env python3
"""
Quality Filter and Adapter Trimming Script using BBDuk
Removes adapters, poly-G tails, and low-quality sequences
"""

import os
import subprocess
import glob
import time
from pathlib import Path

# ============ CONFIGURATION ============
# Path to BBDuk
BBDUK_PATH = ""

# Input directory containing FASTA/FASTQ files
INPUT_FOLDER = ""

# Output folder name (will be created inside INPUT_FOLDER)
OUTPUT_FOLDER_NAME = "Cleaned_Sequences"

# BBDuk parameters
BBDUK_PARAMS = {
    "ref": "adapters",       # Use built-in adapter sequences
    "ktrim": "r",            # Trim adapters from right (3') end
    "k": 23,                 # Kmer length for adapter matching
    "mink": 11,              # Minimum kmer length at read ends
    "hdist": 1,              # Hamming distance (allow 1 mismatch)
    "qtrim": "rl",           # Quality trim both ends (right and left)
    "trimq": 20,             # Phred quality score threshold
    "minlen": 50,            # Minimum read length after trimming
    "entropy": 0.5,          # Filter low-complexity sequences
    "entropywindow": 50,     # Window size for entropy calculation
    "entropyk": 5,           # Kmer size for entropy
    "trimpolya": 8,          # Trim poly-A tails (8+ A's)
    "trimpolyg": 8,          # Trim poly-G tails (8+ G's) - KEY FOR YOUR DATA
    "trimpolyc": 8,          # Trim poly-C tails
    "threads": 4,            # Number of threads
}

# ============ SCRIPT ============

def get_file_size(filepath):
    """Get human-readable file size"""
    if not os.path.exists(filepath):
        return "0 B"
    size = os.path.getsize(filepath)
    for unit in ['B', 'KB', 'MB', 'GB']:
        if size < 1024.0:
            return f"{size:.2f} {unit}"
        size /= 1024.0
    return f"{size:.2f} TB"

def count_sequences(filepath):
    """Count number of sequences in FASTA/FASTQ file"""
    if not os.path.exists(filepath):
        return 0
    
    try:
        if filepath.endswith('.gz'):
            import gzip
            opener = gzip.open
            mode = 'rt'
        else:
            opener = open
            mode = 'r'
        
        with opener(filepath, mode) as f:
            # FASTA files have > at start of header
            # FASTQ files have @ at start of header
            count = sum(1 for line in f if line.startswith('>') or line.startswith('@'))
        return count
    except:
        return 0

def quality_filter(input_file, output_file, stats_file):
    """Quality filter and trim adapters from sequences using BBDuk"""
    
    # Build BBDuk command
    cmd = [
        BBDUK_PATH,
        f"in={input_file}",
        f"out={output_file}",
        f"stats={stats_file}",
    ]
    
    # Add parameters
    for param, value in BBDUK_PARAMS.items():
        cmd.append(f"{param}={value}")
    
    try:
        # Run BBDuk and show output
        result = subprocess.run(cmd, check=True)
        return True
    except subprocess.CalledProcessError as e:
        print(f"  ✗ Error: BBDuk failed with exit code {e.returncode}")
        return False

def main():
    """Main execution function"""
    
    print("Quality Filtering and Adapter Trimming Script")
    print("="*60)
    
    # Verify BBDuk exists
    if not os.path.exists(BBDUK_PATH):
        print(f"Error: BBDuk not found at {BBDUK_PATH}")
        return
    
    # Verify input folder exists
    if not os.path.exists(INPUT_FOLDER):
        print(f"Error: Input folder not found at {INPUT_FOLDER}")
        return
    
    # Create output folder
    output_folder = os.path.join(INPUT_FOLDER, OUTPUT_FOLDER_NAME)
    os.makedirs(output_folder, exist_ok=True)
    
    # Create stats folder
    stats_folder = os.path.join(output_folder, "QC_Stats")
    os.makedirs(stats_folder, exist_ok=True)
    
    print(f"Input folder: {INPUT_FOLDER}")
    print(f"Output folder: {output_folder}")
    print(f"Stats folder: {stats_folder}\n")
    
    print("Quality filtering parameters:")
    print(f"  - Trim adapters: Yes (Illumina)")
    print(f"  - Quality threshold: Q{BBDUK_PARAMS['trimq']}")
    print(f"  - Minimum length: {BBDUK_PARAMS['minlen']} bp")
    print(f"  - Trim poly-G tails: Yes (≥{BBDUK_PARAMS['trimpolyg']} G's)")
    print(f"  - Low complexity filter: Yes (entropy={BBDUK_PARAMS['entropy']})")
    print()
    
    # Find all sequence files (FASTA, FASTQ, compressed)
    patterns = [
        os.path.join(INPUT_FOLDER, "*.fasta"),
        os.path.join(INPUT_FOLDER, "*.fasta.gz"),
        os.path.join(INPUT_FOLDER, "*.fastq"),
        os.path.join(INPUT_FOLDER, "*.fastq.gz"),
        os.path.join(INPUT_FOLDER, "*.fa"),
        os.path.join(INPUT_FOLDER, "*.fa.gz"),
        os.path.join(INPUT_FOLDER, "*.fq"),
        os.path.join(INPUT_FOLDER, "*.fq.gz"),
    ]
    
    input_files = []
    for pattern in patterns:
        input_files.extend(glob.glob(pattern))
    
    # Remove any files from subdirectories
    input_files = [f for f in input_files if os.path.dirname(f) == INPUT_FOLDER]
    input_files = sorted(input_files)
    
    if not input_files:
        print("No sequence files found in the input directory!")
        print("Looking for: *.fasta, *.fastq, *.fa, *.fq (and .gz versions)")
        return
    
    print(f"Found {len(input_files)} file(s) to process\n")
    
    # Process each file
    success_count = 0
    failed_files = []
    total_stats = {
        'input_seqs': 0,
        'output_seqs': 0,
        'input_bases': 0,
        'output_bases': 0
    }
    
    start_time = time.time()
    
    for i, input_file in enumerate(input_files, 1):
        basename = os.path.basename(input_file)
        
        # Generate output filename (keep same extension)
        output_file = os.path.join(output_folder, basename)
        
        # Generate stats filename
        stats_name = basename.replace('.gz', '').replace('.fasta', '').replace('.fastq', '').replace('.fa', '').replace('.fq', '') + '_stats.txt'
        stats_file = os.path.join(stats_folder, stats_name)
        
        # Process
        file_start_time = time.time()
        print(f"[{i}/{len(input_files)}] Processing: {basename}")
        print(f"  Input size: {get_file_size(input_file)}")
        
        input_seqs = count_sequences(input_file)
        if input_seqs > 0:
            print(f"  Input sequences: {input_seqs:,}")
        
        print(f"  Running BBDuk...")
        print("-" * 60)
        
        if quality_filter(input_file, output_file, stats_file):
            success_count += 1
            file_time = time.time() - file_start_time
            
            print("-" * 60)
            if os.path.exists(output_file):
                output_seqs = count_sequences(output_file)
                retained_pct = (output_seqs / input_seqs * 100) if input_seqs > 0 else 0
                
                print(f"  ✓ Completed in {file_time:.1f}s")
                print(f"  Output size: {get_file_size(output_file)}")
                print(f"  Output sequences: {output_seqs:,} ({retained_pct:.1f}% retained)")
                print(f"  Stats saved to: {stats_name}")
                
                # Update totals
                total_stats['input_seqs'] += input_seqs
                total_stats['output_seqs'] += output_seqs
            else:
                print(f"  ✗ Output file not created")
                failed_files.append(basename)
        else:
            print(f"  ✗ Processing failed")
            failed_files.append(basename)
        
        print()
    
    # Summary
    total_time = time.time() - start_time
    print("="*60)
    print("Processing complete!")
    print(f"Successfully processed: {success_count}/{len(input_files)}")
    print(f"Total time: {total_time:.1f} seconds")
    
    if total_stats['input_seqs'] > 0:
        overall_retention = (total_stats['output_seqs'] / total_stats['input_seqs'] * 100)
        print(f"\nOverall statistics:")
        print(f"  Input sequences: {total_stats['input_seqs']:,}")
        print(f"  Output sequences: {total_stats['output_seqs']:,}")
        print(f"  Retention rate: {overall_retention:.1f}%")
        print(f"  Sequences removed: {total_stats['input_seqs'] - total_stats['output_seqs']:,}")
    
    if failed_files:
        print(f"\nFailed files:")
        for filename in failed_files:
            print(f"  - {filename}")
    
    print(f"\nOutput locations:")
    print(f"  Cleaned sequences: {output_folder}")
    print(f"  QC statistics: {stats_folder}")
    print("="*60)

if __name__ == "__main__":
    main()
