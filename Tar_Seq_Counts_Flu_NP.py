Sequence-Counts/Tar_Seq_Counts_Flu_NP.py

import os
import re

def count_sequences_in_fasta(file_path, targets):
    """
    Count occurrences of multiple target sequences in a FASTA file.
    For each target found, sum the counts from the FASTA headers (number after dash).
    Also tracks reads that contain neither target sequence.
    
    Args:
        file_path: Path to the FASTA file
        targets: Dictionary mapping target names to their sequences
    
    Returns:
        Dictionary mapping target names to their summed counts, including "Neither"
    """
    counts = {name: 0 for name in targets.keys()}
    counts["Neither"] = 0
    
    with open(file_path, 'r') as file:
        current_header = None
        current_count = 0
        
        for line in file:
            line = line.strip()
            
            if line.startswith('>'):
                # Parse header to get the count after the dash
                # Format: >ID-COUNT
                current_header = line
                try:
                    count_str = line.split('-')[1]
                    current_count = int(count_str)
                except (IndexError, ValueError):
                    current_count = 0
            else:
                # Process sequence line
                found_target = False
                for name, sequence in targets.items():
                    if sequence in line:
                        counts[name] += current_count
                        found_target = True
                
                # If no target was found in this read, add to "Neither"
                if not found_target:
                    counts["Neither"] += current_count
    
    return counts

def natural_sort_key(filename):
    """
    Generate a sort key for natural sorting (numerical order).
    Converts numeric parts to integers for proper sorting.
    """
    return [int(text) if text.isdigit() else text.lower() 
            for text in re.split('([0-9]+)', filename)]

def count_sequences_in_folder(folder_path, targets):
    """
    Count occurrences of target sequences across all FASTA files in a folder.
    Sums the counts from headers where target sequences are found.
    Reports results for each file separately.
    
    Args:
        folder_path: Path to folder containing .fasta or .fa files
        targets: Dictionary mapping target names to their sequences
    """
    # Get list of FASTA files and sort them numerically
    fasta_files = [f for f in os.listdir(folder_path) 
                   if f.endswith((".fasta", ".fa", ".derep.cut"))
                   and not f.startswith(".")]
    fasta_files.sort(key=natural_sort_key)
    
    for filename in fasta_files:
        file_path = os.path.join(folder_path, filename)
        print("=" * 60)
        print(f"File: {filename}")
        print("=" * 60)
        
        file_counts = count_sequences_in_fasta(file_path, targets)
        
        # Print results for this file
        for name, count in file_counts.items():
            print(f"  {name}: {count} (sum of read counts)")
        print()

if __name__ == "__main__":
    # Configuration: easily add, remove, or modify targets
    targets = {
        "H5N1 Virus": "GCATTCGAGGATCTGAGAGTGTCAAGCTTCATCAGAGGAACAAGAGTAGTCCCAAGAGGACAACTG",
        "H1N1 Virus": "CGCATTTGAAGATCTAAGAGTATTAAGCTTCATCAAAGGGACGAAGGTGCTCCCAAGAGGGAAGCTT"
    }
    
    folder_path = r""
    
    # Display configuration
    print(f"Folder path: {folder_path}\n")
    print("Target sequences:")
    for name, sequence in targets.items():
        print(f"  {name}: {sequence}")
    print("\n" + "=" * 60 + "\n")
    
    # Run the analysis
    count_sequences_in_folder(folder_path, targets)
