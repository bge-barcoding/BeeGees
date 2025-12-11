#!/usr/bin/env python3
import sys
import re
import argparse
from collections import defaultdict
from datetime import datetime

def analyse_fasta(fasta_file, log_file=None, tsv_file=None):
    # Record start time
    start_time = datetime.now()
    
    # Initialise counters
    total_headers = 0
    empty_sequences = 0
    non_empty_sequences = 0
    
    # Counters for different sequence types
    concat_count = 0         # No additional parameters
    merge_only_count = 0     # Has "merge" but not "fcleaner_merge"
    fcleaner_count = 0       # Has "fcleaner" but not "fcleaner_merge"
    fcleaner_merge_count = 0 # Has "fcleaner_merge"
    
    # Create counters for the 24 parameter types
    parameter_counts = defaultdict(int)
    
    # Define the parameter patterns
    r_values = ['r_1', 'r_1.3', 'r_1.5']
    s_values = ['s_50', 's_100']
    
    # Initialise variables to track current sequence
    current_header = None
    current_sequence = ""
    
    # Regular expression to match valid sequence characters
    # Including: uppercase and lowercase nucleotides, N, -, and ~
    valid_seq_chars = re.compile(r'^[GCATNgcatn\-~]+$')
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            
            # Skip empty lines
            if not line:
                continue
                
            # If line starts with >, it's a header
            if line.startswith('>'):
                total_headers += 1
                
                # Process previous sequence if exists
                if current_header:
                    if current_sequence:
                        non_empty_sequences += 1
                        
                        # Check for different header types - ensure mutual exclusivity
                        has_fcleaner_merge = "fcleaner_merge" in current_header
                        has_fcleaner = "fcleaner" in current_header and not has_fcleaner_merge
                        has_merge = "merge" in current_header and not has_fcleaner_merge and not has_fcleaner
                        
                        # Update category counts
                        if has_fcleaner_merge:
                            fcleaner_merge_count += 1
                        elif has_fcleaner:
                            fcleaner_count += 1
                        elif has_merge:
                            merge_only_count += 1
                        else:
                            concat_count += 1
                        
                        # Check for the 24 parameter types
                        for r_val in r_values:
                            for s_val in s_values:
                                param_pattern = f"{r_val}_{s_val}"
                                
                                if param_pattern in current_header:
                                    # Plain parameter
                                    if not has_fcleaner and not has_fcleaner_merge and not has_merge:
                                        parameter_counts[f"{param_pattern}"] += 1
                                    # With merge (but not fcleaner_merge or fcleaner)
                                    elif has_merge:
                                        parameter_counts[f"{param_pattern} (+ \"merge\")"] += 1
                                    # With fcleaner_merge
                                    elif has_fcleaner_merge:
                                        parameter_counts[f"{param_pattern} (+ \"fcleaner_merge\")"] += 1
                                    # With fcleaner (but not fcleaner_merge)
                                    elif has_fcleaner:
                                        parameter_counts[f"{param_pattern} (+ \"fcleaner\")"] += 1
                    else:
                        empty_sequences += 1
                
                # Start new sequence
                current_header = line
                current_sequence = ""
            # Otherwise it's sequence data - verify it contains valid characters
            elif valid_seq_chars.match(line):
                current_sequence += line
            # If it contains invalid characters, log for debugging
            else:
                print(f"Warning: Line with invalid characters found after header {current_header}: {line[:50]}...")
            
        # Process the last sequence
        if current_header:
            if current_sequence:
                non_empty_sequences += 1
                
                # Check for different header types - ensure mutual exclusivity
                has_fcleaner_merge = "fcleaner_merge" in current_header
                has_fcleaner = "fcleaner" in current_header and not has_fcleaner_merge
                has_merge = "merge" in current_header and not has_fcleaner_merge and not has_fcleaner
                
                # Update category counts
                if has_fcleaner_merge:
                    fcleaner_merge_count += 1
                elif has_fcleaner:
                    fcleaner_count += 1
                elif has_merge:
                    merge_only_count += 1
                else:
                    concat_count += 1
                
                # Check for the 24 parameter types for the last sequence
                for r_val in r_values:
                    for s_val in s_values:
                        param_pattern = f"{r_val}_{s_val}"
                        
                        if param_pattern in current_header:
                            # Plain parameter
                            if not has_fcleaner and not has_fcleaner_merge and not has_merge:
                                parameter_counts[f"{param_pattern}"] += 1
                            # With merge (but not fcleaner_merge or fcleaner)
                            elif has_merge:
                                parameter_counts[f"{param_pattern} (+ \"merge\")"] += 1
                            # With fcleaner_merge
                            elif has_fcleaner_merge:
                                parameter_counts[f"{param_pattern} (+ \"fcleaner_merge\")"] += 1
                            # With fcleaner (but not fcleaner_merge)
                            elif has_fcleaner:
                                parameter_counts[f"{param_pattern} (+ \"fcleaner\")"] += 1
            else:
                empty_sequences += 1
    
    # Record end time
    end_time = datetime.now()
    
    # Verify our category counts add up to the total
    category_sum = concat_count + merge_only_count + fcleaner_count + fcleaner_merge_count
    
    # Prepare results as a list of strings
    results = []
    results.append(f"Total sequences: {total_headers}")
    results.append(f"Empty sequences: {empty_sequences}")
    results.append(f"Non-empty sequences: {non_empty_sequences}")
    results.append(f"")
    results.append(f"Non-empty sequences breakdown by mode type:")
    results.append(f"  'concat' mode sequences: {concat_count}")
    results.append(f"  'merge' mode sequences: {merge_only_count}")
    results.append(f"  'fcleaner' sequences: {fcleaner_count}")
    results.append(f"  'fcleaner_merge' sequences: {fcleaner_merge_count}")
    
    results.append("\nParameter Type Counts:")
    # Sort parameter types to ensure consistent output order
    for param_type in sorted(parameter_counts.keys()):
        results.append(f"{param_type}: {parameter_counts[param_type]}")
    
    # Print results to console
    for line in results:
        print(line)
    
    # Write results to log file if specified (with timestamps)
    if log_file:
        with open(log_file, 'w') as f:
            # Write start timestamp
            f.write("=" * 70 + '\n')
            f.write(f"Script: analyse_fasta.py\n")
            f.write(f"Start time: {start_time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Input file: {fasta_file}\n")
            f.write("=" * 70 + '\n\n')
            
            # Write analysis results
            for line in results:
                f.write(line + '\n')
            
            # Write end timestamp
            f.write('\n')
            f.write("=" * 70 + '\n')
            f.write(f"End time: {end_time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Duration: {end_time - start_time}\n")
            f.write("=" * 70 + '\n')
        print(f"\nResults have been written to: {log_file}")
    
    # Write results to TSV file if specified
    if tsv_file:
        with open(tsv_file, 'w') as f:
            # Write summary statistics
            f.write("Metric\tCount\n")
            f.write(f"Total sequences\t{total_headers}\n")
            f.write(f"Empty sequences\t{empty_sequences}\n")
            f.write(f"Non-empty sequences\t{non_empty_sequences}\n")
            f.write("\n")
            
            # Write mode type breakdown
            f.write("Mode Type\tCount\n")
            f.write(f"concat mode sequences\t{concat_count}\n")
            f.write(f"merge mode sequences\t{merge_only_count}\n")
            f.write(f"fcleaner sequences\t{fcleaner_count}\n")
            f.write(f"fcleaner_merge sequences\t{fcleaner_merge_count}\n")
            f.write("\n")
            
            # Write parameter type counts
            f.write("Parameter Type\tCount\n")
            for param_type in sorted(parameter_counts.keys()):
                f.write(f"{param_type}\t{parameter_counts[param_type]}\n")
        
        print(f"Results have been written to TSV: {tsv_file}")

def main():
    # Set up command line argument parsing
    parser = argparse.ArgumentParser(description='analyse FASTA file for sequence parameters')
    parser.add_argument('--input', required=True, help='Input FASTA file path')
    parser.add_argument('--log-file', help='Output log file path (optional)')
    parser.add_argument('--tsv', help='Output TSV file path (optional)')
    
    args = parser.parse_args()
    
    # Run the analysis
    analyse_fasta(args.input, args.log_file, args.tsv)

if __name__ == "__main__":
    main()
