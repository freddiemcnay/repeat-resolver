#!/usr/bin/env python3
"""
Repeat Structure Classifier
Classifies the structure of reads containing repeat expansions by comparing
sliding windows to known unit sequences using minimap2.

Author: Freddie
Date: 2026-02-04
"""

import argparse
import subprocess
import tempfile
import os
from pathlib import Path
from collections import defaultdict


def parse_fasta(fasta_file):
    """
    Parse FASTA file and return dictionary of {header: sequence}
    """
    sequences = {}
    current_header = None
    current_seq = []
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header:
                    sequences[current_header] = ''.join(current_seq)
                current_header = line[1:].split()[0]  # Get first word after >
                current_seq = []
            else:
                current_seq.append(line)
        
        # Don't forget the last sequence
        if current_header:
            sequences[current_header] = ''.join(current_seq)
    
    return sequences


def sanitize_name(name):
    """
    Sanitize a name to be safe for use in filenames
    Replace problematic characters with underscores
    """
    # Replace forward slashes, backslashes, and other problematic characters
    return name.replace('/', '_').replace('\\', '_').replace(' ', '_')


def write_fasta(header, sequence, filename):
    """
    Write a single sequence to a FASTA file
    """
    with open(filename, 'w') as f:
        f.write(f">{header}\n")
        # Write sequence in 80-character lines
        for i in range(0, len(sequence), 80):
            f.write(sequence[i:i+80] + '\n')


def run_minimap2(query_fasta, ref_fasta, preset='map-hifi'):
    """
    Run minimap2 and return PAF output as list of lines
    """
    cmd = [
        'minimap2',
        '-c',  # Output in PAF format
        '-x', preset,
        ref_fasta,
        query_fasta
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        return result.stdout.strip().split('\n') if result.stdout.strip() else []
    except subprocess.CalledProcessError as e:
        print(f"Warning: minimap2 failed with error: {e}")
        return []


def parse_paf_line(paf_line):
    """
    Parse a single PAF line and return relevant information
    Returns: dict with query_name, query_len, query_start, query_end, strand,
             target_name, target_len, target_start, target_end, matches, aln_len, mapq
    """
    if not paf_line:
        return None
    
    fields = paf_line.split('\t')
    if len(fields) < 12:
        return None
    
    return {
        'query_name': fields[0],
        'query_len': int(fields[1]),
        'query_start': int(fields[2]),
        'query_end': int(fields[3]),
        'strand': fields[4],
        'target_name': fields[5],
        'target_len': int(fields[6]),
        'target_start': int(fields[7]),
        'target_end': int(fields[8]),
        'matches': int(fields[9]),
        'aln_len': int(fields[10]),
        'mapq': int(fields[11])
    }


def calculate_identity(matches, aln_len):
    """
    Calculate percent identity from PAF alignment
    """
    if aln_len == 0:
        return 0.0
    return (matches / aln_len) * 100


def is_uppercase_assignment(target_start, target_len, threshold=0.3):
    """
    Determine if assignment should be uppercase based on alignment position
    Returns True (uppercase) if alignment includes first 30% of reference
    """
    # If alignment starts within first 30% of reference, it's uppercase
    return (target_start / target_len) < threshold


def classify_window(window_seq, window_name, unit_seqs, unit_names, 
                   temp_dir, identity_threshold=40.0, uppercase_threshold=0.3,
                   preset='map-hifi'):
    """
    Classify a single window against all unit sequences
    
    Returns:
        - assigned_unit: The assigned unit (uppercase or lowercase) or 'U'
        - scores_dict: Dictionary with alignment info for each unit
        - strand: The strand of the best alignment ('+', '-', or None)
    """
    # Write window to temporary FASTA (sanitize name for filesystem)
    safe_window_name = sanitize_name(window_name)
    window_fasta = os.path.join(temp_dir, f"{safe_window_name}.fa")
    write_fasta(window_name, window_seq, window_fasta)
    
    # Store results for each unit
    scores_dict = {}
    best_identity = 0.0
    best_unit = None
    best_target_start = None
    best_target_len = None
    best_strand = None
    
    # Run minimap2 against each unit separately
    for unit_name in unit_names:
        unit_fasta = os.path.join(temp_dir, f"unit_{unit_name}.fa")
        write_fasta(unit_name, unit_seqs[unit_name], unit_fasta)
        
        paf_output = run_minimap2(window_fasta, unit_fasta, preset)
        
        if paf_output and paf_output[0]:
            # Parse the best alignment (first line)
            aln = parse_paf_line(paf_output[0])
            if aln:
                identity = calculate_identity(aln['matches'], aln['aln_len'])
                scores_dict[unit_name] = {
                    'identity': identity,
                    'target_start': aln['target_start'],
                    'target_end': aln['target_end'],
                    'target_len': aln['target_len'],
                    'strand': aln['strand']
                }
                
                # Track best alignment
                if identity > best_identity:
                    best_identity = identity
                    best_unit = unit_name
                    best_target_start = aln['target_start']
                    best_target_len = aln['target_len']
                    best_strand = aln['strand']
            else:
                scores_dict[unit_name] = {
                    'identity': 0.0,
                    'target_start': 'NA',
                    'target_end': 'NA',
                    'target_len': len(unit_seqs[unit_name]),
                    'strand': 'NA'
                }
        else:
            # No alignment found
            scores_dict[unit_name] = {
                'identity': 0.0,
                'target_start': 'NA',
                'target_end': 'NA',
                'target_len': len(unit_seqs[unit_name]),
                'strand': 'NA'
            }
    
    # Determine final assignment
    if best_identity < identity_threshold or best_unit is None:
        assigned_unit = 'U'
    else:
        # Check if should be uppercase
        if is_uppercase_assignment(best_target_start, best_target_len, uppercase_threshold):
            assigned_unit = best_unit.upper()
        else:
            assigned_unit = best_unit.lower()
    
    return assigned_unit, scores_dict, best_strand


def process_read(read_name, read_seq, unit_seqs, unit_names, temp_dir,
                window_size=1000, identity_threshold=40.0, 
                uppercase_threshold=0.3, preset='map-hifi'):
    """
    Process a single read by sliding window analysis
    
    Returns:
        - structure: String like "A-a-a-B-C-U-D"
        - detailed_scores: List of dictionaries with scores for each window
        - read_strand: Overall strand of the read ('+', '-', or 'mixed')
    """
    structure = []
    detailed_scores = []
    strand_votes = []  # Track strand for each assigned window
    
    read_len = len(read_seq)
    window_num = 0
    
    # Slide window across read
    for start in range(0, read_len, window_size):
        window_num += 1
        end = min(start + window_size, read_len)
        window_seq = read_seq[start:end]
        window_name = f"{read_name}_window{window_num}"
        
        # Classify this window
        assigned_unit, scores, strand = classify_window(
            window_seq, window_name, unit_seqs, unit_names, temp_dir,
            identity_threshold, uppercase_threshold, preset
        )
        
        structure.append(assigned_unit)
        
        # Track strand for non-U assignments
        if assigned_unit != 'U' and strand is not None:
            strand_votes.append(strand)
        
        # Store detailed scores
        score_entry = {
            'read_name': read_name,
            'window_num': window_num,
            'window_start': start,
            'window_end': end,
            'assigned_unit': assigned_unit,
            'strand': strand if strand else 'NA'
        }
        
        # Add scores for each unit
        for unit_name in unit_names:
            score_entry[f'{unit_name}_identity'] = scores[unit_name]['identity']
            score_entry[f'{unit_name}_ref_start'] = scores[unit_name]['target_start']
            score_entry[f'{unit_name}_ref_end'] = scores[unit_name]['target_end']
            score_entry[f'{unit_name}_strand'] = scores[unit_name]['strand']
        
        detailed_scores.append(score_entry)
    
    # Determine overall read strand by majority vote
    if not strand_votes:
        read_strand = 'unknown'
    else:
        plus_count = strand_votes.count('+')
        minus_count = strand_votes.count('-')
        
        if plus_count > minus_count * 2:  # Strongly forward
            read_strand = '+'
        elif minus_count > plus_count * 2:  # Strongly reverse
            read_strand = '-'
        else:
            read_strand = 'mixed'
    
    structure_string = '-'.join(structure)
    return structure_string, detailed_scores, read_strand


def interpret_structure(structure_string, strand='+'):
    """
    Convert structure string like 'b-C-C-C-c-c-c-c-c-c-C-C-C-c-c-c-c-c'
    into simplified format like 'B-C×2'
    
    If strand is '-', the structure is reversed before interpretation to give
    the true 5'->3' structure of the read.
    
    Logic:
    1. If reverse strand: reverse the structure string
    2. Group consecutive windows of same unit
    3. Count tandem copies by looking for lowercase→uppercase transitions
    4. Simplify consecutive Us to single U
    5. Output format: Unit×count (if count > 1)
    """
    # If reverse strand, reverse the structure to get true 5'->3' orientation
    if strand == '-':
        elements = structure_string.split('-')
        elements.reverse()
        structure_string = '-'.join(elements)
    
    elements = structure_string.split('-')
    
    result = []
    i = 0
    
    while i < len(elements):
        current_unit = elements[i].upper()
        
        # Special handling for U (undetermined)
        if current_unit == 'U':
            # Skip all consecutive Us
            j = i
            while j < len(elements) and elements[j].upper() == 'U':
                j += 1
            result.append('U')
            i = j
            continue
        
        # Count tandem copies for this unit
        tandem_count = 1  # Start with 1 (first occurrence)
        j = i
        prev_was_lower = False
        
        while j < len(elements) and elements[j].upper() == current_unit:
            # Count transition from lowercase to uppercase as new tandem
            if elements[j].isupper() and prev_was_lower:
                tandem_count += 1
            
            prev_was_lower = elements[j].islower()
            j += 1
        
        # Format output
        if tandem_count > 1:
            result.append(f"{current_unit}×{tandem_count}")
        else:
            result.append(current_unit)
        
        i = j
    
    return '-'.join(result)


def main():
    parser = argparse.ArgumentParser(
        description='Classify repeat structure in long reads using minimap2',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('-i', '--input', required=True,
                       help='Input FASTA file with unassigned reads')
    parser.add_argument('-u', '--units', required=True,
                       help='FASTA file with unit sequences')
    parser.add_argument('-o', '--output', required=True,
                       help='Output prefix for results files')
    parser.add_argument('-w', '--window-size', type=int, default=1000,
                       help='Window size in base pairs')
    parser.add_argument('-t', '--identity-threshold', type=float, default=40.0,
                       help='Minimum percent identity for assignment (below this = U)')
    parser.add_argument('-c', '--uppercase-threshold', type=float, default=0.3,
                       help='Fraction of unit length for uppercase assignment (0.3 = first 30%%)')
    parser.add_argument('-p', '--preset', default='map-hifi',
                       choices=['map-ont', 'map-hifi', 'map-pb', 'asm20'],
                       help='Minimap2 preset')
    
    args = parser.parse_args()
    
    # Validate inputs
    if not os.path.exists(args.input):
        print(f"Error: Input file {args.input} not found")
        return 1
    if not os.path.exists(args.units):
        print(f"Error: Units file {args.units} not found")
        return 1
    
    # Parse input files
    print("Reading input sequences...")
    reads = parse_fasta(args.input)
    units = parse_fasta(args.units)
    
    print(f"Found {len(reads)} reads to process")
    print(f"Found {len(units)} unit sequences: {', '.join(units.keys())}")
    
    unit_names = sorted(units.keys())
    
    # Create temporary directory for intermediate files
    with tempfile.TemporaryDirectory() as temp_dir:
        # Open output files
        structure_file = f"{args.output}_structures.txt"
        scores_file = f"{args.output}_detailed_scores.txt"
        summary_file = f"{args.output}_structure_summary.txt"
        
        # Dictionary to track structure frequencies
        structure_counts = defaultdict(int)
        strand_distribution = defaultdict(int)
        
        with open(structure_file, 'w') as sf, open(scores_file, 'w') as df:
            # Write headers
            sf.write("read_name\tstrand\tfull_structure\tsimplified_structure\n")
            
            # Header for detailed scores
            header_parts = ['read_name', 'window_num', 'window_start', 'window_end', 'strand']
            for unit_name in unit_names:
                header_parts.extend([
                    f'{unit_name}_identity',
                    f'{unit_name}_ref_start',
                    f'{unit_name}_ref_end',
                    f'{unit_name}_strand'
                ])
            header_parts.append('assigned_unit')
            df.write('\t'.join(header_parts) + '\n')
            
            # Process each read
            for idx, (read_name, read_seq) in enumerate(reads.items(), 1):
                print(f"Processing read {idx}/{len(reads)}: {read_name}")
                
                structure, detailed, read_strand = process_read(
                    read_name, read_seq, units, unit_names, temp_dir,
                    args.window_size, args.identity_threshold,
                    args.uppercase_threshold, args.preset
                )
                
                # Write structure with strand information
                simplified = interpret_structure(structure, read_strand)
                sf.write(f"{read_name}\t{read_strand}\t{structure}\t{simplified}\n")
                
                # Track structure frequency and strand distribution
                structure_counts[simplified] += 1
                strand_distribution[read_strand] += 1
                
                # Write detailed scores
                for score_entry in detailed:
                    row_parts = [
                        score_entry['read_name'],
                        str(score_entry['window_num']),
                        str(score_entry['window_start']),
                        str(score_entry['window_end']),
                        str(score_entry['strand'])
                    ]
                    
                    for unit_name in unit_names:
                        row_parts.extend([
                            f"{score_entry[f'{unit_name}_identity']:.2f}",
                            str(score_entry[f'{unit_name}_ref_start']),
                            str(score_entry[f'{unit_name}_ref_end']),
                            str(score_entry[f'{unit_name}_strand'])
                        ])
                    
                    row_parts.append(score_entry['assigned_unit'])
                    df.write('\t'.join(row_parts) + '\n')
        
        # Generate structure summary table
        total_reads = len(reads)
        with open(summary_file, 'w') as sumf:
            sumf.write("# Structure Summary\n")
            sumf.write(f"# Total reads processed: {total_reads}\n")
            sumf.write(f"# Window size: {args.window_size} bp\n")
            sumf.write(f"# Identity threshold: {args.identity_threshold}%\n")
            sumf.write("#\n")
            
            # Strand distribution summary
            sumf.write("# Strand Distribution:\n")
            for strand in sorted(strand_distribution.keys()):
                count = strand_distribution[strand]
                pct = (count / total_reads) * 100
                sumf.write(f"#   {strand}: {count} reads ({pct:.1f}%)\n")
            sumf.write("#\n")
            
            # Structure frequency table
            sumf.write("simplified_structure\tcount\tpercentage\n")
            
            # Sort by frequency (most common first), then alphabetically for ties
            sorted_structures = sorted(
                structure_counts.items(),
                key=lambda x: (-x[1], x[0])  # -count (descending), then structure (ascending)
            )
            
            for structure, count in sorted_structures:
                percentage = (count / total_reads) * 100
                sumf.write(f"{structure}\t{count}\t{percentage:.2f}\n")
    
    print(f"\nComplete! Results written to:")
    print(f"  Structures: {structure_file}")
    print(f"  Detailed scores: {scores_file}")
    print(f"  Structure summary: {summary_file}")
    print(f"\nStructure Summary:")
    print(f"  Total reads: {total_reads}")
    print(f"  Unique structures: {len(structure_counts)}")
    if sorted_structures:
        top_structure, top_count = sorted_structures[0]
        top_pct = (top_count / total_reads) * 100
        print(f"  Most common: {top_structure} ({top_count} reads, {top_pct:.1f}%)")
    
    return 0


if __name__ == '__main__':
    exit(main())
