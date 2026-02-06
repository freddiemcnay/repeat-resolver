# Repeat Structure Classifier

A bioinformatics tool for analyzing long sequencing reads containing complex repeat expansions. Identifies and quantifies tandem repeat structures by aligning sliding windows against known repeat unit sequences.

## Features

- **Strand-aware analysis**: Correctly interprets reverse complement alignments
- **Tandem repeat detection**: Automatically identifies and counts tandem copies using uppercase/lowercase logic
- **Quality metrics**: Per-window alignment scores for troubleshooting
- **Population summary**: Frequency table showing dominant repeat configurations across your dataset
- **Flexible parameters**: Adjustable window size, identity threshold, and minimap2 presets

## Quick Start

```bash
# Basic usage
./repeat_classifier_final.py \
    -i reads.fasta \
    -u units.fasta \
    -o output_prefix

# With custom parameters
./repeat_classifier_final.py \
    -i reads.fasta \
    -u units.fasta \
    -o output_prefix \
    -w 2000 \
    -t 45.0 \
    -p map-ont
```

## Requirements

- Python 3.6+
- minimap2 (must be in your PATH)

```bash
# Verify minimap2 is available
minimap2 --version

# Make script executable
chmod +x repeat_classifier_final.py
```

## Input Files

### 1. Reads File (-i, --input)
FASTA file with your long reads to classify:

```
>read001
ATCGATCGATCG...
>read002
GCTAGCTAGCTA...
```

**Supported read types**: PacBio HiFi, PacBio CLR, Oxford Nanopore

### 2. Units File (-u, --units)
FASTA file with known repeat unit sequences:

```
>A
ATCGATCGATCGATCG...
>B
GCTAGCTAGCTAGCTA...
>C
TTAATTAATTAATTAA...
```

**Tips**: 
- Use simple single-letter names (A, B, C, D, etc.)
- Include full consensus sequence for each unit
- Units can be different lengths

## Output Files

### 1. Structures File (`*_structures.txt`)

Tab-delimited file with one line per read:

| Column | Description |
|--------|-------------|
| read_name | Unique read identifier |
| strand | Overall orientation (+, -, mixed, unknown) |
| full_structure | Window-by-window assignments (e.g., A-a-a-B-C-U-D) |
| simplified_structure | Collapsed with tandem counts (e.g., A-B-C-U-D) |

### 2. Detailed Scores File (`*_detailed_scores.txt`)

One line per window with alignment metrics for each unit. Useful for:
- Understanding assignment decisions
- Identifying ambiguous regions
- Troubleshooting unexpected results
- Quality control

### 3. Structure Summary File (`*_structure_summary.txt`)

Population-level frequency table:

```
# Structure Summary
# Total reads processed: 3774
# Window size: 1000 bp
# Identity threshold: 40.0%
#
# Strand Distribution:
#   +: 1850 reads (49.0%)
#   -: 1824 reads (48.3%)
#   mixed: 100 reads (2.7%)
#
simplified_structure    count   percentage
C-D                     1250    33.12
C×2-D                   850     22.52
C-D×2                   620     16.43
U                       340     9.01
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `-i, --input` | *Required* | Input FASTA with reads to classify |
| `-u, --units` | *Required* | FASTA with known repeat unit sequences |
| `-o, --output` | *Required* | Output prefix for results files |
| `-w, --window-size` | 1000 | Window size in base pairs |
| `-t, --identity-threshold` | 40.0 | Min % identity for assignment (below = U) |
| `-c, --uppercase-threshold` | 0.3 | Fraction of unit for uppercase (0.3 = first 30%) |
| `-p, --preset` | map-hifi | Minimap2 preset: map-hifi, map-ont, map-pb, asm20 |


### Minimap2 Presets

- `map-hifi`: PacBio HiFi reads (default)
- `map-ont`: Oxford Nanopore reads
- `map-pb`: Older PacBio CLR reads
- `asm20`: Highly similar sequences

## Understanding Results

### Structure Notation

**Uppercase (A, B, C, D)**: Start of tandem repeat or single-copy unit
- Window aligns to first 30% of reference unit

**Lowercase (a, b, c, d)**: Continuation of tandem copy
- Window aligns to last 70% of reference unit

**U**: Undetermined
- Window identity below threshold
- May represent novel sequence, structural variants, or poor quality

**×n**: Tandem notation in simplified structures
- Example: `C×3` = 3 tandem copies of unit C

### Examples

| Full Structure | Simplified | Interpretation |
|----------------|------------|----------------|
| `C-c-c-D-d` | `C-D` | One copy of C, then one copy of D |
| `C-c-c-C-c-c-C-c-c` | `C×3` | Three tandem copies of C |
| `A-a-U-U-U-B-b` | `A-U-B` | Unit A, unassigned region, then unit B |
| `C-c-D-d-C-c-D-d` | `C-D-C-D` | Alternating C and D (not tandem) |
| `C-c-c-C-c-c-D-d-d-D-d` | `C×2-D×2` | Two tandem copies of C, then two of D |

### How Tandem Detection Works

Consider unit C (3000 bp) present as 2 tandem copies:

```
Window 1 (0-1000):     Aligns to C:0-1000     → first 30% → C (uppercase)
Window 2 (1000-2000):  Aligns to C:1000-2000  → past 30%  → c (lowercase)
Window 3 (2000-3000):  Aligns to C:2000-3000  → past 30%  → c (lowercase)
Window 4 (3000-4000):  Aligns to C:0-1000     → first 30% → C (2nd copy!)
Window 5 (4000-5000):  Aligns to C:1000-2000  → past 30%  → c (lowercase)
Window 6 (5000-6000):  Aligns to C:2000-3000  → past 30%  → c (lowercase)

Result: C-c-c-C-c-c (clearly shows 2 tandem repeats)
Simplified: C×2
```

### Strand Awareness

**Why it matters**: Reverse-aligned reads produce windows in 3'→5' order. Without correction, this appears as incorrect tandem repeats.

**Example problem** (without strand correction):
- True structure: `C-D` (single copies)
- Reverse alignment produces: `d-d-D-D-c-c-C-C`
- Would incorrectly appear as: `D×2-C×2`

**Solution**: Tool detects reverse alignment and reverses the structure before interpretation to give correct `C-D` result.

**Strand determination**:
- `+`: Forward strand (>2:1 ratio of forward alignments)
- `-`: Reverse strand (>2:1 ratio of reverse alignments)
- `mixed`: No clear majority (may indicate chimeric reads)
- `unknown`: No assigned windows

## Workflow Example

```bash
# 1. Prepare input files
# - reads.fasta: your long reads
# - units.fasta: your A, B, C, D unit sequences

# 2. Run classifier with default settings
./repeat_classifier_final.py -i reads.fasta -u units.fasta -o results

# 3. Check outputs
head results_structures.txt
head results_structure_summary.txt

# 4. Find specific patterns
grep "C×3-D×2" results_structures.txt

# 5. If results look off, adjust parameters
./repeat_classifier_final.py \
    -i reads.fasta \
    -u units.fasta \
    -o results_v2 \
    -w 2000 \
    -t 45.0
```

## Common Use Cases

### Identify dominant structures
```bash
# View most common patterns
head -20 results_structure_summary.txt
```

### Find all reads with specific structure
```bash
# Get read names with structure "C×3-D"
awk '$4 == "C×3-D" {print $1}' results_structures.txt
```

### Extract high-quality assignments
```bash
# Find reads with no U regions
grep -v "U" results_structures.txt
```

### Count reads by complexity
```bash
# Simple structures (1-2 units)
awk '{split($4, a, "-"); if (length(a) <= 2) print}' results_structures.txt | wc -l

# Complex structures (3+ units)
awk '{split($4, a, "-"); if (length(a) >= 3) print}' results_structures.txt | wc -l
```

### Analyze strand distribution
```bash
# Count by strand
awk 'NR>1 {print $2}' results_structures.txt | sort | uniq -c
```

## Troubleshooting

### Too many U (undetermined) assignments

**Solutions**:
- Lower identity threshold: `-t 35` or `-t 30`
- Verify all expected units are in your units file
- Check if unit sequences are outdated/diverged
- Review detailed scores to see actual identity values

### Window size seems inappropriate

**Solutions**:
- If windows span multiple units → decrease window size
- If alignments are poor quality → increase window size
- Rule of thumb: window = 1/3 to 1/2 of unit length

### Tandem counts seem incorrect

**Solutions**:
- Check `strand` column in structures file
- Adjust uppercase threshold: `-c 0.2` or `-c 0.4`
- Review detailed scores to see alignment positions
- Verify minimap2 preset matches your data type

### minimap2 errors

**Solutions**:
- Verify minimap2 installation: `which minimap2`
- Try different preset: `-p map-ont` or `-p map-pb`
- Check FASTA formatting (no blank lines, proper headers)

## Performance

Processing time scales with:
- Number of reads
- Read length
- Number of units
- Window size (smaller = more windows = slower)

For large datasets (>10,000 reads or very long reads), consider running on a computing cluster. Temporary files are automatically cleaned up.

## Quality Control

**Check structure summary for**:
- High U percentage (>20%) → may need parameter adjustment
- Unexpected dominant structures → review detailed scores
- High mixed strand count → possible chimeric reads

**Use detailed scores to**:
- Identify ambiguous windows (multiple units with similar identity)
- Check alignment quality across reads
- Troubleshoot specific assignments

## Citation

If you use this tool in your research, please cite:

[Your citation information here]

## Author

Freddie  
Transmissible Cancer Group, Murchison Lab  
Department of Veterinary Medicine, University of Cambridge

## License

[Your license here]

## Version History

- **v2.0** (2026-02-06): Added strand awareness and structure summary
- **v1.0** (2026-02-04): Initial release

## Contact

For questions, issues, or feature requests, contact Freddie in the Transmissible Cancer Group.
