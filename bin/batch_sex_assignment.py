#!/usr/bin/env python3
"""
Batch sex assignment runner.
Reads a CSV file where each row is a sample with parameters for assign_sex.py
"""

import csv
import subprocess
import argparse
import sys
from pathlib import Path


def run_sex_assignment(source_organism_id, specimen_id, coverage, ref_version, dry_run=False):
    """
    Run assign_sex.py for a single sample.
    
    Args:
        source_organism_id: Source organism identifier
        specimen_id: Specimen/BAM identifier
        coverage: Path to coverage TSV file
        ref_version: Reference genome version
        dry_run: If True, print command without executing
    
    Returns:
        bool: True if successful, False otherwise
    """
    cmd = [
        "./assign_sex.py",
        "--source_organism_id", source_organism_id,
        "--bam_id", specimen_id,
        "--coverage", coverage,
        "--ref-version", ref_version
    ]
    
    print(f"Processing: {specimen_id} ({source_organism_id})")
    print(f"  Command: {' '.join(cmd)}")
    
    if dry_run:
        print("  [DRY RUN - Command not executed]")
        return True
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        print(f"  ✓ Success")
        return True
    except subprocess.CalledProcessError as e:
        print(f"  ✗ Error: {e.stderr}")
        return False
    except FileNotFoundError:
        print(f"  ✗ Error: assign_sex.py not found in PATH")
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Batch sex assignment from CSV file"
    )
    parser.add_argument(
        "csv_file",
        help="CSV file with samples (columns: source_organism_id, specimen_id, coverage, ref_version)"
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands without executing"
    )
    parser.add_argument(
        "--delimiter",
        default=",",
        help="CSV delimiter (default: comma)"
    )
    parser.add_argument(
        "--skip-header",
        action="store_true",
        help="Skip first row (header)"
    )
    
    args = parser.parse_args()
    
    # Check if CSV file exists
    csv_path = Path(args.csv_file)
    if not csv_path.exists():
        print(f"Error: CSV file not found: {args.csv_file}", file=sys.stderr)
        sys.exit(1)
    
    # Read and process CSV
    total_samples = 0
    successful = 0
    failed = 0
    
    try:
        with open(csv_path, 'r') as f:
            reader = csv.reader(f, delimiter=args.delimiter)
            
            # Skip header if requested
            if args.skip_header:
                next(reader)
            
            for row_num, row in enumerate(reader, start=1):
                # Skip empty rows
                if not row or all(cell.strip() == '' for cell in row):
                    continue
                
                # Expect at least 4 columns
                if len(row) < 4:
                    print(f"Warning: Row {row_num} has fewer than 4 columns, skipping")
                    continue
                
                source_organism_id = row[0].strip()
                specimen_id = row[1].strip()
                coverage = row[2].strip()
                ref_version = row[3].strip()
                
                total_samples += 1
                
                if run_sex_assignment(source_organism_id, specimen_id, coverage, ref_version, args.dry_run):
                    successful += 1
                else:
                    failed += 1
        
        # Summary
        print("\n" + "="*60)
        print(f"Summary: {total_samples} samples processed")
        print(f"  Successful: {successful}")
        print(f"  Failed: {failed}")
        print("="*60)
        
        sys.exit(0 if failed == 0 else 1)
    
    except Exception as e:
        print(f"Error reading CSV file: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
