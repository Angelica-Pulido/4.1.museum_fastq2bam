#!/usr/bin/env python3

"""
Script: rg_parser.py

Description:
    This script processes FASTQ files to extract read group (RG) information
    for downstream analysis (e.g., GATK, SAMtools). It parses sequencing metadata
    (run ID, flow cell, lane, index) from the first read of the FASTQ file and 
    writes the read group string to an output file.

Usage example:
    python3 rg_parser.py \
        --source_organism_id HG001 \
        --specimen_id Sample123 \
        --run_id Run42 \
        --lane 1 \
        --r1 /path/to/read1.fastq.gz
"""


# Import necessary libraries
import os.path
import argparse
import logging
import sys
import gzip

# Parse input arguments
parser = argparse.ArgumentParser(description="Get FASTQ files")
parser.add_argument('--source_organism_id', required=True, help='RingID.')
parser.add_argument('--specimen_id', required=True, help='Specimen identifier, defined as the library name.')
parser.add_argument('--run_id', required=True, help='Sequencing run ID provided by the sequencing platform.')
parser.add_argument('--lane', required=True, help='Flow cell lane the FASTQ was obtained from.')
parser.add_argument('--r1', required=True, help='Path to Read 1 FASTQ file')
args = parser.parse_args()


organism_id = args.source_organism_id
specimen_id = args.specimen_id
run_id = args.run_id
lane  = args.lane
R1 = args.r1

# Validate essential options
if specimen_id is None:
    raise RuntimeError("Error! no specimen id provided (--specimen_id)")
if run_id is None:
    raise RuntimeError("Error! no GTF run id provided (--run_id)")
if R1 is None:
    raise RuntimeError("Error! no forward read provided (--r1)")


def get_rg(sample, specimen, run, lane, read):
    """
    Extract read group (RG) information from a gzipped FASTQ file.

    Args:
        sample (str): Organism/sample identifier (e.g., RingID).
        specimen (str): Specimen identifier (library name).
        run (str): Sequencing run ID.
        lane (str): Flow cell lane number.
        read (str): Path to the Read 1 FASTQ file (.fastq.gz).

    Returns:
        str: A formatted string representing the read group information:
             "<specimen>|<seqrun>.<flowcell_id>.<lane>|<index>|<sample>"

    Behavior:
        - Reads the first header line of the FASTQ file.
        - Extracts sequencer run ID, flow cell ID, flow cell lane, and index.
        - Writes the RG string into an output file:
          "<specimen>_<run>_<lane>_RG_info.txt".
    """
    with gzip.open(read, "rt") as handle:
        header = handle.readline().strip()
        fields = header.split(":")

        if len(fields) < 10:
            raise ValueError(f"Unexpected FASTQ header format: {header}")

        seqrun = fields[1]      # Sequencer run ID
        fc_id = fields[2]       # Flow cell ID
        fc_lane = fields[3]     # Flow cell lane
        idx = fields[9]         # Index sequence

    read_group_info = f"{specimen}|{seqrun}.{fc_id}.{fc_lane}|{idx}|{sample}"

    out_file = f"{specimen}_{run}_{lane}_RG_info.txt"
    with open(out_file, "w") as f:
        f.write(read_group_info + "\n")

    return read_group_info


def main(sample, specimen, run, lane, read):
    """
    Main entry point: Extracts and saves read group information.

    Args:
        sample (str): Organism/sample identifier.
        specimen (str): Specimen/library identifier.
        run (str): Sequencing run ID.
        lane (str): Flow cell lane.
        read (str): Path to Read 1 FASTQ file (.fastq.gz).

    Behavior:
        - Calls get_rg() to extract metadata from the FASTQ header.
        - Saves the read group information into a text file.
    """
    rg_info = get_rg(sample, specimen, run, lane, read)
    print(f"Read group info extracted: {rg_info}")


# -----------------------------
# Run script
# -----------------------------
if __name__ == "__main__":
    main(organism_id, specimen_id, run_id, lane, R1)
