#!/usr/bin/env python3

"""
Determine individual's sex based on genomic coverage stats.
"""

# Import necessary libraries

import pandas as pd
import os.path
import argparse
import logging
import sys

# Parse input arguments
parser = argparse.ArgumentParser(description='Determine individual\'s sex based on genomic coverage stats.')

parser.add_argument('--source_organism_id', required=True, help='RingID or source organism ID.')
parser.add_argument('--bam_id', required=True, help='BAM identifier.')
parser.add_argument('--coverage', type=str, required=True, help='Path to individual coverage stats')
parser.add_argument('--ref-version', type=str, dest= 'refv', required=True, help='Reference genome version; e.g. "2020"')
args = parser.parse_args()

source_organism_id = args.source_organism_id
bam_id = args.bam_id
cov = args.coverage
refg = args.refv


# Validate essential options
if source_organism_id is None:
    raise RuntimeError("Error! no Source Organism ID provided (--source_organism_id)")
if bam_id is None:
    raise RuntimeError("Error! no BAM ID provided (--bam_id)")
if cov is None:
    raise RuntimeError("Error! no coverage file provided (--coverage)")
if refg is None:
    raise RuntimeError("Error! no reference genome version provided (--refref-version)")


def ref2020(cov):
    """
    Determine sex based on coverage depth using Tyto alba's 2020 reference genome.
    
    Parameters:
        cov (DataFrame): Coverage stats loaded from file.
        
    Returns:
        str: 'M' for male, 'F' for female.
    """
    # Get mean depth for the sex chromosome (Super-Scaffold_13)
    scaffold13_depth = cov[cov['#rname'] == 'Super-Scaffold_13']['meandepth'].mean()
    # Filter out sex chromosome and PAR from autosomal chromosomes
    cov_autosomal = cov[~cov['#rname'].isin(['Super-Scaffold_13', 'Super-Scaffold_42'])]
    # Compute average autosomal depth
    autosomal_depth = cov_autosomal['meandepth'].mean()
    # Compute depth ratio between sex chromosome and autosomes
    depth_ratio = scaffold13_depth/autosomal_depth
    if depth_ratio > 0.65:
        indsex = 'M'
    else:
        indsex = 'F'
    return(indsex)


def main(ringID, bamid, cov_file, ref):
    """
    Main function to process coverage file and determine sex.
    
    Parameters:
        cov_file (str): Path to input coverage file.
        ref (str): Reference genome version.
    """
    # Load coverage data
    cov = pd.read_csv(cov_file, sep = "\t", header = 0)
  
    # Determine sex based on reference genome
    if ref == "2020":
        sex = ref2020(cov)
    else:
        # Placeholder for handling other reference genomes
        sex = "NA"
    # Write output to a TSV file
    with open(bamid + "_ind_sex.tsv", "w") as f:
        f.write(f"{ringID} {sex}\n")



if __name__ == "__main__":
    main(source_organism_id, bam_id, cov, refg)
