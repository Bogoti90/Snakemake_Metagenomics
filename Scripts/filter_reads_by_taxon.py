#!/usr/bin/env python3
"""
Script to filter reads by taxon from Centrifuge output.
Usage: python filter_reads_by_taxon.py <centrifuge_report> <taxon_id> <output_file>
"""

import sys
import pandas as pd
from snakemake.shell import shell
from snakemake.utils import logger

def parse_centrifuge_report(report_file):
    """Parse the Centrifuge report file."""
    try:
        df = pd.read_csv(report_file, sep='\t')
        return df
    except Exception as e:
        logger.error(f"Error reading Centrifuge report: {e}")
        sys.exit(1)

def filter_by_taxon(df, taxon_id):
    """Filter reads by the specified taxon ID."""
    try:
        # Convert taxon_id to integer for comparison
        taxon_id = int(taxon_id)
        
        # Filter rows where taxID matches
        filtered_df = df[df['taxID'] == taxon_id]
        
        if filtered_df.empty:
            logger.warning(f"No reads found for taxon ID {taxon_id}")
            return None
            
        return filtered_df['readID'].unique()
        
    except ValueError:
        logger.error(f"Invalid taxon ID: {taxon_id}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Error filtering reads: {e}")
        sys.exit(1)

def write_read_ids(read_ids, output_file):
    """Write filtered read IDs to output file."""
    try:
        with open(output_file, 'w') as f:
            for read_id in read_ids:
                f.write(f"{read_id}\n")
    except Exception as e:
        logger.error(f"Error writing output file: {e}")
        sys.exit(1)

def main():
    """Main function to process the filtering."""
    # Get input and output files from snakemake
    report_file = snakemake.input.report
    output_file = snakemake.output.ids
    taxon_id = snakemake.params.taxon_of_interest

    logger.info(f"Processing Centrifuge report: {report_file}")
    logger.info(f"Filtering for taxon ID: {taxon_id}")

    # Parse the report
    df = parse_centrifuge_report(report_file)
    
    # Filter reads
    filtered_reads = filter_by_taxon(df, taxon_id)
    
    if filtered_reads is not None:
        # Write results
        write_read_ids(filtered_reads, output_file)
        logger.info(f"Found {len(filtered_reads)} reads for taxon {taxon_id}")
    else:
        # Create empty output file if no reads found
        open(output_file, 'w').close()
        logger.warning("Created empty output file")

if __name__ == "__main__":
    main()