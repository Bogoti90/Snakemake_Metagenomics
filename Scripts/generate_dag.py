#!/usr/bin/env python3

"""
Generate a clear and well-formatted DAG visualization for the metagenomics pipeline.
Usage: python generate_dag.py
"""

import os
import subprocess
from textwrap import dedent

# Colors for different pipeline stages
colors = {
    'host_removal': '#FFB6C1',      # Light pink
    'taxonomic_classification': '#98FB98',  # Pale green
    'assembly': '#87CEEB',          # Sky blue
    'mapping': '#DDA0DD',           # Plum
    'cyb_analysis': '#F0E68C',      # Khaki
    'visualization': '#FFA07A',      # Light salmon
    'stats': '#B8860B',             # Dark golden rod
    'extraction': '#E6E6FA'         # Lavender
}

# Define node attributes with descriptions
node_attrs = {
    # Host removal related rules
    'build_bowtie2_db': {
        'fillcolor': colors['host_removal'],
        'label': 'Build Host\nDatabase'
    },
    'map_reads': {
        'fillcolor': colors['host_removal'],
        'label': 'Map Reads to\nHost'
    },
    'sam_to_bam': {
        'fillcolor': colors['host_removal'],
        'label': 'Convert SAM\nto BAM'
    },
    'extract_unmapped_reads': {
        'fillcolor': colors['host_removal'],
        'label': 'Extract Non-host\nReads'
    },
    'sort_reads': {
        'fillcolor': colors['host_removal'],
        'label': 'Sort Reads'
    },
    'bam_to_fastq': {
        'fillcolor': colors['host_removal'],
        'label': 'Convert to\nFASTQ'
    },
    
    # Taxonomic classification rules
    'run_kaiju': {
        'fillcolor': colors['taxonomic_classification'],
        'label': 'Kaiju\nClassification'
    },
    'kaiju2table_summary': {
        'fillcolor': colors['taxonomic_classification'],
        'label': 'Kaiju\nSummary'
    },
    'run_centrifuge': {
        'fillcolor': colors['taxonomic_classification'],
        'label': 'Centrifuge\nClassification'
    },
    
    # Assembly related rules
    'run_metaspades': {
        'fillcolor': colors['assembly'],
        'label': 'MetaSPAdes\nAssembly'
    },
    'build_assembly_database': {
        'fillcolor': colors['assembly'],
        'label': 'Build Assembly\nDatabase'
    },
    
    # Mapping related rules
    'map_reads_to_assembly': {
        'fillcolor': colors['mapping'],
        'label': 'Map Reads to\nAssembly'
    },
    'index_assembly': {
        'fillcolor': colors['mapping'],
        'label': 'Index\nAssembly'
    },
    'sam2_to_bam2': {
        'fillcolor': colors['mapping'],
        'label': 'Convert Assembly\nSAM to BAM'
    },
    'sort_bam2': {
        'fillcolor': colors['mapping'],
        'label': 'Sort Assembly\nBAM'
    },
    'index_bam': {
        'fillcolor': colors['mapping'],
        'label': 'Index BAM'
    },
    'index_bam2': {
        'fillcolor': colors['mapping'],
        'label': 'Index Filtered\nBAM'
    },
    
    # CyB analysis rules
    'extract_aligned_reads_for_mapped_reads': {
        'fillcolor': colors['cyb_analysis'],
        'label': 'Extract CytB\nReads'
    },
    'bam_to_fastq_for_cyb': {
        'fillcolor': colors['cyb_analysis'],
        'label': 'Convert CytB\nto FASTQ'
    },
    'run_metaspades_take_two': {
        'fillcolor': colors['cyb_analysis'],
        'label': 'CytB\nAssembly'
    },
    
    # Visualization and stats rules
    'centrifuge_to_krona': {
        'fillcolor': colors['visualization'],
        'label': 'Generate\nKrona Chart'
    },
    'generate_assembly_stats': {
        'fillcolor': colors['stats'],
        'label': 'Generate\nStats'
    },
    
    # Additional extraction rules
    'filter_reads_by_taxon': {
        'fillcolor': colors['extraction'],
        'label': 'Filter Reads\nby Taxon'
    },
    'extract_reads': {
        'fillcolor': colors['extraction'],
        'label': 'Extract\nFiltered Reads'
    }
}

def create_legend():
    """Create a legend subgraph."""
    legend = dedent("""
    subgraph cluster_legend {
        label="Pipeline Stages";
        fontsize=12;
        fontname=sans;
        style=rounded;
        color=gray;
        margin=20;
        
        node[shape=box, style=filled, margin=0.2];
        
        legend_host[label="Host Removal", fillcolor="%(host_removal)s"];
        legend_tax[label="Taxonomic Classification", fillcolor="%(taxonomic_classification)s"];
        legend_assembly[label="Assembly", fillcolor="%(assembly)s"];
        legend_mapping[label="Mapping", fillcolor="%(mapping)s"];
        legend_cyb[label="CytB Analysis", fillcolor="%(cyb_analysis)s"];
        legend_viz[label="Visualization", fillcolor="%(visualization)s"];
        legend_stats[label="Statistics", fillcolor="%(stats)s"];
        legend_extract[label="Read Extraction", fillcolor="%(extraction)s"];
        
        {rank=same; legend_host legend_tax legend_assembly legend_mapping}
        {rank=same; legend_cyb legend_viz legend_stats legend_extract}
    }
    """) % colors
    return legend

def combine_dag_with_style(raw_dag, styled_dag):
    """Combine raw DAG structure with styling."""
    # Extract edges from raw DAG
    with open(raw_dag, 'r') as f:
        raw_content = f.read()
    
    # Extract edges (lines containing '->') from raw DAG
    edges = [line.strip() for line in raw_content.split('\n') if '->' in line]
    
    # Read styled DAG
    with open(styled_dag, 'r') as f:
        styled_content = f.readlines()
    
    # Find where to insert edges (before the last '}')
    insert_position = -1
    for i, line in enumerate(styled_content):
        if line.strip() == '}':
            insert_position = i
            break
    
    # Combine content
    final_content = (
        ''.join(styled_content[:insert_position]) +
        '\n    # Edges\n    ' +
        '\n    '.join(edges) +
        '\n' +
        ''.join(styled_content[insert_position:])
    )
    
    return final_content

def generate_dag():
    """Generate the DAG visualization."""
    print("Generating DAG visualization...")
    
    # Create raw DAG using Snakemake
    print("1. Generating raw DAG structure...")
    subprocess.run(["snakemake", "--dag"], stdout=open("dag_raw.dot", "w"))
    
    # Create styled DOT file
    print("2. Creating styled DOT file...")
    dot_content = """digraph snakemake_dag {
    graph[bgcolor=white, margin=0.5, pad=0.5];
    node[shape=box, style=filled, fontname=sans, fontsize=10, penwidth=2, margin=0.15];
    edge[color=gray50, arrowsize=0.7, arrowhead=vee, penwidth=1];
    
    # Node styling
"""
    
    # Add node style configurations
    for node, attrs in node_attrs.items():
        attrs_str = ",".join(f'{k}="{v}"' for k, v in attrs.items())
        dot_content += f'    {node}[{attrs_str}];\n'
    
    # Add legend
    dot_content += create_legend()
    dot_content += "}\n"
    
    with open("dag_styled.dot", "w") as f:
        f.write(dot_content)
    
    # Combine raw DAG with styling
    print("3. Combining structure with styling...")
    final_content = combine_dag_with_style("dag_raw.dot", "dag_styled.dot")
    
    with open("dag_final.dot", "w") as f:
        f.write(final_content)
    
    # Generate final visualizations
    print("4. Rendering final visualizations...")
    subprocess.run(["dot", "-Tpng", "-Gdpi=300", "dag_final.dot", "-o", "dag.png"])
    subprocess.run(["dot", "-Tsvg", "dag_final.dot", "-o", "dag.svg"])
    
    # Cleanup temporary files
    print("5. Cleaning up temporary files...")
    for temp_file in ["dag_raw.dot", "dag_styled.dot", "dag_final.dot"]:
        if os.path.exists(temp_file):
            os.remove(temp_file)
    
    print("\nGenerated DAG visualizations:")
    print("- dag.png (PNG format)")
    print("- dag.svg (SVG format)")

if __name__ == "__main__":
    generate_dag()