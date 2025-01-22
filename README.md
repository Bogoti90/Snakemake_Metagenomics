# Wildlife Metagenomics Analysis Pipeline

A Snakemake workflow for analyzing metagenomic data from wildlife samples collected in the Afar region of Ethiopia.

## Overview

This pipeline performs:
- Host read removal
- Taxonomic classification (Kaiju, Centrifuge)
- Metagenomic assembly (MetaSPAdes)
- Read mapping and coverage analysis
- Cytochrome B gene analysis
- Taxonomic visualization (Krona)

## Requirements

### Software Dependencies
- Python ≥3.9
- Snakemake ≥8.11.6
- Mamba/Conda
- Required tools (automatically installed via conda environments):
  - Bowtie2
  - Samtools
  - SPAdes
  - Kaiju
  - Centrifuge
  - Krona
  - seqtk

### Input Data Requirements
- Paired-end Illumina reads (FASTQ format)
- Host genome reference
- Kaiju database
- Centrifuge database

## Installation

### 1. Install Mamba (Recommended)
```bash
# Install Mambaforge if you don't have it
wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
bash Mambaforge-Linux-x86_64.sh

# Or install mamba into existing conda
conda install -n base -c conda-forge mamba
```

### 2. Create Snakemake Environment
```bash
# Create and activate environment
mamba create -c conda-forge -c bioconda -n snakemake snakemake
mamba activate snakemake
```

### 3. Clone Repository
```bash
git clone https://github.com/your-username/wildlife-metagenomics.git
cd wildlife-metagenomics
```

## Usage

### 1. Configure Pipeline
Edit `config.yaml` to set your parameters:
```yaml
threads:
  bowtie2: 8
  spades: 48
  # ... other parameters
```

### 2. Prepare Input Data
- Place paired-end FASTQ files in the working directory
- Ensure reference databases are available and paths are correctly specified

### 3. Run Pipeline
```bash
# Dry run to check workflow
snakemake -n

# Run pipeline with 4 cores
snakemake --cores 4 --use-conda

# Run on a cluster (e.g., SLURM)
snakemake --cluster "sbatch --mem=100g" --jobs 100
```

## Pipeline Steps

1. **Host Read Removal**
   - Maps reads against host genome
   - Filters unmapped reads for further analysis

2. **Taxonomic Classification**
   - Kaiju classification of non-host reads
   - Centrifuge parallel classification
   - Krona visualization of results

3. **Assembly and Analysis**
   - MetaSPAdes assembly of filtered reads
   - Read mapping to assemblies
   - Coverage analysis

4. **Cytochrome B Analysis**
   - Extraction of CytB reads
   - Targeted assembly
   - Taxonomic assignment

## Output Directory Structure
```
├── filtered_reads/           # Host-filtered reads
├── kaiju_outputs/           # Kaiju classification results
├── centrifuge_reports/      # Centrifuge classification
├── metaspades_output/       # Assemblies
├── krona_charts/           # Interactive visualizations
└── assembly_stats/         # Assembly statistics
```

## Troubleshooting

Common issues and solutions:
1. **Memory Issues**: Adjust memory in config.yaml
2. **Database Errors**: Verify database paths and formats
3. **Missing Files**: Check input file naming conventions

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.


## Contact

* **Brian M Ogoti**
  * Center for Epidemiological Modelling and Analysis (CEMA)
  * Twitter: [@diyobraz2](https://twitter.com/diyobraz2)
  * Web: [CEMA Website](URL-to-CEMA)

