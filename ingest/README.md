# Ingest Workflow

The `ingest` directory contains scripts and workflow files for automated data ingestion and preparation. Sequences and corresponding metadata are fetched from GenBank, then cleaned and curated. This process restructures and formats the data into formats suitable for downstream phylogenetic analyses.

---

## Directory Overview

- **`Snakefile`** — Snakemake workflow file that orchestrates the ingest process
- **`bin/`** — Custom scripts used in the ingest workflow
  - [`generate_from_genbank.py`](bin/generate_from_genbank.py) — Downloads and parses GenBank files into required formats
- **`config/`** — Configuration files (e.g., `config.yaml`)
- **`data/`** — Input and output data for the ingest process
- **`source-data/`** — Annotations and geo-location rules
- **`vendored/`** — Vendored scripts from the central Nextstrain ingest repository
- **`workflow/`** — Snakemake rules for the ingest workflow

---

## Getting Started

### Prerequisites

Ensure you have the following installed:
- Python 3.x
- Snakemake
- Required packages: `csvtk`, `nextclade`, `tsv-utils`, `seqkit`, `zip`, `unzip`, `entrez-direct`, `ncbi-datasets-cli` (installable via conda-forge/bioconda)

---

## Workflow Setup

### Step 1: Prepare Reference Files

[Nextclade](https://clades.nextstrain.org/) is used during ingest to align sequences to a reference and assign them to clades. It requires reference files in specific formats (`reference.fasta` and `genome_annotation.gff3`).

#### 1.1 Verify Configuration

Open [`config/config.yaml`](config/config.yaml) and confirm that the `ncbi_taxon_id` is correct for your pathogen.

For EV-A71, this should be:
```yaml
ncbi_taxon_id: 39054
```

#### 1.2 Generate Reference Files from GenBank

Run the `generate_from_genbank.py` script to create properly formatted reference files:

```bash
python3 bin/generate_from_genbank.py --reference "AY426531.1" --output-dir data/references/
```

**During execution**, you will be prompted to specify CDS annotations. Use these inputs for automatic annotation:
   - `[0]` — Select the first option for feature type
   - `[product]` or press Enter for manual selection — Choose how to identify proteins
   - `[2]` — Select the CDS annotation method

**Outputs:**
- `data/references/reference.fasta` — Reference sequence in FASTA format
- `data/references/genome_annotation.gff3` — Genome annotation in GFF3 format
- `data/references/pathogen.json` — Nextclade pathogen configuration

> [!NOTE]  
> For some pathogens, common fields may not be available in the CDS set. In such cases, leave the input blank and manually assign CDS names when prompted.

---

### Step 2: Run the Ingest Workflow

#### 2.1 Make Scripts Executable (if needed)

On Unix-like systems, you may need to make the scripts executable:

```bash
chmod +x ./vendored/* ./bin/*
```

#### 2.2 Execute the Workflow

Run the Snakemake workflow from within the `ingest/` directory:

```bash
snakemake --cores 9 all
```

**What this does:**
1. Fetches sequences and metadata from NCBI Virus based on the taxon ID
2. Runs Nextclade to align sequences and assign clades
3. Cleans and curates metadata
4. Outputs formatted files to `data/`:
   - `sequences.fasta` — Curated sequences
   - `metadata.tsv` — Cleaned and standardized metadata

**Outputs are placed in:**
- `../data/sequences.fasta` (for use by the main workflow)
- `../data/metadata.tsv` (for use by the main workflow)

---

## Workflow Outputs

After successful completion, the following files will be available in the parent `data/` directory:

| File | Description |
|------|-------------|
| `sequences.fasta` | Curated EV-A71 sequences downloaded from NCBI |
| `metadata.tsv` | Standardized metadata including collection dates, locations, and authors |

These files are used as inputs for the main Nextclade dataset workflow.

---

## Customization

### Modifying Metadata Fields

Edit the [config.yaml](config/config.yaml#L100) to customize which metadata fields are extracted and how they are formatted.

### Adjusting Sequence Filters

Modify the Snakemake rules in `workflow/` and `config/config.yaml` to change filtering criteria (e.g., minimum length, date ranges, geographic regions).

### Updating Geo-location Rules

Edit files in `source-data/` to improve geographic location standardization and parsing.

---

## Updating Vendored Scripts

This repository uses [`git subrepo`](https://github.com/ingydotnet/git-subrepo) to manage vendored scripts in the `vendored/` directory. These scripts are maintained in the central Nextstrain ingest repository.

### To Update Vendored Scripts:

1. **Install `git subrepo`:**  
   Follow the [installation guide](https://github.com/ingydotnet/git-subrepo#installation).

2. **Pull the latest changes:**  
   Follow the instructions in [`vendored/README.md`](vendored/README.md#vendoring) to sync with the upstream repository.

3. **Test after updating:**  
   Run the ingest workflow to ensure compatibility with the updated scripts.

---

## Troubleshooting

### Common Issues

**Problem:** Script fails with "permission denied"  
**Solution:** Ensure scripts are executable: `chmod +x ./vendored/* ./bin/*`

**Problem:** No sequences downloaded  
**Solution:** Check that `ncbi_taxon_id` in `config/config.yaml` is correct and that NCBI Virus has sequences available

**Problem:** Metadata parsing errors  
**Solution:** Review geo-location rules in `source-data/` and ensure they cover the locations in your dataset

---

## Additional Resources

- [Nextstrain Ingest Documentation](https://docs.nextstrain.org/projects/ingest/)
- [Nextclade Documentation](https://docs.nextstrain.org/projects/nextclade/)
- [Snakemake Documentation](https://snakemake.readthedocs.io/)

For questions or issues specific to this workflow, please [open an issue](https://github.com/enterovirus-phylo/dataset-template-inferred-root/issues).