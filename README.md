# Nextclade Workflow for Enterovirus A71

This repository contains a robust, reproducible workflow for building a custom [Nextclade](https://github.com/nextstrain/nextclade) dataset for Enterovirus A71 (EV-A71). It enables you to generate reference and annotation files, download and process sequence data, infer an ancestral sequence, and create all files needed for Nextclade analyses and visualization.

---

## Quick Start

```bash
# 1. Set up folders
mkdir -p dataset data ingest resources results scripts

# 2. Generate reference files
python3 scripts/generate_from_genbank.py --reference "U22521.1" --output-dir dataset/

# 3. Configure pathogen.json (edit manually)

# 4. If first time, enable inference in Snakefile:
# Set INFERRENCE_RERUN = True

# 5. Run workflow
snakemake --cores 9 all --config static_inference_confirmed=true
```

See detailed instructions below for each step.

---

## Folder Structure

Follow the [Nextclade example workflow](https://github.com/nextstrain/nextclade_data/tree/master/docs/example-workflow) or use the structure below:

```bash
mkdir -p dataset data ingest resources results scripts
```

---

## Workflow Overview

This workflow is composed of several modular steps:

1. **Reference Generation**  
   Extracts relevant reference and annotation files from GenBank.
2. **Dataset Ingest**  
   Downloads and processes sequences and metadata from NCBI Virus.
3. **Inferred Ancestral Root (Recommended)**  
   Uses outgroup rooting to infer a dataset-specific ancestral sequence. This is rooted on a *Static Inferred Ancestor* — a phylogenetically reconstructed sequence at the MRCA (most recent common ancestor) of the ingroup, which provides a stable, biologically accurate reference point for mutation and clade assignments. This approach addresses the issue that the BrCr reference differs substantially from currently circulating strains.
4. **Augur Phylogenetics & Nextclade Preparation**  
   Builds trees rooted on the inferred ancestor, prepares multiple sequence alignments, and generates all files required for Nextclade and Auspice.
5. **Visualization & Analysis**  
   Enables both command-line and web-based Nextclade analyses, including local dataset hosting.

---

## Setup Instructions

### 1. Generate Reference Files

Run the script to extract the reference FASTA and genome annotation from GenBank:

```bash
python3 scripts/generate_from_genbank.py --reference "AY426531.1" --output-dir dataset/
```

During the script execution, follow the prompts for CDS annotation selection.

- `[0]`
- `[product]` or `[leave empty for manual choice]` to select proteins.
- `[2]`.

**Outputs:**

- `dataset/reference.fasta`
- `dataset/genome_annotation.gff3`

---

### 2. Configure `pathogen.json`

Edit `pathogen.json` to:

- Reference your generated files (`reference.fasta`, `genome_annotation.gff3`)
- Update metadata and QC settings as needed  

> [!WARNING]  
> If QC is not set, Nextclade will skip quality checks.

See the [Nextclade pathogen config documentation](https://docs.nextstrain.org/projects/nextclade/en/latest/user/input-files/05-pathogen-config.html) for details.

---

### 3. Prepare GenBank Reference

Copy your GenBank file to `resources/reference.gb` and edit it to ensure compatibility with the workflow.

**Important requirements:**

- Each coding sequence (CDS) must have either a `product` or `gene` name present
- The annotation keys must **match exactly** between `reference.gb` and `genome_annotation.gff3`
- Use simple, consistent names (e.g., `product="VP1"` instead of `product="VP1_protein"`)
- Remove any genes that are not relevant for your dataset

> [!WARNING]  
> Mismatched or inconsistent gene names will cause `augur ancestral` to fail, as it cannot match features across files. Ensure your protein names match those defined in the `GENES` list in the [Snakefile](/Snakefile#L4).

---

### 4. Update the `Snakefile`

- Adjust the workflow parameters and file paths as needed for your dataset.
- Ensure required files are available:
  - `data/sequences.fasta`
  - `data/metadata.tsv`
  - `resources/auspice_config.json`

Sequences and metadata can be downloaded automatically via the ingest process (see below).

---

## Subprocesses

### Ingest

Automates downloading of EV-A71 sequences and metadata from NCBI Virus.  
See [ingest/README.md](ingest/README.md) for specifics.

**Required packages:**  
`csvtk, nextclade, tsv-utils, seqkit, zip, unzip, entrez-direct, ncbi-datasets-cli` (installable via conda-forge/bioconda)

---

### Inferred Ancestral Root with Outgroup Rooting (Recommended)

The `inferred-root/` directory contains a reproducible pipeline that uses **outgroup rooting** to infer a dataset-specific ancestral sequence for EV-A71. This method:

- **Builds a phylogenetic tree** including both EV-A71 sequences (ingroup) and related enterovirus sequences (outgroup)
- **Roots the tree on the outgroup** to establish correct evolutionary directionality
- **Extracts the ancestral sequence** at the MRCA of all EV-A71 sequences
- **Fills gaps** with reference nucleotides to ensure a complete, biologically plausible genome

This **Static Inferred Ancestor** serves as the root of your Nextclade dataset, providing:

- More accurate mutation calls relative to a realistic EV-A71 ancestor
- A stable reference that better represents EV-A71 diversity than the distant BrCr sequence

#### Configuration

The workflow has two key parameters in the main `Snakefile`:

- `STATIC_ANCESTRAL_INFERRENCE = True` — enables using the inferred root (default: `True`)
- `INFERRENCE_RERUN = False` — controls whether to regenerate the inferred root (default: `False`)

#### For Regular Dataset Builds

Use the existing inferred root:

```bash
snakemake --cores 9 all
```

#### To Regenerate the Inferred Root

When you need to regenerate with new data or updated outgroups:

1. Set `INFERRENCE_RERUN = True` in the Snakefile
2. Run the workflow:

   ```bash
   snakemake --cores 9 all --config static_inference_confirmed=true
   ```

3. The workflow will:
   - Clean previous results in `inferred-root/results/`
   - Run the full inference pipeline with your current sequences
   - Generate a new `resources/inferred-root.fasta`
   - Incorporate it into the dataset build
4. After successful regeneration, set `INFERRENCE_RERUN = False` for future runs

> [!WARNING]  
> Setting `INFERRENCE_RERUN = True` will **overwrite** your existing `resources/inferred-root.fasta` file and clear `inferred-root/results/`. Only use this when you want to regenerate the root with updated data.

> [!NOTE]  
> - **First-time users:** If `resources/inferred-root.fasta` doesn't exist, you must set `INFERRENCE_RERUN = True` initially.
> - **To disable this feature:** Set `STATIC_ANCESTRAL_INFERRENCE = False` and change `ROOTING` parameter (e.g., `ROOTING="mid_point"`).
> - **Outgroup configuration:** Sequences are in `resources/outgroup/`; update the `OUTGROUP` list in `inferred-root/Snakefile` to modify which species are used.

**See:** [`inferred-root/README.md`](inferred-root/README.md) for technical details and the complete workflow.

---

## Star-like root

The reference tree for the Nextclade dataset is not a plain phylogeny. Before it reaches `augur refine`, it is deliberately restructured into a **star-like tree** by [scripts/star_like_rooter.py](scripts/star_like_rooter.py).

**Why.** Recombinant clades — and any strain flagged as a known recombinant — have no single meaningful position in a bifurcating tree: their genomes are mosaics of multiple parents, so their placement is ambiguous and can distort branch lengths and topology for the rest of the tree. Each recombinant clade is therefore cut out and reattached directly to the root, radiating from it like a star. The non-recombinant backbone, including reference and outgroup strains, keeps its original structure and position.

**How it works.** On each run:

1. **Label the tips.** Every accession in `resources/clades_metadata.tsv` is labeled with its clade (e.g. `C4`). Strains listed in the [recombinants file](resources/recombinants.tsv) are labeled `RFs`; tips with no clade recorded at all are labeled `unassigned`.
2. **Find the reattachment points.** The tree is traversed bottom-up. For every recombinant-classified tip (`RFs`, plus any clade passed to `--recombinant_clades`), the script climbs to the *highest* ancestor whose descendants are all recombinant-classified. This merges adjacent and nested recombinant clades into a single reattachment point rather than scattering many small ones.
3. **Reattach.** Each top-level recombinant ancestor is detached from its parent — pruning any parent left childless — and reattached as a direct child of the root, with its original root-to-node distance preserved as the new branch length.
4. **Log the result**, e.g.:

```text
   Star-like restructuring: 6 recombinant clades reattached to root, 224 strains moved, backbone preserved
```

> [!WARNING]
> **An unexpectedly high number of reattached clades or strains is usually a metadata problem, not a recombination-heavy dataset.** It most often means the tree's tips and `resources/clades_metadata.tsv` have fallen out of sync — for example, the tree was rebuilt from a different or newly subsampled strain set, but the clade metadata was not regenerated. Tips with no matching row default to `unassigned` (non-recombinant), which can silently **fragment a recombinant clade** into several pieces that are each reattached separately instead of merged into one.
>
> The script logs how many tips fell into this state:
>
> ```text
> N/M (x.x%) tree tips have NO row in --input_clades and will be treated as non-recombinant 'unassigned' strains...
> Missing accessions (up to 10 shown): [...]
> ```
>
> Whenever the reattachment count looks high, **check the log for this warning first**, look up the accessions it lists, and add the missing rows to `resources/clades_metadata.tsv` with correct clade assignments before trusting the resulting star tree.

### Template for Other Enteroviruses

If you want to apply this approach to other enterovirus types (e.g., EV-A71, CVA16), a [Nextclade Dataset Template for Inferred Root](https://github.com/enterovirus-phylo/dataset-template-inferred-root) is available and recommended for reuse.

---

## Running the Workflow

To generate the Auspice JSON and Nextclade dataset:

```bash
snakemake --cores 9 all
```

This will use the existing inferred root (see [Inferred Ancestral Root](#inferred-ancestral-root-with-outgroup-rooting-recommended) section above for regeneration instructions).

The workflow will:
- Build the reference tree rooted on the inferred ancestor
- Produce the Nextclade dataset in `out-dataset/`
- Run Nextclade on example sequences
- Output results to `test_out/` (alignment, translations, summary TSV)

**Key Snakefile parameters:**
- `ROOTING = "ancestral_sequence"` — roots tree on the inferred ancestor
- `STATIC_ANCESTRAL_INFERRENCE = True` — enables inferred root in the dataset (default)
- `INFERRENCE_RERUN = False` — set to `True` only when regenerating the root (default: `False`)

### Labeling Mutations of Interest

To label mutations of interest, execute the `mutLabels` rule as a standalone instance. They will be added to the `out-dataset/pathogen.json` file.

---

## Visualizing Your Custom Nextclade Dataset

To use the dataset in Nextclade Web, serve it locally:

```bash
serve --cors out-dataset -l 3000
```

Then open:

```
https://master.clades.nextstrain.org/?dataset-url=http://localhost:3000
```

- Click "Load example", then "Run"
- You may want to reduce "Max. nucleotide markers" to 500 under "Settings" → "Sequence view" to optimize performance

---

## Author & Contact

- Maintainers: Nadia Neuner-Jehle, Alejandra González-Sánchez and Emma B. Hodcroft ([eve-lab.org](https://eve-lab.org/))
- For questions or suggestions, please [open an issue](https://github.com/enterovirus-phylo/dataset-template-inferred-root/issues) or email: eve-group[at]swisstph.ch

## Troubleshooting and Further Help

- For issues, see the [official Nextclade documentation](https://docs.nextstrain.org/projects/nextclade/en/stable/index.html#) or [open an issue](https://github.com/enterovirus-phylo/dataset-template-inferred-root/issues).
- For details on the inferred root workflow, see [`inferred-root/README.md`](inferred-root/README.md).
- For adapting to other enteroviruses, see the [dataset-template-inferred-root](https://github.com/enterovirus-phylo/dataset-template-inferred-root).

---

This guide provides a structured, scalable approach to building and using high-quality Nextclade datasets for EV-A71 — and can be adapted for other enterovirus types as well.

---

## Task List

**Completed:**
- [x] Integrate ancestral inferred-root into workflow
- [x] Integrate epitope mutation information as tree coloring and/or display in the Nextclade results table
- [x] Generate the inferred ancestral sequence from the outgroup-rooted tree (see [mpox](https://github.com/corneliusroemer/alignment-ref-construction/tree/e7ae8c7ed51d9754212a113c8d795180f1b80410) for technical details)

**Documentation & Visualization:**
- [ ] Document outgroup selection and validation — explain which enterovirus species are used as outgroups and phylogenetic justification
- [ ] Add workflow diagram — visual representation showing when INFERRENCE_RERUN triggers the inferred-root sub-workflow
- [ ] Add troubleshooting for INFERRENCE_RERUN — common errors when regenerating (missing outgroups, alignment failures, etc.)
- [] Star-like root: documentation


**Analysis & Validation:**
- [] Star-like root: use tide-tree instead?
- [] Validate clade assignment of fragmented sequences in Nextclade (`testing/`)
- [] Review and validate EV-A71 nomenclature, including robustness with recombinant sequences
- [] Create test dataset — small example demonstrating the full inferred-root workflow end-to-end
- [ ] Document when to regenerate inferred root — guidelines on how often to rerun with new data
- [ ] Validate rooting stability — test sensitivity of inferred root to outgroup selection and subsampling
