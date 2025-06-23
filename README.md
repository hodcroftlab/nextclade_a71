# Nextclade Setup for Enterovirus A71

## Folder Structure
First, create the necessary folder structure as shown in the [example workflow](https://github.com/nextstrain/nextclade_data/tree/master/docs/example-workflow):

```
dataset/
profiles/
resources/
rules/
scripts/
results/
```

You can create these directories using the following command:
```bash
mkdir -p dataset profiles resources rules scripts results
```

---
## Steps to Set Up The Workflow

### 1. Run `generate_from_genbank.py`
This script (located in `scripts/`) generates reference files from GenBank.

Run the following command:
```bash
python3 scripts/generate_from_genbank.py --reference "U22521.1" --output-dir dataset/
```

During execution, you may be asked to provide CDS annotations. You can use the following codes to specify the CDS automatically:
   - `[0]`
   - `[product]` or `[leave empty for manual choice]` to select proteins.
   - `[2]`.

The script will generate:
- `dataset/reference.fasta`
- `dataset/genome_annotation.gff3`

---
### 2. Update `pathogen.json`
Modify `pathogen.json` to:
- Ensure file names match the generated reference files.
- Update attributes as needed.
- Adjust the Quality Control (QC) settings if necessary. If QC is not configured, Nextclade will not perform any checks.

For more details on configuration, refer to the [Nextclade documentation](https://docs.nextstrain.org/projects/nextclade/en/latest/user/input-files/05-pathogen-config.html).

---
### 3. Prepare `reference.gb`
- Copy the `reference.gb` file into the `resources/` directory.
- Modify protein names as needed to match your requirements.

---
### 4. Update the `Snakefile`
- Modify lines 1-29 to adjust paths and parameters.
- Ensure all necessary files for the Augur pipeline are present, including:
  - `sequences.fasta` & `metadata.tsv` 
    - can be downloaded from NCBI Virus via ingest: `FETCH_SEQUENCES==True`
  - [`auspice_config.json`](resources/auspice_config.json)
- These files are essential for building the reference tree and running Nextclade.
- Turning the parameter `STAR_ROOT=True` creates a star-like root after augur refine. This was developed for highly recombinant viruses. 
    - For the first run, run it with `STAR_ROOT=False`.
    - It needs the [clade_map](resources/clade_map.tsv) file: Download the auspice metadata (Download Data > METADATA (TSV)). Extract `accession` and `clade_membership` and save it as .tsv file with the headers ["accession","clade"].

---

## Runnning the `Snakefile`
To create the auspice JSON and a Nextclade example dataset:
```bash
snakemake --cores 9 all
```
This runs Nextclade on the example sequences in [`out-dataset/sequences.fasta`](out-dataset/sequences.fasta) using the dataset in `dataset`. The results are saved to the `test_out` directory and contain alignment, aligned translations and a summary TSV file.

## Visualizing the Nextclade build

One can also use the dataset in Nextclade Web by hosting the dataset through a local web server. For example, after having installed `node` and run `npm install -g serve`, one can host the dataset via:

```bash
serve --cors out-dataset -l 3000
```

And open Nextclade Web with a URL parameter `dataset-url` pointing to the local web server:

```bash
https://master.clades.nextstrain.org/?dataset-url=http://localhost:3000
```

Once the web page loads, you can click "Load example" and click run to test. You may want to reduce the maximum number of nucleotide markers to 500 to prevent Nextclade from freezing (click "Settings" at the top right, then select the "Sequence view" tab and reduce "Max. nucleotide markers to 500).


---
This guide provides a structured workflow for setting up Nextclade for Enterovirus A71. If you encounter issues, refer to the [official documentation](https://docs.nextstrain.org/projects/nextclade/en/stable/index.html#).