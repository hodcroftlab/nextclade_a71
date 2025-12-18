# Inferred Ancestral Root Sub-Workflow

This sub-workflow generates a **Static Inferred Ancestor** for Enterovirus A71 (EV-A71) using outgroup rooting. This ancestral sequence represents the MRCA (most recent common ancestor) of all EV-A71 sequences and serves as a stable, biologically accurate reference for the Nextclade dataset.

> [!NOTE]  
> This README is for users who want to **regenerate** the inferred root (i.e., `INFERRENCE_RERUN = True`). For general information about how the inferred root is used in the main workflow, see the [main README](../README.md#inferred-ancestral-root-with-outgroup-rooting-recommended).

---

## When to Regenerate the Inferred Root

You should regenerate the inferred root when:
- **First-time setup:** `resources/inferred-root.fasta` doesn't exist yet
- **New data:** You've added significant new sequences that may shift the root position
- **Updated outgroups:** You've changed which enterovirus species are used as outgroups
- **Different subsampling:** You want to test sensitivity with a different sequence subset

---

## Configuration

### Snakefile Parameters

Edit the first 8 lines of [`inferred-root/Snakefile`](Snakefile) to configure:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `REFERENCE_ACCESSION` | `"U22521"` | GenBank accession for the reference (Fermon) |
| `MIN_DATE` | `"1990-01-01"` | Earliest collection date to include |
| `MIN_LENGTH` | `"7000"` | Minimum sequence length (bp) |
| `MAX_SEQS` | `"1000"` | Maximum sequences to subsample for tree building |
| `ID_FIELD` | `"accession"` | Metadata column for sequence IDs |
| `OUTGROUP` | List of accessions | Enterovirus species used to root the tree |

### Outgroup Sequences

Outgroup sequences must be placed in `resources/outgroup/` as individual FASTA files, named by their accession ID:

```
resources/outgroup/
├── AF326750.2.fasta
├── AY421769.1.fasta
```

**Current default outgroups:**
- Other enterovirus D species (e.g., EV-D70, EV-D94, EV-D111)
- Closely related enteroviruses that root the EV-A71 tree

To add or modify outgroups:
1. Add FASTA files to `resources/outgroup/`
2. Update the `OUTGROUP` list in `inferred-root/Snakefile`

---

## How It Works

The workflow performs the following steps:

### 1. Subsample & Filter
- Filters sequences by date (`MIN_DATE`) and length (`MIN_LENGTH`)
- Subsamples to `MAX_SEQS` sequences using `augur filter`
- The rule `join_fastas` merges the filtered sequences with the outgroup sequences 

### 2. Align
- Aligns the combined ingroup + outgroup sequences using MAFFT
- Produces a multiple sequence alignment for phylogenetic analysis

### 3. Build Tree
- Constructs a maximum-likelihood tree with `augur tree` (IQ-TREE)
- Tree includes both EV-A71 sequences and outgroup species

### 4. Root & Extract Ancestor
- [`pick_ancestral_sequence.py`](../scripts/pick_ancestral_sequence.py) reroots the tree on the outgroup(s)
- Identifies the **MRCA of the ingroup** (all EV-A71 sequences)
- Extracts the reconstructed ancestral sequence at this node
- **Fills gaps** with nucleotides from the reference to ensure a complete genome

### 5. Export
- Saves the inferred root as `resources/inferred-root.fasta`
- Creates a visualization of the rooted tree: `results/nwk_tree_outgroup.png`

---

## Running the Sub-Workflow

### Standalone Execution

To run the inferred-root sub-workflow independently:

```bash
cd inferred-root
snakemake --cores 9 all_sub
```

---

## Requirements

- **Snakemake** — Workflow management
- **Augur** (nextstrain/augur) — Phylogenetic tools (`filter`, `tree`, `merge`)
- **MAFFT** — Multiple sequence alignment
- **IQ-TREE** (via augur) — Maximum-likelihood tree inference
- **TreeTime** — Ancestral reconstruction (via custom script)
- **Python packages:** Biopython, pandas, typer, matplotlib

---

## Validation & Quality Control

After running, verify the output:

### 1. Check the Inferred Root Sequence

```bash
# View the sequence header
grep ">" ../resources/inferred-root.fasta

# Check sequence length
seqkit stats ../resources/inferred-root.fasta
```

Expected: Should be ~7,500 bp (full EV-A71 genome) with ID `ancestral_sequence`

### 2. Inspect the Rooted Tree

Open `results/nwk_tree_outgroup.png` to visually confirm:
- Outgroup sequences are at the base of the tree
- EV-A71 sequences form a monophyletic clade
- The root is positioned correctly between outgroup and ingroup

### 3. Verify Metadata Consistency

```bash
# Check that the ancestral sequence ID matches metadata
cat ../resources/static_inferred_metadata.tsv
```

The `accession` (or `strain`) column should contain `ancestral_sequence` matching the FASTA header.

---

## Troubleshooting

### Common Issues

**Problem:** `augur filter` produces no sequences  
**Solution:** Check that `MIN_DATE`, `MIN_LENGTH`, and `MAX_SEQS` parameters are appropriate for your dataset. Verify that `data/sequences.fasta` and `data/metadata.tsv` exist.

**Problem:** MAFFT alignment fails  
**Solution:** Ensure outgroup sequences are valid FASTA files and compatible with EV-A71 (same genomic region). Outgroups should be closely related enteroviruses.

**Problem:** Tree rooting fails with "outgroup not found"  
**Solution:** Verify that the outgroup accessions in the `OUTGROUP` list exactly match the FASTA headers in `resources/outgroup/*.fasta`.

**Problem:** Inferred sequence has many gaps  
**Solution:** This is expected for divergent regions. The script fills gaps with reference nucleotides. If excessive, consider exluding divergent sequences in `../resources/exclude.txt`.

**Problem:** `pick_ancestral_sequence.py` fails  
**Solution:** Ensure TreeTime is installed and that BioPython is up-to-date.

---

## Tips and Best Practices

- **Keep metadata synchronized:** Ensure `resources/static_inferred_metadata.tsv` has all required columns matching your main metadata
- **Document outgroup choices:** Add comments in the Snakefile explaining why specific outgroups were selected
- **Test sensitivity:** Try different `MAX_SEQS` values to confirm the inferred root is stable
- **Version control:** Commit the generated `resources/inferred-root.fasta` so others can reproduce your dataset without re-running inference
- **Regenerate periodically:** When major new sequences are added (e.g., new clades, geographic regions), consider regenerating the root

---

## Additional Resources

- **Template for other pathogens:** [dataset-template-inferred-root](https://github.com/enterovirus-phylo/dataset-template-inferred-root)
- **Augur documentation:** https://docs.nextstrain.org/projects/augur/
- **TreeTime documentation:** https://treetime.readthedocs.io/

---

## Author & Contact

- Maintainers: Nadia Neuner-Jehle, Alejandra González-Sánchez and Emma B. Hodcroft ([eve-lab.org](https://eve-lab.org/))
- For questions or suggestions, please [open an issue](https://github.com/enterovirus-phylo/dataset-template-inferred-root/issues) or email: eve-group[at]swisstph.ch