# Nextclade Dataset for "Enterovirus A71" based on reference "U22521.1"


| Key              | Value                                                                 |
|------------------|-----------------------------------------------------------------------|
| authors          | Nadia Neuner-Jehle, Emma B. Hodcroft                                  |
| name             | Enterovirus A71                                                       |
| reference        | U22521.1                                                              |
| dataset path     | ...                                     |

The dataset represents the species *Enterovirus A71* (EV-A71), a major cause of hand, foot, and mouth disease and neurological complications. It is based on the complete genome of the BrCr strain (U22521.1), one of the first EV-A71 isolates sequenced. The dataset supports genogroup classification, mutation QC, and phylogenetic placement.

The example sequences include a representative global subsample across known genogroups A–G.

## Scope of this dataset

This dataset is designed for analyzing full genomes or full VP1 sequences of *Enterovirus A71*. Other enterovirus types are not supported and may be flagged as outliers or "Outgroup."

## Features

This dataset was created from publicly available EV-A71 sequences using a phylogenetic workflow adapted for enteroviruses. Clade and genogroup annotations follow literature references (e.g., van der Sanden et al., 2009) and community naming conventions.

Quality control is customized to EV-A71 features:
- `frameShifts` strictly penalized (`scoreWeight: 100`)
- `divergence.maxDivergence` set to 0.15
- Suitable for full genomes or full VP1 gene sequences

## What is a Nextclade dataset?

Read more about Nextclade datasets in the [Nextclade documentation](https://docs.nextstrain.org/projects/nextclade/en/stable/user/datasets.html).