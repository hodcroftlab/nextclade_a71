from Bio import Phylo
from Bio.Phylo.BaseTree import Clade
import pandas as pd
from collections import defaultdict
import ipdb
from copy import deepcopy
import matplotlib.pyplot as plt
import sys
from Bio.Phylo import draw

# Load inputs
tree = Phylo.read(snakemake.input.tree, "newick")
df = pd.read_csv(snakemake.input.clades, sep="\t")
id_field = snakemake.params.strain_id_field

# tree = Phylo.read("results/star_tree.nwk", "newick")
# df = pd.read_csv("resources/clade_map.tsv", sep="\t")
# id_field = "accession"

# Assign custom groups
def assign_group(clade):
    if clade.startswith("B"):
        return "B"
    elif clade == "C2.r":
        return "C2.r"
    elif clade == "C1-like":
        return "C1-like"
    elif clade.startswith("C"):
        return "C"
    # elif clade.startswith("A"):
    #     return "A"
    else:
        return "A/E/F"

df["group"] = df["clade"].apply(assign_group)

# Group tips by group
clade_to_strains = defaultdict(list)
for _, row in df.iterrows():
    clade_to_strains[row["group"]].append(str(row[id_field]))

# Create a new root
new_root = Clade(name="star_root")
new_tree = Phylo.BaseTree.Tree(root=new_root)

# Move MRCA nodes from the original tree to the new root
def move_clade(tree, clade, new_parent):
    for parent in tree.find_clades():
        if clade in parent.clades:
            parent.clades.remove(clade)
            break
    new_parent.clades.append(clade)

for group, strains in clade_to_strains.items():
    try:
        mrca = tree.common_ancestor(strains)
        mrca.name = group
        move_clade(tree, mrca, new_root)
    except:
        print(f"No MRCA found for group {group}", file=sys.stderr)

# ipdb.set_trace()
## Plot the tree with colors
## Build tip → group mapping
# tip_colors = {}
# group_to_color = {}
# palette = plt.cm.tab10.colors
# for i, group in enumerate(df["group"].unique()):
#     group_to_color[group] = palette[i % len(palette)]

# for _, row in df.iterrows():
#     tip_colors[row[id_field]] = group_to_color[row["group"]]

# # Map colors directly using label_func and label_colors
# label_colors = {tip: color for tip, color in tip_colors.items()}

# fig = plt.figure(figsize=(10, 20))
# ax = fig.add_subplot(1, 1, 1)

# Phylo.draw(
#     tree,
#     label_func=lambda group: group.name,
#     do_show=False,
#     axes=ax,
#     label_colors=label_colors
# )

# plt.tight_layout()
# plt.show()

# Write output
Phylo.write(new_tree, snakemake.output.tree, "newick")
