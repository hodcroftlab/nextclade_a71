"""
Star-like tree restructuring with recombinant clade handling.

Functional approach (no classes) for interactive testing in VS Code.
Each function is independently testable and returns intermediate objects.

Usage:
    python star_like_rooter.py \
        --input_tree "results/tree.nwk" \
        --input_clades "resources/clade_map.tsv" \
        --recombinant_accessions "resources/recombinants.tsv" \
        --output_tree "results/star_tree.nwk" \
        --ancestral_name "ancestral_sequence" \
        --strain_id_field "accession" \
        --recombinant_clades "C2.2 C2.3 C1.2 C1.3 E F A C6" \
        --root_name "NODE_0000000"
"""

from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree
import pandas as pd
import sys
import logging
import argparse
from typing import Dict, List, Optional, Set
import matplotlib.pyplot as plt

#%%
# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# ============================================================================
# Data loading
# ============================================================================


def load_recombinants(recombinants_file: Optional[str]) -> Set[str]:
    """
    Load recombinant accession IDs from TSV file.
    
    Args:
        recombinants_file: Path to recombinants TSV. If None, return empty set.
        
    Returns:
        Set of recombinant accession IDs (stripped, converted to str).
    """
    if recombinants_file is None:
        logger.info("No recombinants file provided; proceeding without recombinant accessions.")
        return set()    
    try:
        df = pd.read_csv(recombinants_file, sep="\t")
        accession_col = df.columns[0]
        return set(
            df[accession_col]
            .astype(str)
            .str.strip()
        )
    except Exception as e:
        logger.warning(f"Could not load recombinants file: {e}")
        return set()


# ============================================================================
# Group assignment
# ============================================================================

def assign_group(
    clade: str,
    strain_id: str,
    recombinant_accessions: Set[str],
    recombinant_clades: List[str],
) -> str:
    """
    Assign a strain to a group based on clade and recombinant status.
           
    Args:
        clade: Clade name from metadata.
        strain_id: Strain identifier.
        recombinant_accessions: Set of known recombinant accessions.
        recombinant_clades: List of recombinant clade names.
        
    Returns:
        Group name (str).
    """
    # Check if strain is in recombinant list
    strain_id = str(strain_id).strip()
    if strain_id in recombinant_accessions:
        return "RF"

    # Handle missing or invalid clade
    if pd.isna(clade):
        return "unassigned"

    clade = str(clade).strip()

    # Check if clade is explicitly recombinant
    if clade in recombinant_clades:
        return clade
    else:
        return "keep"

def is_recombinant_group(group: str) -> bool:
    """A strain counts as recombinant for re-rooting purposes if it was
    matched to a known recombinant accession or an explicitly recombinant
    clade. "keep" and "unassigned" strains are not recombinant."""
    return group not in ("keep", "unassigned")


def compute_strain_groups(
    df: pd.DataFrame,
    strain_id_field: str,
    recombinant_accessions: Set[str],
    recombinant_clades: List[str],
) -> Dict[str, str]:
    """
    Assign every strain in the metadata to a group (clade name, "RF",
    "keep", or "unassigned").

    Args:
        df: DataFrame with clade assignments.
        strain_id_field: Column name for strain IDs.
        recombinant_accessions: Set of recombinant accessions.
        recombinant_clades: List of recombinant clade names.

    Returns:
        Dictionary mapping strain ID → group name.
    """
    # Accept either a pre-split list (e.g. from argparse nargs="*") or a
    # single space-separated string.
    if isinstance(recombinant_clades, str):
        recombinant_clades = recombinant_clades.split()

    strain_to_group: Dict[str, str] = {}

    for _, row in df.iterrows():
        strain = str(row[strain_id_field]).strip()
        if not strain or strain == "nan":
            continue

        strain_to_group[strain] = assign_group(
            row["clade"],
            strain,
            recombinant_accessions,
            recombinant_clades,
        )

    return strain_to_group


# ============================================================================
# Recombinant ancestor finding
# ============================================================================

def find_recombinant_ancestors(
    tree: Tree,
    strain_to_group: Dict[str, str],
) -> List[Clade]:
    """
    For every recombinant-classified strain, climb from its tip up through
    ancestors for as long as the ancestor's *entire* set of descendant
    terminals is also recombinant-classified (regardless of which specific
    recombinant group they belong to). This finds the highest node in the
    tree that is still fully recombinant, so an isolated recombinant tip
    and a whole recombinant subtree (or a recombinant clade nested inside
    another recombinant clade) are all collapsed to a single ancestor.

    Args:
        tree: Tree to search (not modified).
        strain_to_group: Mapping of strain ID → group name.

    Returns:
        Deduplicated list of highest fully-recombinant ancestor clades.
        Never includes the tree root itself.
    """
    parent_of: Dict[int, Clade] = {}
    fully_recombinant: Dict[int, bool] = {}

    for clade in tree.find_clades(order="postorder"):
        for child in clade.clades:
            parent_of[id(child)] = clade

        if clade.is_terminal():
            fully_recombinant[id(clade)] = clade.name is not None and is_recombinant_group(
                strain_to_group.get(clade.name, "unassigned")
            )
        else:
            fully_recombinant[id(clade)] = all(
                fully_recombinant[id(child)] for child in clade.clades
            )

    root = tree.root
    top_ancestors: Dict[int, Clade] = {}

    for terminal in tree.get_terminals():
        if not fully_recombinant.get(id(terminal), False):
            continue

        node = terminal
        while True:
            parent = parent_of.get(id(node))
            if parent is None or parent is root or not fully_recombinant.get(id(parent), False):
                break
            node = parent

        top_ancestors[id(node)] = node

    return list(top_ancestors.values())


def name_for_ancestor(node: Clade, strain_to_group: Dict[str, str]) -> str:
    """Derive a human-readable name for a recombinant ancestor from the
    groups of the strains it subtends."""
    groups = sorted({
        strain_to_group.get(t.name, "unassigned") for t in node.get_terminals()
    })
    label = "_".join(groups) if len(groups) <= 3 else f"{groups[0]}_and_{len(groups) - 1}_more"
    return f"{label}_root"


# ============================================================================
# Distance calculation (for branch length preservation)
# ============================================================================

def distance_to_root(tree: Tree, clade: Clade) -> float:
    """
    Compute sum of branch lengths from tree root to given clade.
    
    This preserves the original root-to-node distance when a clade is
    reattached to a new root, maintaining phylogenetic signal.
    
    Args:
        tree: Tree containing the clade.
        clade: Target clade.
        
    Returns:
        Sum of branch lengths from root to clade. Returns 0.0 if path not found.
    """
    try:
        path = tree.get_path(clade)
    except Exception:
        return 0.0
    
    total = 0.0
    for p in path:
        branch_length = getattr(p, "branch_length", None)
        if branch_length is not None:
            total += float(branch_length)
    
    return total


# ============================================================================
# Tree manipulation
# ============================================================================

def move_clade_to_new_parent(
    tree: Tree,
    clade: Clade,
    new_parent: Clade,
) -> bool:
    """
    Move a clade from its current parent to a new parent in the tree.
    
    Steps:
        1. Compute distance from original tree root to clade.
        2. Set clade.branch_length = distance (preserve root-to-node distance).
        3. Remove clade from its current parent.
        4. If parent becomes childless, remove parent from grandparent.
        5. Attach clade to new_parent.
    
    Args:
        tree: Original tree (used to compute distances).
        clade: Clade to move.
        new_parent: Target parent clade.
        
    Returns:
        True if successful, False if error.
    """
    try:
        # Compute original distance from root to this clade
        orig_dist = distance_to_root(tree, clade)
        
        # Preserve branch length: set it to original root distance
        # so that when attached to new root, the distance is maintained
        if orig_dist > 0.0:
            clade.branch_length = orig_dist
        
        # Find and remove from current parent
        for parent in tree.find_clades():
            if clade in parent.clades:
                parent.clades.remove(clade)

                # If parent now has no children, remove it from grandparent
                if not parent.clades:
                    for grandparent in tree.find_clades():
                        if parent in grandparent.clades:
                            grandparent.clades.remove(parent)
                            break

                break

        # Attach to new parent
        new_parent.clades.append(clade)
        return True
    
    except Exception as e:
        logger.error(f"Error moving clade {clade.name}: {e}")
        return False


# ============================================================================
# Tree creation
# ============================================================================

def create_star_tree(
    tree: Tree,
    strain_to_group: Dict[str, str],
    root_name: str = "NODE_0000000",
) -> Tree:
    """
    Restructure `tree` in place so every fully-recombinant subtree is
    reattached directly to the tree's existing root, while the rest of the
    tree (the non-recombinant backbone, including reference/outgroup
    strains) keeps its original structure and position.

    This walks the tree itself: for each recombinant-classified strain it
    climbs to the highest ancestor whose entire descendant set is also
    recombinant-classified, so adjacent or nested recombinant clades
    collapse into a single reattached subtree instead of being processed
    (and potentially destroyed) independently. Recombinant clades that
    aren't phylogenetically adjacent to one another climb to distinct
    ancestors and are reattached to the root independently.

    Args:
        tree: Tree to restructure in place.
        strain_to_group: Mapping of strain ID → group name.
        root_name: Name to (re)assign to the tree's root node.

    Returns:
        The same tree, restructured.
    """
    top_ancestors = find_recombinant_ancestors(tree, strain_to_group)

    if not top_ancestors:
        logger.warning("No recombinant strains found; tree left unchanged.")

    root = tree.root
    root.name = root_name

    attached = 0
    total_strains_attached = 0

    for node in top_ancestors:
        if not node.name or node.name.startswith("NODE_"):
            node.name = name_for_ancestor(node, strain_to_group)

        n_terminals = node.count_terminals()

        if move_clade_to_new_parent(tree, node, root):
            attached += 1
            total_strains_attached += n_terminals
        else:
            logger.error(f"Failed to attach recombinant clade rooted at '{node.name}'")

    logger.info(
        f"Star-like restructuring: {attached} recombinant clades reattached to root, "
        f"{total_strains_attached} strains moved, backbone preserved"
    )

    return tree

#%%
# ============================================================================
# Main workflow
# ============================================================================

def main():
    """Parse arguments and run workflow."""
    parser = argparse.ArgumentParser(
        description="Create star-like tree with recombinant clade handling"
    )
    parser.add_argument("--input_tree", required=True, default="results/tree.nwk" ,help="Input tree (Newick)")
    parser.add_argument("--input_clades", required=True, default="resources/clade_map.tsv" ,help="Clade assignments (TSV)")
    parser.add_argument("--recombinant_accessions", default=None, help="Recombinant accessions (TSV)")
    parser.add_argument("--output_tree", required=True, default="star_tree.nwk", help="Output tree (Newick)")
    parser.add_argument("--strain_id_field",required=True, default="accession", help="Column name for strain IDs)")
    parser.add_argument("--recombinant_clades",nargs="*",default=[],help="Recombinant clade names")
    parser.add_argument("--root_name", default="NODE_0000000", help="Name for new root node")
    args = parser.parse_args()

    input_tree = args.input_tree
    input_clades = args.input_clades
    recombinants_file = args.recombinant_accessions
    output_tree = args.output_tree
    strain_id_field = args.strain_id_field
    recombinant_clades = args.recombinant_clades
    root_name = args.root_name

    # %%
    # input_tree = "results/tree.nwk" 
    # input_clades = "results/clades_metadata.tsv" 
    # recombinants_file = "resources/recombinants.tsv" 
    # output_tree = "results/star_tree.nwk" 
    # strain_id_field = "accession" 
    # recombinant_clades = "C2.2 C2.3 C1.2 C1.3 E F A C6" 
    # root_name = "NODE_0000000"


    
    try:
        #%%
        # Load data
        logger.info("Loading input files...")
        tree = Phylo.read(input_tree, "newick")
        df = pd.read_csv(input_clades, sep="\t")
        recombinant_accessions = load_recombinants(recombinants_file)
        #%%
        terminal_names = {t.name for t in tree.get_terminals()}
        # df: drop rows with NODES in accession column
        df = df[~df[strain_id_field].str.startswith("NODE_")]

        logger.info(
            f"Loaded tree with {len(terminal_names)} terminals "
            f"and {len(df)} clade assignments"
        )
        #%%
        # Assign groups and filter
        logger.info("Assigning groups...")
        strain_to_group = compute_strain_groups(
            df,
            strain_id_field,
            recombinant_accessions,
            recombinant_clades,
        )

        #%%
        # Create star tree
        logger.info("Creating star tree...")
        star_tree = create_star_tree(
            tree,
            strain_to_group,
            root_name=root_name,
        )
        #%%
        # Plot tree
        plt.figure(figsize=(12, 8))
        Phylo.draw(star_tree, do_show=False)
        plt.savefig("results/star_tree.png", dpi=300)
        #%%
        # Write output
        logger.info(f"Writing output to {output_tree}...")
        Phylo.write(star_tree, output_tree, "newick")
        logger.info("Successfully completed star tree creation")
        #%%
    except Exception as e:
        logger.error(f"Script failed: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
