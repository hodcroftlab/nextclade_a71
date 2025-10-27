from Bio import Phylo
from Bio.Phylo.BaseTree import Clade
import pandas as pd
from collections import defaultdict
import sys
import logging
import argparse
from typing import Dict, List, Optional
import ipdb

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class CladeGroupAssigner:
    """Class to handle clade group assignment and tree restructuring."""
    
    def __init__(self, tree_file: str, clades_file: str, strain_id_field: str, 
                 recombinant_clades: Optional[List[str]] = None, 
                 recombinants_file: Optional[str] = None):
        self.tree = Phylo.read(tree_file, "newick")
        self.df = pd.read_csv(clades_file, sep="\t")
        self.strain_id_field = strain_id_field
        self.recombinant_clades = recombinant_clades or []
        
        # Load recombinant accessions
        self.recombinant_accessions = set()
        if recombinants_file:
            try:
                rec_df = pd.read_csv(recombinants_file, sep="\t")
                accession_col = rec_df.columns[0]
                self.recombinant_accessions = set(rec_df[accession_col].astype(str).str.strip())
                # logger.info(f"Loaded {len(self.recombinant_accessions)} recombinant accessions")
            except Exception as e:
                logger.warning(f"Could not load recombinants file: {e}")
        
        # Cache terminal names for faster lookup
        self.terminal_names = {tip.name for tip in self.tree.get_terminals()}
        # logger.info(f"Loaded tree with {len(self.terminal_names)} terminals and clade mapping with {len(self.df)} entries")
        
    def assign_group(self, clade: str, strain_id: str) -> str:
        """
        Assign custom groups based on clade names.
        
        Args:
            clade: The clade name to classify
            strain_id: The strain ID to check against recombinant list
            
        Returns:
            The assigned group name
        """
        # Check if in recombinant list first
        if strain_id in self.recombinant_accessions:
            return "recombinant"
        
        if pd.isna(clade) or not isinstance(clade, str):
            return "Unassigned"
            
        clade = clade.strip()
        
        # Check if clade is recombinant
        if clade in self.recombinant_clades:
            return clade
        elif clade.startswith("B"):
            return "B"
        elif clade.startswith("C"):
            return "C"
        else:
            return "recombinant"
    
    def get_clade_to_strains_mapping(self) -> Dict[str, List[str]]:
        """
        Create mapping from groups to strain lists.
        
        Returns:
            Dictionary mapping group names to lists of strain IDs
        """
        # Apply group assignment (modified to pass strain_id)
        self.df["group"] = self.df.apply(
            lambda row: self.assign_group(row["clade"], str(row[self.strain_id_field]).strip()),
            axis=1
        )
        
        # Group strains by assigned group
        clade_to_strains = defaultdict(list)
        
        for _, row in self.df.iterrows():
            strain_id = str(row[self.strain_id_field]).strip()
            group = row["group"]
            
            if strain_id and strain_id != "nan":  # Skip empty or NaN values
                clade_to_strains[group].append(strain_id)
        
        # Log group statistics
        for group, strains in clade_to_strains.items():
            valid_strains = [s for s in strains if s in self.terminal_names]
            # logger.info(f"Group '{group}': {len(strains)} strains, {len(valid_strains)} valid in tree")
            
        return dict(clade_to_strains)
    
    def find_mrca_safely(self, strains: List[str]) -> Optional[Clade]:
        """
        Safely find MRCA for a list of strains.
        
        Args:
            strains: List of strain names
            
        Returns:
            MRCA clade or None if not found
        """
        if len(strains) == 1:
            # For single strain, find the terminal node
            for terminal in self.tree.get_terminals():
                if terminal.name == strains[0]:
                    return terminal
            return None
        
        try:
            return self.tree.common_ancestor(strains)
        except Exception as e:
            logger.warning(f"Could not find MRCA for strains {strains[:5]}{'...' if len(strains) > 5 else ''}: {e}")
            return None
    
    def distance_to_root(self, clade: Clade) -> float:
        """
        Compute sum of branch lengths from the current tree root to the given clade.
        This is used to preserve original root-to-node distance when reattaching clades
        to a new root.
        """
        try:
            path = self.tree.get_path(clade)  # path is list of clades from root->...->clade (excluding root)
        except Exception:
            # if get_path fails, return 0
            return 0.0
        total = 0.0
        for p in path:
            if getattr(p, "branch_length", None):
                total += float(p.branch_length)
        return total

    def move_clade_safely(self, clade: Clade, new_parent: Clade, problematic_parents: list) -> bool:
        """
        Safely move a clade to a new parent.
        
        Args:
            clade: The clade to move
            new_parent: The new parent clade
            problematic_parents: A list of problematic parent node names to check and remove if they have no children.
            
        Returns:
            True if successful, False otherwise
        """
        try:
            # compute original distance from root to this clade (before we remove it)
            orig_dist = self.distance_to_root(clade)
            # If orig_dist is zero or missing, keep whatever branch_length clade had
            if orig_dist > 0.0:
                # set branch_length on the clade so that when attached directly to new root
                # the root->clade distance equals the original root->clade distance.
                clade.branch_length = orig_dist

            # Find and remove from current parent
            for parent in self.tree.find_clades():
                if clade in parent.clades:
                    parent.clades.remove(clade)
                    
                    # Add to problematic_parents if not already present
                    if parent.name not in problematic_parents:
                        problematic_parents.append(parent.name)

                    # Check if the parent has no remaining children
                    if not parent.clades:
                        # Remove the parent from its grandparent
                        for grandparent in self.tree.find_clades():
                            if parent in grandparent.clades:
                                grandparent.clades.remove(parent)
                                break
                    break

            # Add the clade to the new parent
            new_parent.clades.append(clade)
            return True
    
        except Exception as e:
            logger.error(f"Error moving clade {clade.name}: {e}")
            return False
    
    def clean_tree_for_augur(self, tree: Phylo.BaseTree.Tree) -> Phylo.BaseTree.Tree:
        """
        Clean the tree to ensure compatibility with augur ancestral.
        Remove any nodes that don't have sequence data.
        
        Args:
            tree: The tree to clean
            
        Returns:
            Cleaned tree
        """
        # Get all terminal names from the original alignment/tree
        valid_terminals = self.terminal_names
        
        # Find and remove any terminals that shouldn't be there
        terminals_to_remove = []
        for terminal in tree.get_terminals():
            if terminal.name not in valid_terminals:
                logger.warning(f"Removing problematic node: {terminal.name}")
                terminals_to_remove.append(terminal)

        # Remove problematic terminals
        for terminal in terminals_to_remove:
            # Find parent and remove
            for parent in tree.find_clades():
                if terminal in parent.clades:
                    parent.clades.remove(terminal)
                    break
                    
        return tree

    def create_star_tree(self, root_name: str = "NODE_0000000") -> Phylo.BaseTree.Tree:
        """Create a star-like tree with major clades attached to root."""
        clade_to_strains = self.get_clade_to_strains_mapping()

        # Create new root
        new_root = Clade(name=root_name)
        new_tree = Phylo.BaseTree.Tree(root=new_root)
        
        successful_groups = 0
        total_strains_processed = 0
        problematic_parents= []
        
        for group, strains in clade_to_strains.items():
            valid_strains = [strain for strain in strains if strain in self.terminal_names]
            
            if not valid_strains:
                logger.warning(f"No valid strains found for group '{group}'")
                continue
            
            mrca = self.find_mrca_safely(valid_strains)
            if mrca is None:
                logger.warning(f"Could not find MRCA for group '{group}'")
                continue
            
            # Ensure the clade has a proper name
            if not mrca.name or mrca.name.startswith("NODE_"):
                mrca.name = f"{group}_root"
            
            if self.move_clade_safely(mrca, new_root, problematic_parents):
                successful_groups += 1
                total_strains_processed += len(valid_strains)
        
        # Clean the tree before returning
        cleaned_tree = self.clean_tree_for_augur(new_tree)
        
        # logger.info(f"Successfully processed {successful_groups}/{len(clade_to_strains)} groups")
        return cleaned_tree

    
def main():
    """Main execution function with argument parsing."""
    parser = argparse.ArgumentParser(description="Create star-like tree with recombinant handling")
    parser.add_argument("--input_tree", required=True, help="Input tree file (Newick format)")
    parser.add_argument("--input_clades", required=True, help="Input clades TSV file")
    parser.add_argument("--recombinant_accessions", required=False, help="Recombinant accessions TSV file")
    parser.add_argument("--output_tree", required=True, help="Output tree file (Newick format)")
    parser.add_argument("--strain_id_field", required=True, help="Column name for strain IDs")
    parser.add_argument("--recombinant_clades", nargs="*", default=[], help="List of recombinant clade names")
    parser.add_argument("--root_name", default="NODE_0000000", help="Name for root node")
    
    args = parser.parse_args()
    
    try:
        # Initialize assigner with arguments
        assigner = CladeGroupAssigner(
            tree_file=args.input_tree,
            clades_file=args.input_clades,
            strain_id_field=args.strain_id_field,
            recombinant_clades=args.recombinant_clades,
            recombinants_file=args.recombinant_accessions
        )
        
        # Create star tree
        star_tree = assigner.create_star_tree(root_name=args.root_name)

        Phylo.draw(star_tree)
        
        # Write output
        Phylo.write(star_tree, args.output_tree, "newick")
        logger.info(f"Successfully wrote star tree to {args.output_tree}")
        
    except Exception as e:
        logger.error(f"Script failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()