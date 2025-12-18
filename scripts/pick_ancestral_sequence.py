import typer

def main(tree_path: str, alignment: str, sequences: str, outgroup:str, reference:str, output: str):

    from Bio import Phylo, Seq, SeqIO
    from treetime import TreeAnc
    import matplotlib.pyplot as plt

    tree = Phylo.read(tree_path, "newick")
    tree.root_with_outgroup(outgroup.split(" "))

    tt = TreeAnc(tree, alignment, fill_overhangs=False)

    tree_dict = tt.get_tree_dict()
    td = {}
    for sequence in tree_dict:
        td[sequence.id] = sequence.seq

    sequences = SeqIO.to_dict(SeqIO.parse(sequences, "fasta"))
    ingroup_terminals = []
    for terminal in tree.get_terminals():
        if terminal.name in sequences.keys() and terminal.name not in outgroup and terminal.name not in reference:
           ingroup_terminals.append(terminal)

    ancestor_name = tt.tree.common_ancestor(ingroup_terminals).name
    
    print("\n Reference: ", reference.split(".")[0], 
    "\n Ancestor: ", ancestor_name)

    anc = str(td[ancestor_name])
    coord = str(td[reference.split(".")[0]])

    output_list = []
    for (ref,coord) in zip(anc,coord):
        if coord != "-":
            if ref == "-":
                output_list.append(coord)
            else:
                output_list.append(ref)
    output_seq = "".join(output_list)
  
    SeqIO.write(SeqIO.SeqRecord(Seq.Seq(output_seq), id="ancestral_sequence", description=""), output, "fasta")

    fig = plt.figure(figsize=(8, 10))     # define figure size
    ax = fig.add_subplot(1, 1, 1)
    Phylo.draw(tree, do_show=False, axes=ax)
    plt.savefig("results/nwk_tree_outgroup.png", dpi=300, bbox_inches="tight")
    plt.close(fig) 

if __name__ == "__main__":
    typer.run(main)