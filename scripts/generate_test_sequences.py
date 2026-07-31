# scripts/generate_test_sequences.py
"""
Generate all test sequences: fragments, recombinants, non-EV species controls.
"""
import random
import argparse
from pathlib import Path
from Bio import SeqIO, Entrez
import pandas as pd
import time
import sys
import os

def fetch_sequences_from_entrez(taxid, virus_name, num_seqs=5, email=None, api_key=None, min_length=7000):
    """
    Fetch random sequences from NCBI Entrez for a given taxon ID.
    Excludes sequences matching the dataset virus name.
    Only fetches complete, full-length genomes.
    
    Args:
        taxid: NCBI taxonomy ID for the species group
        virus_name: Name of the dataset virus to exclude (e.g., "EV-D68", "EV-D71")
        num_seqs: Number of sequences to fetch (default 5)
        email: Email for NCBI Entrez access
        api_key: Optional NCBI API key for faster access
        min_length: Minimum sequence length in bp (default 7000)
    
    Returns:
        List of (accession, virus_short, sequence) tuples
    """
    if email:
        Entrez.email = email
    else:
        from dotenv import load_dotenv, find_dotenv
        import os
        
        load_dotenv(find_dotenv())
        Entrez.email = os.environ.get("EMAIL")

    if api_key:
        Entrez.api_key = api_key
    
    print(f"Fetching {num_seqs} complete sequences from taxid {taxid}, excluding {virus_name}...", file=sys.stderr)
    
    # Search for complete genomes in the taxon
    search_term = (
        f'txid{taxid}[Organism] AND '
        f'(complete genome[Title] OR complete sequence[Title]) AND '
        f'biomol_genomic[PROP]'
    )
    
    try:
        # First, find out how many records match so we don't silently truncate
        # the pool: esearch's default ordering is newest-submission-first, so a
        # small retmax can bury rarer/older serotypes behind hundreds of recent
        # accessions and they'd never be sampled at all.
        count_handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=0)
        total_count = int(Entrez.read(count_handle)["Count"])
        count_handle.close()

        handle = Entrez.esearch(
            db="nucleotide",
            term=search_term,
            retmax=max(total_count, 1),  # fetch the full matching set, not just the newest 100
            idtype="acc"
        )
        record = Entrez.read(handle)
        handle.close()

        if not record["IdList"]:
            print(f"No complete sequences found for taxid {taxid}", file=sys.stderr)
            return []
        
        all_accs = record["IdList"]
        print(f"Found {len(all_accs)} complete sequences, fetching details...", file=sys.stderr)
        
        fetched_seqs = []
        attempted = 0
        seen_viruses = set()  # Track viruses we've already fetched
        max_attempts = len(all_accs)
        
        while len(fetched_seqs) < num_seqs and attempted < max_attempts:
            # Pick random accession
            if not all_accs:
                break
            
            acc = all_accs.pop(random.randint(0, len(all_accs) - 1))
            attempted += 1
            
            try:
                # Fetch GenBank record
                handle = Entrez.efetch(db="nucleotide", id=acc, rettype="gb", retmode="text")
                seq_record = SeqIO.read(handle, "genbank")
                handle.close()
               
                # Check sequence length
                seq_len = len(seq_record.seq)
                if seq_len < min_length:
                    continue
                
                # Extract virus info from description or organism
                organism = seq_record.annotations.get("organism", "")
                description = seq_record.description.upper()
                virus_short = seq_record.annotations.get("source", organism)
                # source is often formatted like "enterovirus D68 (EV-D68)" - prefer the
                # parenthetical short name when present, otherwise use it as-is (many
                # records, e.g. generic "Enterovirus D", have no parenthetical at all)
                if "(" in virus_short and ")" in virus_short:
                    virus_short = virus_short.split("(")[1].split(")")[0]

                # Skip if we've already fetched this virus species
                if organism in seen_viruses:
                    continue
                
                # Verify it's complete
                if "COMPLETE" not in description:
                    continue
                
                # Skip if it matches our dataset virus
                if virus_name.lower() in organism.lower() or virus_name.lower() in description.lower():
                    continue
                               
                fetched_seqs.append((acc, virus_short, str(seq_record.seq)))
                seen_viruses.add(organism)
                print(f"✓ {acc} | {virus_short} ({seq_len} bp, complete)", file=sys.stderr)
                
                # Rate limiting
                time.sleep(0.35)  # 3 req/sec without API key
                
            except Exception as e:
                continue
        
        print(f"Successfully fetched {len(fetched_seqs)}/{num_seqs} sequences", file=sys.stderr)
        return fetched_seqs
        
    except Exception as e:
        print(f"Error searching NCBI: {e}", file=sys.stderr)
        return []


def write_entrez_sequences(sequences, output_file):
    """Write fetched sequences to FASTA with formatted headers."""
    with open(output_file, "w") as f:
        for accession, virus_short, seq in sequences:
            header = f"{accession} | {virus_short}"
            f.write(f">{header}\n{seq}\n")
    return len(sequences)

def eligible_seqs(records, min_length):
    """Filter records by minimum length."""
    return [r for r in records if len(r) >= min_length]

def parse_gff_gene_coords(gff_file):
    """Parse gene/CDS coordinates from a GFF3 file into 0-based half-open (start, end) tuples."""
    gene_coords = {}
    if not gff_file or not gff_file.exists():
        return gene_coords
    with open(gff_file) as gff:
        for line in gff:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] not in ("gene", "CDS"):
                continue
            start, end, attributes = parts[3], parts[4], parts[8]
            attrs = dict(item.split("=", 1) for item in attributes.split(";") if "=" in item)
            gene_name = attrs.get("Name") or attrs.get("gene") or attrs.get("ID")
            if gene_name:
                gene_coords[gene_name] = (int(start) - 1, int(end))
    return gene_coords


def random_ungapped_fragment(seq, length):
    """Return a random subsequence of `length` non-gap/N bases, or None if `seq` is too short."""
    ungapped_len = len(seq) - seq.count("-") - seq.count("N")
    if ungapped_len <= length:
        return None
    start = random.randint(0, ungapped_len - length)
    return seq[start:start + length]


def generate_fragments(sequences_file, rivm_file, clades_file, nextstrain_file ,output_file,
                       lengths=None, genes=None, gff_file=None):
    """Generate fragment test sequences."""
    if lengths is None:
        lengths = range(100, 3100, 100)
    if genes is None:
        genes = ["VP1", "3D"]

    gene_coords = parse_gff_gene_coords(gff_file)

    # Read all sequences from the input file
    records = list(SeqIO.parse(sequences_file, "fasta"))

    # Read Nextstrain Clades
    nextstrain_df = pd.read_csv(nextstrain_file, sep="\t")
    meta_ids = set(nextstrain_df["accession"] if "accession" in nextstrain_df.columns else nextstrain_df.iloc[:, 0])
    records = [r for r in records if r.id in meta_ids]

    ns_map = nextstrain_df.set_index("accession")["clade_membership"].to_dict()

    # filter records NOT in clades file
    acc = pd.read_csv(clades_file, sep="\t").accession
    records = [r for r in records if r.id not in acc]

    rivm_map = {}
    if rivm_file:
        rivm = pd.read_csv(rivm_file, sep=",")
        rivm.dropna(axis=0, inplace=True)
        if "accession" in rivm:
            rivm_map = rivm.set_index("accession")["RIVM"].to_dict()
        elif "name" in rivm:
            rivm_map = rivm.set_index("name")["VP1 subgenogroup"].to_dict()
        else: print("neither accession nor name in ", rivm_file,". Please check!")

    # Most records are already partial (e.g. VP1-only), so only keep, per gene,
    # the records long enough to actually contain that gene's coordinates.
    gene_records = {
        gene: [r for r in records if len(r.seq) >= end]
        for gene, (start, end) in gene_coords.items()
    }

    with open(output_file, "w") as out_handle:
        for length in lengths:
            for gene in genes:
                if gene not in gene_coords:
                    print(f"Gene {gene} not found in GFF3 annotation, skipping.", file=sys.stderr)
                    continue
                eligible = gene_records[gene]
                if not eligible:
                    continue
                record = random.choice(eligible)
                start, end = gene_coords[gene]
                fragment = random_ungapped_fragment(record.seq[start:end], length)
                if fragment is not None:
                    cl = rivm_map.get(record.id, "NA")
                    ns = ns_map.get(record.id, "NA")
                    header = f"{record.id}_partial_{length}_{gene}|{ns}|{cl}"
                    out_handle.write(f">{header}\n{fragment}\n")

            record = random.choice(records)
            seq_len = len(record.seq)
            while seq_len < length:
                record = random.choice(records)
                seq_len = len(record.seq)
            cl = rivm_map.get(record.id, "NA")
            ns = ns_map.get(record.id, "NA")
            start = random.randint(0, seq_len - length)
            fragment_seq = record.seq[start:start+length]
            header = f"{record.id}_partial_{length}|{ns}|{cl}"
            out_handle.write(f">{header}\n{fragment_seq}\n")

    return sum(1 for line in open(output_file) if line.startswith(">"))

def generate_recombinants(sequences_file, metadata, clades_file, ev_file, 
                          output_file, inter_count=10, intra_count=10, min_length=3500):
    """Generate inter- and intra-recombinant test sequences."""
    seqs = list(SeqIO.parse(sequences_file, "fasta"))
    ns_ids = set(pd.read_csv(metadata, sep="\t").accession)
    seqs = eligible_seqs([r for r in seqs if r.id in ns_ids], min_length)
    
    clade_map = pd.read_csv(clades_file, sep="\t").set_index("accession")["clade"].to_dict()
    clade2seqs = {}
    for r in seqs:
        clade = clade_map.get(r.id, "NA")
        if clade != "NA":
            clade2seqs.setdefault(clade, []).append(r)
    
    clades = [c for c in clade2seqs if len(clade2seqs[c]) > 0]
    ev_seqs = eligible_seqs(list(SeqIO.parse(ev_file, "fasta")), min_length)
    
    recombinants = []
    
    # Intra-typic recombinants (between clades)
    for i in range(intra_count):
        if len(clades) < 2:
            break
        c1, c2 = random.sample(clades, 2)
        p1, p2 = random.choice(clade2seqs[c1]), random.choice(clade2seqs[c2])
        minlen = min(len(p1.seq), len(p2.seq))
        if minlen < min_length:
            continue
        breakpoint = random.randint(1, minlen - 1)
        header = f"intra_{p1.id}_{c1}_{breakpoint}_{p2.id}_{c2}"
        recomb_seq = str(p1.seq[:breakpoint]) + str(p2.seq[breakpoint:])
        recombinants.append((header, recomb_seq))
    
    # Inter-typic recombinants (target virus x related species)
    for j in range(inter_count):
        p1 = random.choice(seqs)
        p2 = random.choice(ev_seqs)
        minlen = min(len(p1.seq), len(p2.seq))
        if minlen < min_length:
            continue
        breakpoint = random.randint(1, minlen - 1)
        header = f"inter_{p1.id}_ virus_{breakpoint}_{p2.id}_{p2.description.split('|')[1].strip()}"
        recomb_seq = str(p1.seq[:breakpoint]) + str(p2.seq[breakpoint:])
        recombinants.append((header, recomb_seq))
    
    # Write recombinants
    with open(output_file, "w") as f:
        for header, seq in recombinants:
            f.write(f">{header}\n{seq}\n")
    
    return len(recombinants)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate test sequences")
    parser.add_argument("--sequences", required=True, help="Input sequences file")
    parser.add_argument("--nextstrain", required=True, help="Nextstrain metadata TSV")
    parser.add_argument("--clades", required=True, help="Clades metadata TSV")
    parser.add_argument("--rivm", required=False, help = "RIVM Clade Assignment")
    parser.add_argument("--output-fragments", required=True, help="Output fragments FASTA")
    parser.add_argument("--output-recombinants", required=True, help="Output recombinants FASTA")
    parser.add_argument("--output-ev", required=False, help="Output EV species sequences FASTA (if fetching from NCBI)")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--ev", required=False, help="EV species sequences file (skip Entrez if provided)")
    parser.add_argument("--taxid", type=str, required=False, help="NCBI Taxon ID for EV species (if not provided as file)")
    parser.add_argument("--virus", type=str, required=False, default="EV", help="Name of dataset virus to exclude from Entrez fetch")
    parser.add_argument("--genes", required=False, nargs='+', default=None, help="Genes for which fragments are created")
    parser.add_argument("--email", required=False, help="Email for NCBI Entrez access")
    
    args = parser.parse_args()
    random.seed(args.seed)
    
    Path(args.output_fragments).parent.mkdir(parents=True, exist_ok=True)
    Path(args.output_recombinants).parent.mkdir(parents=True, exist_ok=True)

    # Handle EV species sequences: either use provided file or fetch from NCBI
    if args.ev:
        ev_file = args.ev
        print(f"Using provided EV species file: {ev_file}", file=sys.stderr)
    elif args.taxid:
        ev_file = args.output_ev or Path(args.output_fragments).parent / "EV_fetched.fasta"
        fetched = fetch_sequences_from_entrez(
            args.taxid,
            args.virus,
            num_seqs=10,
            email=args.email,
            min_length=7000
        )
        
        if fetched:
            num_ev = write_entrez_sequences(fetched, ev_file)
            print(f"Fetched and saved {num_ev} EV species sequences to {ev_file}", file=sys.stderr)
        else:
            print("Failed to fetch EV species sequences, skipping recombinant generation", file=sys.stderr)
            ev_file = None
    else:
        print("Error: Either --ev file or --taxid must be provided", file=sys.stderr)
        sys.exit(1)

    gff_file = Path("dataset/genome_annotation.gff3")

    if args.rivm:
        rivm = Path(args.rivm)
    else:
        rivm = None

    n_frags = generate_fragments(args.sequences, rivm, args.clades,args.nextstrain, args.output_fragments,
                                  genes=args.genes, gff_file=gff_file)
    n_recomb = generate_recombinants(args.sequences, args.nextstrain, args.clades, ev_file, args.output_recombinants)
    
    print(f"Generated {n_frags} fragments and {n_recomb} recombinants")