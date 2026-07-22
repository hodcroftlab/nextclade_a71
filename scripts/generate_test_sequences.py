# scripts/generate_test_sequences.py
"""
Generate all test sequences: fragments, recombinants, non-EV-A controls.
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
        virus_name: Name of the dataset virus to exclude (e.g., "CVA16", "EV-A71")
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
        handle = Entrez.esearch(
            db="nucleotide",
            term=search_term,
            retmax=100,  # Get more results to filter
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
                virus_short = seq_record.annotations.get("source",organism)
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

def generate_fragments(sequences_file, metadata_file, clades_file, output_file, 
                       lengths=None, genes=None):
    """Generate fragment test sequences."""
    if lengths is None:
        lengths = range(100, 3100, 100)
    if genes is None:
        genes = ["VP1", "3D"]
    
    # Read all sequences from the input file
    records = list(SeqIO.parse(sequences_file, "fasta"))

    # filter records in metadata file
    metadata_df = pd.read_csv(metadata_file, sep="\t")
    meta_ids = set(metadata_df["accession"] if "accession" in metadata_df.columns else metadata_df.iloc[:, 0])
    records = [r for r in records if r.id in meta_ids]

    clade_map = pd.read_csv(clades_file, sep="\t").set_index("accession")["clade"].to_dict()
    
    with open(output_file, "w") as out_handle:
        for length in lengths:
            record = random.choice(records)
            seq_len = len(record.seq)
            cl = clade_map.get(record.id, "NA")
            
            if "VP1" in genes or "3D" in genes:
                if "VP1" in genes:
                    seq1 = record.seq[2389:3315]
                    l = len(seq1) - seq1.count("-") - seq1.count("N")
                    if l > length:
                        s = random.randint(0, l - length)
                        seq1 = seq1[s:s+length]
                        header = f"{record.id}_partial_{length}_VP1_{cl}"
                        out_handle.write(f">{header}\n{seq1}\n")
                if "3D" in genes:
                    seq2 = record.seq[5926:7296]
                    l = len(seq2) - seq2.count("-") - seq2.count("N")
                    if l > length:
                        s = random.randint(0, l - length)
                        seq2 = seq2[s:s+length]
                        header = f"{record.id}_partial_{length}_3D_{cl}"
                        out_handle.write(f">{header}\n{seq2}\n")
            else:
                print(f"Gene {genes} not recognized.")
            
            while seq_len < length:
                record = random.choice(records)
                seq_len = len(record.seq)
            start = random.randint(0, seq_len - length)
            fragment_seq = record.seq[start:start+length]
            header = f"{record.id}_partial_{length}_{cl}"
            out_handle.write(f">{header}\n{fragment_seq}\n")
    
    return sum(1 for line in open(output_file) if line.startswith(">"))

def generate_recombinants(sequences_file, metadata, clades_file, evA_file, 
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
    evA_seqs = eligible_seqs(list(SeqIO.parse(evA_file, "fasta")), min_length)
    
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
    
    # Inter-typic recombinants (e.g., virus x EV-A)
    for j in range(inter_count):
        p1 = random.choice(seqs)
        p2 = random.choice(evA_seqs)
        minlen = min(len(p1.seq), len(p2.seq))
        if minlen < min_length:
            continue
        breakpoint = random.randint(1, minlen - 1)
        header = f"inter_{p1.id}_A16_{breakpoint}_{p2.id}_{p2.description.split('|')[1].strip()}"
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
    parser.add_argument("--metadata", required=True, help="Metadata TSV")
    parser.add_argument("--clades", required=True, help="Clades metadata TSV")
    parser.add_argument("--output-fragments", required=True, help="Output fragments FASTA")
    parser.add_argument("--output-recombinants", required=True, help="Output recombinants FASTA")
    parser.add_argument("--output-evA", required=False, help="Output EV-A sequences FASTA (if fetching from NCBI)")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--evA", required=False, help="EV-A sequences file (skip Entrez if provided)")
    parser.add_argument("--taxid", type=str, required=False, help="NCBI Taxon ID for EV-A (if not provided as file)")
    parser.add_argument("--virus", type=str, required=False, default="EV", help="Name of dataset virus to exclude from Entrez fetch")
    parser.add_argument("--genes", required=False, nargs='+', default=None, help="Genes for which fragments are created")
    parser.add_argument("--email", required=False, help="Email for NCBI Entrez access")
    
    args = parser.parse_args()
    random.seed(args.seed)
    
    Path(args.output_fragments).parent.mkdir(parents=True, exist_ok=True)
    Path(args.output_recombinants).parent.mkdir(parents=True, exist_ok=True)

    # Handle EV-A sequences: either use provided file or fetch from NCBI
    if args.evA:
        eva_file = args.evA
        print(f"Using provided EV-A file: {eva_file}", file=sys.stderr)
    elif args.taxid:
        eva_file = args.output_evA or Path(args.output_fragments).parent / "EV_A_fetched.fasta"
        fetched = fetch_sequences_from_entrez(
            args.taxid,
            args.virus,
            num_seqs=10,
            email=args.email,
            min_length=7000
        )
        
        if fetched:
            num_eva = write_entrez_sequences(fetched, eva_file)
            print(f"Fetched and saved {num_eva} EV-A sequences to {eva_file}", file=sys.stderr)
        else:
            print("Failed to fetch EV-A sequences, skipping recombinant generation", file=sys.stderr)
            eva_file = None
    else:
        print("Error: Either --evA file or --taxid must be provided", file=sys.stderr)
        sys.exit(1)
    
    n_frags = generate_fragments(args.sequences, args.metadata, args.clades, args.output_fragments, args.genes)
    n_recomb = generate_recombinants(args.sequences, args.metadata, args.clades, eva_file, args.output_recombinants)
    
    print(f"Generated {n_frags} fragments and {n_recomb} recombinants")