#!/usr/bin/env python3

import argparse
import csv
from Bio import SeqIO
from Bio.SeqIO import FastaIO

parser = argparse.ArgumentParser(description="Remove duplicated sequences by strain name and sequence content")
parser.add_argument("--sequences", required=True, help="Input sequences in FASTA format")
parser.add_argument("--metadata", required=True, help="Metadata TSV with strain names")
parser.add_argument("--metadata-id-columns", required=True, help="Column in metadata that matches FASTA headers")
parser.add_argument("--metadata-strain-columns", required=True, help="Column with strain names to use")
parser.add_argument("--output", required=True, help="Output FASTA with unique sequences")

args = parser.parse_args()

id = args.metadata_id_columns
strain_id = args.metadata_strain_columns

# Load metadata and map IDs to strain names
id_to_strain = {}
with open(args.metadata) as meta_fh:
    reader = csv.DictReader(meta_fh, delimiter='\t')
    for row in reader:
        if id not in row or strain_id not in row:
            continue
        strain = row[strain_id].strip()
        seq_id = row[id].strip()
        id_to_strain[seq_id] = strain

# Deduplicate sequences
seen = set()
unique_records = []

for record in SeqIO.parse(args.sequences, "fasta"):
    seq_id = record.id
    strain = id_to_strain.get(seq_id)
    if not strain:
        continue

    key = (strain, str(record.seq))
    if key not in seen:
        seen.add(key)
        unique_records.append(record)

# Write unique sequences to output
with open(args.output, "w") as out_fh:
    writer = FastaIO.FastaWriter(out_fh, wrap=None)
    writer.write_file(unique_records)