#!/usr/bin/env python3

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def main():
    parser = argparse.ArgumentParser(description='Fill gaps in root sequence with reference')
    parser.add_argument('--ancestral', required=True, help='Ancestral sequences file')
    parser.add_argument('--reference', required=True, help='Reference sequence file')
    parser.add_argument('--output', required=True, help='Output root sequence file')
    
    args = parser.parse_args()
    
    # Find root sequence
    root_seq = None
    with open(args.ancestral, 'r') as f:
        for record in SeqIO.parse(f, 'fasta'):
            if 'NODE_0000000' in record.id:
                root_seq = record
                break
    
    if not root_seq:
        raise ValueError('NODE_0000000 not found in ancestral sequences')
    
    # Read reference
    with open(args.reference, 'r') as f:
        ref_seq = next(SeqIO.parse(f, 'fasta'))
    
    print(f'Root: {root_seq.id} ({len(root_seq.seq)} bp)')
    print(f'Reference: {ref_seq.id} ({len(ref_seq.seq)} bp)')
    
    # Fill gaps position by position
    root_str = str(root_seq.seq)
    ref_str = str(ref_seq.seq)
    
    min_len = min(len(root_str), len(ref_str))
    if len(root_str) != len(ref_str):
        print(f'WARNING: Length mismatch - Root: {len(root_str)}, Ref: {len(ref_str)}')
    
    filled_seq = []
    gaps_filled = 0
    ns_filled = 0
    
    for i in range(min_len):
        root_nt = root_str[i]
        ref_nt = ref_str[i]
        
        if root_nt == '-':
            filled_seq.append(ref_nt)
            gaps_filled += 1
        elif root_nt.upper() == 'N':
            filled_seq.append(ref_nt)
            ns_filled += 1
        else:
            filled_seq.append(root_nt)
    
    # Handle remaining sequence
    if len(root_str) > min_len:
        filled_seq.extend(root_str[min_len:])
    
    filled_sequence = ''.join(filled_seq)
    
    # Create and write output
    filled_record = SeqRecord(
        Seq(filled_sequence),
        id="ancestral_sequence",
        description=""
    )
    
    with open(args.output, 'w') as f:
        SeqIO.write(filled_record, f, 'fasta')
    
    print(f'Successfully filled {gaps_filled} gaps and {ns_filled} Ns')
    print(f'Final root sequence: {len(filled_sequence)} bp')

if __name__ == '__main__':
    main()