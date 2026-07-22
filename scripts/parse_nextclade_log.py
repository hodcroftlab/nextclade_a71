"""
Parse Nextclade log to extract failed sequences and seed alignment coverage stats.
"""
import re
from pathlib import Path
from collections import Counter
import matplotlib.pyplot as plt
from Bio import SeqIO
import json
import numpy as np

def parse_nextclade_log(log_file):
    """
    Extract failed sequence info and seed alignment percentages from Nextclade log.
    
    Returns:
        failed_sequences (list): List of sequence IDs that failed alignment
        coverage_pcts (list): List of seed alignment coverage percentages
    """
    failed_sequences = []
    coverage_pcts = []
    
    # Pattern: "seed alignment covers XX.XX% of the query sequence"
    coverage_pattern = r"seed alignment covers ([\d.]+)% of the query sequence"
    
    with open(log_file) as f:
        for line in f:
            # Extract sequence ID from warning lines
            if "[W]" in line and "Unable to align" in line:
                # Pattern: "In sequence #XXXX 'SEQID'" — extract just the ID before any space/pipe
                seq_match = re.search(r"In sequence #\d+ '([^']+)'", line)
                if seq_match:
                    full_id = seq_match.group(1)
                    # Take only the accession part before pipe/space
                    clean_id = full_id.split()[0].split('|')[0]
                    failed_sequences.append(clean_id)
            
            # Extract coverage percentage
            cov_match = re.search(coverage_pattern, line)
            if cov_match:
                coverage_pcts.append(float(cov_match.group(1)))
    
    return failed_sequences, coverage_pcts

def write_failed_sequences_fasta(failed_sequences, fasta_file, output_dir):
    """
    Write all failed sequences to a new FASTA file.
    """
    output_dir = Path(output_dir)
    output_fasta = output_dir / "failed_sequences.fasta"
    
    # Read original FASTA and extract failed sequences
    failed_records = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id in failed_sequences and "|" not in record.description and "int" not in record.id: # remove non-EV sequences (otherwise alignment will not make sense)
            failed_records.append(record)   
    
    # Write to output file
    SeqIO.write(failed_records, output_fasta, "fasta")
    print(f"Failed sequences written to: {output_fasta}\n")

def categorize_test_sequences(failed_sequences, fasta_file, qc_status):
    """
    Categorize sequences into: EV, non-EV-A, fragments, inter-recombinants, intra-recombinants.
    Returns failure counts and QC stats for each category.
    """
    categories = {
        short_name: [],
        'non_EV_A': [],
        'EV_A': [],
        'fragments': [],
        'inter_recombinants': [],
        'intra_recombinants': []
    }
    
    # Read all sequences and categorize them
    all_seqs = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        all_seqs[record.id] = record.description
    
    for seq_id in all_seqs:
        description = all_seqs[seq_id]

        # Determine category based on sequence ID (check in order of specificity)
        if seq_id.startswith('inter_'):  # Inter-recombinants first
            categories['inter_recombinants'].append(seq_id)
        elif seq_id.startswith('intra_'):  # Intra-recombinants
            categories['intra_recombinants'].append(seq_id)
        elif '_partial_' in seq_id:  # Fragments
            categories['fragments'].append(seq_id)
        elif 'EV-A' in description or 'CVA' in description:  # EV-A sequences (but not EV)
            if virus_name not in description and short_name not in description:
                categories['EV_A'].append(seq_id)
        elif '|' in description:  # Non-EV (has pipe symbol)
            categories['non_EV_A'].append(seq_id)
        else:  # EV
            categories[short_name].append(seq_id)
    
    # Count failures per category
    results = {}
    for cat_name, seq_list in categories.items():
        failed_in_cat = [s for s in seq_list if s in failed_sequences]
        
        # Get QC status by matching on description (since qc_status uses full description as key)
        qc_in_cat = []
        for seq_id in seq_list:
            desc = all_seqs[seq_id]
            # Try matching by id first, then by description
            if seq_id in qc_status:
                qc_in_cat.append(qc_status[seq_id])
            elif desc in qc_status:
                qc_in_cat.append(qc_status[desc])
        
        results[cat_name] = {
            'total': len(seq_list),
            'failed': len(failed_in_cat),
            'passed': len(seq_list) - len(failed_in_cat),
            'qc_stats': Counter(qc_in_cat),
            'seq_ids': seq_list
        }
    
    return results

def print_test_summary(test_results, output_dir):
    """
    Print and save test sequence summary.
    """
    output_dir = Path(output_dir)
    
    print(f"\n{'='*70}")
    print(f"TEST SEQUENCE SUMMARY")
    print(f"{'='*70}")
    print(f"{'Category':<20} {'Total':>8} {'Passed':>8} {'Failed':>8} {'Pass %':>10}")
    print(f"-" * 70)
    
    table_path = output_dir / "test_sequences_summary.txt"
    with open(table_path, 'w') as f:
        f.write(f"{'Category':<20} {'Total':>8} {'Passed':>8} {'Failed':>8} {'Pass %':>10}\n")
        f.write(f"-" * 70 + "\n")
        
        for cat_name, stats in test_results.items():
            total = stats['total']
            failed = stats['failed']
            passed = stats['passed']
            pass_pct = 100 * passed / total if total > 0 else 0
            
            print(f"{cat_name:<20} {total:>8} {passed:>8} {failed:>8} {pass_pct:>9.1f}%")
            f.write(f"{cat_name:<20} {total:>8} {passed:>8} {failed:>8} {pass_pct:>9.1f}%\n")
    
    print(f"{'='*70}\n")
    print(f"Test sequence summary saved to: {table_path}\n")
    
    return table_path

def plot_test_qc_distribution(test_results, output_dir):
    """
    Create bar plot of QC status by test sequence category.
    """
    output_dir = Path(output_dir)
    
    statuses = ['good', 'mediocre', 'bad', 'failed']
    categories = list(test_results.keys())
    
    # Prepare data for grouped bar chart
    data_by_status = {status: [] for status in statuses}
    for cat in categories:
        qc_counter = test_results[cat]['qc_stats']
        for status in statuses:
            data_by_status[status].append(qc_counter.get(status, 0))
    
    # Create grouped bar chart
    fig, ax = plt.subplots(figsize=(12, 6))
    
    x = range(len(categories))
    width = 0.2
    colors = ['green', 'orange', 'red', 'gray']
       
    for xi, cat in enumerate(categories):
        qc_counter = test_results[cat]['qc_stats']
        
        # Keep only statuses that actually have values
        present_statuses = [s for s in statuses if qc_counter.get(s, 0) > 0]
        n_present = len(present_statuses)
        
        for j, status in enumerate(present_statuses):
            value = qc_counter.get(status, 0)
            
            # Recompute centered offsets per category
            offset = (j - (n_present - 1) / 2) * width
            
            ax.bar(
                xi + offset,
                value,
                width,
                label=status if xi == 0 else "",  # avoid duplicate legend
                color=colors[statuses.index(status)],
                alpha=0.7,
                edgecolor='black'
            )   

    ax.set_xlabel('Sequence Category', fontsize=12)
    ax.set_ylabel('Log Scale of Number of Sequences', fontsize=12)
    ax.set_title('QC Status Distribution by Test Sequence Category', fontsize=13)
    ax.set_xticks(x)
    ax.set_xticklabels(categories, rotation=15, ha='right')
    handles, labels = ax.get_legend_handles_labels()
    unique = dict(zip(labels, handles))
    ax.legend(unique.values(), unique.keys())
    ax.grid(axis='y', alpha=0.3)
    ax.set_yscale('log')
    
    plt.tight_layout()
    plot_path = output_dir / "test_sequences_qc_distribution.png"
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    print(f"Test QC distribution plot saved to: {plot_path}\n")

def summarize_results(failed_sequences, coverage_pcts, total_sequences, seq_lengths, qc_status, fasta_file, output_dir):
    """
    Print summary stats and generate histograms.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    passed = total_sequences - len(failed_sequences)
    fail_pct = 100 * len(failed_sequences) / total_sequences if total_sequences > 0 else 0
    pass_pct = 100 * passed / total_sequences if total_sequences > 0 else 0
    
    # Get lengths of failed sequences
    failed_lengths = [seq_lengths.get(seq_id, None) for seq_id in failed_sequences]
    failed_lengths = [l for l in failed_lengths if l is not None]
    
    # QC stats for ALL sequences
    all_qc_status = list(qc_status.values())
    all_qc_counter = Counter(all_qc_status)
    
    # QC stats for only failed sequences (if they have QC data)
    failed_with_qc = [qc_status.get(seq_id) for seq_id in failed_sequences if seq_id in qc_status]
    failed_qc_counter = Counter(failed_with_qc) if failed_with_qc else Counter()
    
    # Summary stats
    print(f"\n{'='*60}")
    print(f"NEXTCLADE LOG SUMMARY")
    print(f"{'='*60}")
    print(f"Total sequences: {total_sequences}")
    print(f"Passed: {passed} ({pass_pct:.1f}%)")
    print(f"Failed: {len(failed_sequences)} ({fail_pct:.1f}%)")
    
    if coverage_pcts:
        print(f"\nSeed alignment coverage statistics (failed sequences only):")
        print(f"  Min:  {min(coverage_pcts):.2f}%")
        print(f"  Max:  {max(coverage_pcts):.2f}%")
        print(f"  Mean: {sum(coverage_pcts)/len(coverage_pcts):.2f}%")
        sorted_cov = sorted(coverage_pcts)
        median_idx = len(sorted_cov) // 2
        print(f"  Median: {sorted_cov[median_idx]:.2f}%")
    else:
        print(f"\nNo failed sequences to analyze.")
    
    if failed_lengths:
        print(f"\nSequence length statistics (failed sequences only):")
        print(f"  Min:  {min(failed_lengths)} bp")
        print(f"  Max:  {max(failed_lengths)} bp")
        print(f"  Mean: {sum(failed_lengths)/len(failed_lengths):.0f} bp")
        sorted_len = sorted(failed_lengths)
        median_len_idx = len(sorted_len) // 2
        print(f"  Median: {sorted_len[median_len_idx]} bp")
    
    print(f"\nQC Overall Status (all sequences):")
    for status in ['good', 'mediocre', 'bad', 'failed']:
        count = all_qc_counter.get(status, 0)
        pct = 100 * count / len(all_qc_status) if all_qc_status else 0
        print(f"  {status:10}: {count:5} ({pct:5.1f}%)")
    
    print(f"{'='*60}\n")
    
    # Test sequence analysis
    test_results = categorize_test_sequences(failed_sequences, fasta_file, qc_status)
    print_test_summary(test_results, output_dir)
    plot_test_qc_distribution(test_results, output_dir)
    
    # Write failed sequences to FASTA
    write_failed_sequences_fasta(failed_sequences, fasta_file, output_dir)
    
    # Histograms
    if coverage_pcts:
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        
        # Coverage histogram
        axes[0, 0].hist(coverage_pcts, bins=range(int(min(coverage_pcts)), int(max(coverage_pcts)) + 2), edgecolor='black', alpha=0.7, color='steelblue')
        axes[0, 0].set_xlabel('Seed Alignment Coverage (%)', fontsize=11)
        axes[0, 0].set_ylabel('Frequency', fontsize=11)
        axes[0, 0].set_title('Seed Alignment Coverage Distribution', fontsize=12)
        axes[0, 0].grid(axis='y', alpha=0.3)
        
        # Length histogram
        if failed_lengths:
            axes[0, 1].hist(failed_lengths, bins=100, edgecolor='black', alpha=0.7, color='coral')
            axes[0, 1].set_xlabel('Sequence Length (bp)', fontsize=11)
            axes[0, 1].set_ylabel('Frequency', fontsize=11)
            axes[0, 1].set_title('Failed Sequence Length Distribution', fontsize=12)
            axes[0, 1].grid(axis='y', alpha=0.3)
        
        # QC status bar chart (all sequences)
        statuses = ['good', 'mediocre', 'bad', 'failed']
        counts = [all_qc_counter.get(s, 0) for s in statuses]
        colors = ['green', 'orange', 'red', 'gray']
        axes[1, 0].bar(statuses, counts, color=colors, alpha=0.7, edgecolor='black')
        axes[1, 0].set_ylabel('Frequency', fontsize=11)
        axes[1, 0].set_title('QC Overall Status (all sequences)', fontsize=12)
        axes[1, 0].grid(axis='y', alpha=0.3)
        
        # Length by QC status (box plot for all sequences)
        qc_by_length = {status: [] for status in statuses}
        for seq_id, length in seq_lengths.items():
            status = qc_status.get(seq_id, 'unknown')
            if status in statuses:  # Only include known statuses
                qc_by_length[status].append(length)
        
        data_to_plot = [qc_by_length[s] for s in statuses if qc_by_length[s]]
        labels = [s for s in statuses if qc_by_length[s]]
        axes[1, 1].boxplot(data_to_plot, tick_labels=labels)
        axes[1, 1].set_ylabel('Sequence Length (bp)', fontsize=11)
        axes[1, 1].set_title('Length Distribution by QC Status (all sequences)', fontsize=12)
        axes[1, 1].grid(axis='y', alpha=0.3)
        
        plt.tight_layout()
        histogram_path = output_dir / "nextclade_failed_sequences.png"
        plt.savefig(histogram_path, dpi=300, bbox_inches='tight')
        print(f"Histograms saved to: {histogram_path}\n")
        
        # Coverage table
        coverage_counter = Counter([f"{int(pct)}-{int(pct)+1}" for pct in coverage_pcts])
        table_path = output_dir / "coverage_table.txt"
        with open(table_path, 'w') as f:
            f.write("Coverage Range (%)  | Count\n")
            f.write("-" * 30 + "\n")
            for bin_label in sorted(coverage_counter.keys()):
                f.write(f"{bin_label:18} | {coverage_counter[bin_label]}\n")
        print(f"Coverage table saved to: {table_path}\n")
        
        # QC status table (all sequences)
        qc_table_path = output_dir / "qc_status_table.txt"
        with open(qc_table_path, 'w') as f:
            f.write("QC Status   | Count | Percentage\n")
            f.write("-" * 40 + "\n")
            for status in statuses:
                count = all_qc_counter.get(status, 0)
                pct = 100 * count / len(all_qc_status) if all_qc_status else 0
                f.write(f"{status:11} | {count:5} | {pct:6.1f}%\n")
        print(f"QC status table saved to: {qc_table_path}\n")

def extract_mutation_stats_from_tree(tree_file):
    """
    Extract private mutation statistics from phylogenetic tree JSON.
    
    Counts mutations on each branch (private mutations relative to parent node).
    Mutations are stored in branch_attrs.mutations["nuc"] as a list of mutation strings.
    
    Returns:
        typical (float): Median of mutations per branch
        cutoff (float): 97th percentile of mutations per branch
    """
    
    with open(tree_file) as f:
        tree_data = json.load(f)
    
    # Extract mutation counts from all nodes
    mutation_counts = []
    
    def traverse_tree(node):
        # Count mutations on branch leading to this node
        if 'branch_attrs' in node and 'mutations' in node['branch_attrs']:
            branch_muts = node['branch_attrs']['mutations']
            # Count only nucleotide mutations ("nuc" key)
            nuc_muts = branch_muts.get('nuc', [])
            count = len(nuc_muts) if isinstance(nuc_muts, list) else 0
            mutation_counts.append(count)
        
        # Recursively traverse children
        if 'children' in node:
            for child in node['children']:
                traverse_tree(child)
    
    # Start from tree root
    if 'tree' in tree_data:
        traverse_tree(tree_data['tree'])
    
    if not mutation_counts:
        print("Warning: No mutations found in tree")
        return None, None
    
    typical = float(np.percentile(mutation_counts, 50))  # median
    cutoff = float(np.percentile(mutation_counts, 97))   # 97th percentile
    
    return typical, cutoff

if __name__ == "__main__":
    import sys
    import pandas as pd

    log_file = sys.argv[1] if len(sys.argv) > 1 else "test_out/test.log"
    fasta_file = sys.argv[2] if len(sys.argv) > 2 else "data/sequences.fasta"
    tsv_file = sys.argv[3] if len(sys.argv) > 3 else "test_out/nextclade.tsv"
    output_dir = sys.argv[4] if len(sys.argv) > 4 else "test_out"
    virus_name = sys.argv[5] if len(sys.argv) > 5 else "Enterovirus"
    tree_file = sys.argv[6] if len(sys.argv) > 6 else "out-dataset/tree.json"
    short_name = sys.argv[7] if len(sys.argv) > 7 else "EV"

    # Extract mutation statistics from tree
    typical, cutoff = extract_mutation_stats_from_tree(tree_file)
    if typical is not None and cutoff is not None:
        print(f"\nPrivate Mutations Statistics from Tree:")
        print(f"  Typical (median): {typical:.1f}")
        print(f"  Cutoff (97th percentile): {cutoff:.1f}\n")
    
    seq_lengths = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_lengths[record.id] = len(record.seq)
    
    total_seqs = len(seq_lengths)
    failed_seqs, coverage_vals = parse_nextclade_log(log_file)

    df = pd.read_csv(tsv_file, sep='\t', low_memory=False)
    qc_status = dict(zip(df['seqName'], df['qc.overallStatus'].fillna('failed')))

    summarize_results(failed_seqs, coverage_vals, total_seqs, seq_lengths, qc_status, fasta_file, output_dir)