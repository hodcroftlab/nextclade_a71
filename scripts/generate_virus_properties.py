## adapted from https://github.com/neherlab/nextclade_data_workflows/blob/v3-sc2/sars-cov-2/scripts/generate_virus_properties.py
import json
from collections import defaultdict
from os import rename
import pandas as pd
from tqdm import tqdm
import os
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Annotate sequences using a genbank reference')
    parser.add_argument('--clade_meta', required=True, type=str, help='Metadata with clade and substitutions')
    parser.add_argument('--properties', required=True, type=str, help='Path to output virus properties JSON file')
    parser.add_argument('--min-prop', required=True, type=float, default=0.2, help='minimum proportion of sequences in a clade that must have a mutation for it to be considered relevant')
    parser.add_argument('--high-min-prop', required=True, type=float, default=0.7, help='minimum proportion for clades with high threshold')
    parser.add_argument('--high-prop-clades', required=True, type=str,nargs="+",  help='list of clades that should use the high proportion threshold')
    parser.add_argument('--exclude-clades', required=True, type=str, nargs="+", default="unassigned", help='list of clades to drop from mutation list')
    return parser.parse_args()

def accumulate_mutations(acc: defaultdict(int), row) -> defaultdict(int):
    try:
        for mutation in str(row).split(","):
            acc[mutation] += 1
    except:
        print(row)
        raise
    return acc


def aggregate_mutations(series) -> defaultdict(int):
    mutations = defaultdict(int)
    for row in tqdm(series):
        mutations = accumulate_mutations(mutations, row)
    return mutations

def reverse_a_dict(dict_of_lists):
    result = {}
    for k, v in dict_of_lists.items():
        for x in v:
            result.setdefault(x, []).append(k)
    return result


def main():
    args = parse_args()

    # Defines minimum proportions and counts for relevant mutations to keep
    MIN_PROPORTION = args.min_prop
    clades_with_high_proportion_threshold = args.high_prop_clades
    clades_to_drop = args.exclude_clades
    HIGH_THRESHOLD_PROPORTION = args.high_min_prop
    json_path = args.properties

    df = pd.read_csv(
        args.clade_meta,
        sep="\t",
        usecols=["clade", "substitutions"],
        parse_dates=False,
        dtype={"clade": str, "substitutions": str},
    ).dropna()

    df = df[~df["clade"].isin(clades_to_drop)]

    clade_muts = (
        df.groupby("clade")
        .substitutions.apply(aggregate_mutations)
        .dropna()
        .astype(int)
    )

    clade_muts.rename("mut_count", inplace=True)
    clade_muts

    clade_count = df.groupby("clade").count()
    clade_count.rename(columns={"substitutions": "clade_count"}, inplace=True)
    clade_count

    mutations = pd.DataFrame(clade_muts).reset_index()
    mutations.rename(columns={"level_1": "mutation"}, inplace=True)
    mutations = mutations.join(clade_count, on="clade")
    mutations

    mutations["proportion"] = mutations["mut_count"] / mutations["clade_count"]
    mutations.sort_values(
        by=["clade", "proportion"], ascending=False, inplace=True
    )
    mutations["genotype"] = mutations["mutation"].apply(lambda x: x[1:])
    mutations


    min_proportion = mutations["clade"].apply(lambda x: HIGH_THRESHOLD_PROPORTION if x in clades_with_high_proportion_threshold else MIN_PROPORTION)

    # Choose which mutations to keep
    relevant = mutations[
        (mutations["proportion"] > min_proportion)
    ]
    relevant

    #  newly relevant
    newly_relevant = relevant[mutations["proportion"] < 0.7]
    # print(newly_relevant[['clade', 'mutation', 'proportion', 'mut_count']].to_csv())

    mut_dict = {}
    for mutation, row in relevant.groupby("genotype"):
        mut_dict[mutation] = (
            row[["clade", "mut_count"]]
            .sort_values(by="mut_count", ascending=False)["clade"]
            .to_list()
        )
    mut_dict = dict(sorted(mut_dict.items(), key=lambda item: int(item[0][:-1])))
    mut_dict

    with open(json_path, "w") as f_out:

        virus_json = {
            "schemaVersion": "1.10.0",
            "nucMutLabelMap": mut_dict,
            "nucMutLabelMapReverse": dict(sorted(reverse_a_dict(mut_dict).items())),
        }

        json.dump(virus_json, f_out, indent=2)

if __name__ == "__main__":
    main()