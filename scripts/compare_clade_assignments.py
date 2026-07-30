"""
Compare clade assignments for fragment test sequences across three sources:
Nextstrain (embedded in the fragment header), Nextclade (from a nextclade run
on the fragments), and RIVM (from the RIVM typing tool's exported results).
"""
import sys
import pandas as pd

RFS_LABEL = "RFs"


def parse_fragment_name(name):
    """Extract accession, length, gene and Nextstrain clade from a fragment header.

    Headers look like "{accession}_partial_{length}_{gene}|{ns_clade}|{rivm_clade}"
    for gene fragments, or "{accession}_partial_{length}|{ns_clade}|{rivm_clade}"
    for whole-genome fragments (no gene).
    """
    prefix, ns_clade, rivm_clade = name.split("|")
    accession, rest = prefix.split("_partial_", 1)
    tokens = rest.split("_", 1)
    length = tokens[0]
    gene = tokens[1] if len(tokens) > 1 else "genome"
    return accession, length, gene, ns_clade, rivm_clade


if __name__ == "__main__":
    nextclade_tsv,rivm_file,report_tsv, recombinant_clades_arg, summary_file = sys.argv[1:6]
    #%%
    # parse recombinant clade list
    recombinant_clades = recombinant_clades_arg.split(",") if recombinant_clades_arg else []
    # get nextclade tsv
    nc_df = pd.read_csv(nextclade_tsv, sep="\t", low_memory=False)
    # get RIVM df
    rivm_df = pd.read_csv(rivm_file)

    # get virus
    virus = (rivm_df.type.value_counts().index)[0]
    
    ## get report
    rivm_clade = dict(zip(rivm_df["name"], rivm_df["VP1 subgenogroup"].fillna("NA")))

    rows = []
    for _, row in nc_df.iterrows():
        name = row["seqName"]
        if "_partial_" not in name:
            continue
        accession, length, gene, ns_clade, rivm_og = parse_fragment_name(name)
        rows.append({
            "fragment": name,
            "accession": accession,
            "length": int(length),
            "gene": gene,
            "nextstrain_clade": ns_clade,
            "og_rivm": rivm_og,
            "nextclade_clade": row.get("clade", "NA") if pd.notna(row.get("clade")) else "NA",
            "rivm_clade": rivm_clade.get(name, "NA"),
        })

    report_df = pd.DataFrame(rows)
    report_df.to_csv(report_tsv, sep="\t", index=False)
    print(f"Wrote {len(report_df)} rows to {report_tsv}")

    # get summary file
    nc_df["genome_position"] = (nc_df["alignmentStart"] + nc_df["alignmentEnd"]) / 2
    position = dict(zip(nc_df["seqName"], nc_df["genome_position"]))

    recombinant_labels = set(recombinant_clades) | {RFS_LABEL}

    report_df["genome_position"] = report_df["fragment"].map(position)
    report_df["recombinant"] = report_df["nextclade_clade"].isin(recombinant_labels).map({True: 1, False: 0})

    report_df.replace("Could not assign","NA", inplace=True)

    report_df["nextclade_success"] = report_df["nextclade_clade"].map({"NA": 0}).fillna(1)
    report_df["RIVM_success"] = report_df["rivm_clade"].map({"NA": 0}).fillna(1)

    # Exclude fragments with no ground truth from "correct": otherwise a row where
    # neither the tool nor the reference has an answer (NA == NA) is spuriously
    # counted as a correct call.
    report_df["nextclade_correct"] = (
        (report_df["nextclade_clade"] == report_df["nextstrain_clade"]) & (report_df["nextstrain_clade"] != "NA")
    ).map({True: 1, False: 0})
    report_df["RIVM_correct"] = (
        (report_df["rivm_clade"] == report_df["og_rivm"]) & (report_df["og_rivm"] != "NA")
    ).map({True: 1, False: 0})

    report_df.sort_values("length", inplace=True)

    #%%
    summary = pd.DataFrame([
        {
            "virus": virus,
            "method": "nextclade",
            "n_evaluable": (report_df["nextstrain_clade"] != "NA").sum(),
            "n_assigned": report_df["nextclade_success"].sum(),
            "n_correct": report_df["nextclade_correct"].sum(),
        },
        {
            "virus": virus,
            "method": "RIVM",
            "n_evaluable": (report_df["nextstrain_clade"] != "NA").sum(),
            "n_assigned": report_df["RIVM_success"].sum(),
            "n_correct": report_df["RIVM_correct"].sum(),
        },
    ])
    summary["pct_assigned"] = 100 * summary["n_assigned"] / summary["n_evaluable"]
    summary["pct_correct"] = 100 * summary["n_correct"] / summary["n_evaluable"]
    summary["pct_correct_given_assigned"] = 100 * summary["n_correct"] / summary["n_assigned"]

    summary.to_csv(summary_file, sep="\t", index=False)
    print(f"Wrote {len(summary)} rows to {summary_file}")

# %%
