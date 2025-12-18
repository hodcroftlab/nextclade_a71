# Set the parameters
REFERENCE_ACCESSION =   "U22521"
TAXON_ID =              39054
GENES =                 ["VP4", "VP2", "VP3", "VP1", "2A", "2B", "2C", "3A", "3B", "3C", "3D"]
ALLOWED_DIVERGENCE =    "3000" # was 
MIN_DATE =              "1960-01-01"
MIN_LENGTH =            "6000" # was 6000 for whole genome build on Nextstrain
MAX_SEQS =              "1000" # tree will be subsampled
ROOTING =               "ancestral_sequence"  # mid_point, outgroup, reference, ancestral sequence
ID_FIELD=               "accession" # either accession or strain, used for meta-id-column in augur

# Set the paths
SEQUENCES =             "data/sequences.fasta"
METADATA =              "data/metadata.tsv"

GFF_PATH =              "dataset/genome_annotation.gff3" 
PATHOGEN_JSON =         "dataset/pathogen.json"
README_PATH =           "dataset/README.md"
CHANGELOG_PATH =        "dataset/CHANGELOG.md"
REFERENCE_PATH =        "dataset/reference.fasta"

GENBANK_PATH =          "resources/reference.gbk"
AUSPICE_CONFIG =        "resources/auspice_config.json"
EXCLUDE =               "resources/exclude.txt"
CLADES =                "resources/clades.tsv"
ACCESSION_STRAIN =      "resources/accession_strain.tsv"
EXTRA_META =            "resources/meta_public.tsv"
INCLUDE_EXAMPLES =      "resources/include_examples.txt"
COLORS =                "resources/colors.tsv"
COLORS_SCHEMES =        "resources/color_schemes.tsv"
INFERRED_ANCESTOR =     "resources/inferred-root.fasta"

STAR_ROOT = True                   # whether to use star-like rooting method
FETCH_SEQUENCES = True              # whether to fetch sequences from NCBI Virus via ingest workflow
STATIC_ANCESTRAL_INFERRENCE = True  # whether to use the static inferred ancestral sequence
INFERRENCE_RERUN = False            # whether to rerun the inference of the ancestral sequence worfkflow (inferred-root)

INFERRED_SEQ_PATH = "results/sequences_with_ancestral.fasta" if STATIC_ANCESTRAL_INFERRENCE else SEQUENCES
INFERRED_META_PATH = "results/metadata_with_ancestral.tsv" if STATIC_ANCESTRAL_INFERRENCE else "results/metadata.tsv"
TREE = "results/tree.nwk" if not STAR_ROOT else "results/star_tree.nwk"

include: "scripts/workflow_messages.snkm"
configfile: PATHOGEN_JSON

rule all:
    input:
        auspice = "results/auspice.json",
        augur_jsons = "test_out/",
        data = "dataset.zip",
        seqs = "results/example_sequences.fasta",
        json = "out-dataset/pathogen.json",
        **({"root": INFERRED_ANCESTOR} if STATIC_ANCESTRAL_INFERRENCE else {})


rule testing:
    input:
        "testing/EV-A71_fragments.fasta",
        "testing/EV-A71_recombinants.fasta"


if FETCH_SEQUENCES == True:
    rule fetch:
        input:
            dir = "ingest"
        output:
            sequences=SEQUENCES,
            metadata=METADATA
        threads: workflow.cores
        shell:
            """
            cd {input.dir} 
            snakemake --cores {threads} all
            cd ../
            """

rule curate:
    message:
        """
        Cleaning up metadata with augur merge & augur curate
        """
    input:
        meta = METADATA,  # Path to input metadata file
        strains = ACCESSION_STRAIN,  # Strain - accession lookup table
        public = EXTRA_META,
    params:
        strain_id_field = ID_FIELD,
        date_format = ['%Y-%m-%d','%Y','XX-%m-%Y','%Y-%m-%dT%H:%M:%SZ','201X-XX-XX', 
            'XX-XX-%Y', 'XX-XX-XXXX','%m.%Y','%m-%Y', '%d.%m.%Y', "%b-%Y", "%d-%b-%Y"],  # Date format for metadata
        date_fields = ["date","date_released","date_updated"]
    output:
        metadata = "results/metadata.tsv",  # Final output file for publications metadata
    shell:
        """
        augur merge --metadata metadata={input.meta} strains={input.strains} public={input.public}\
            --metadata-id-columns {params.strain_id_field} \
            --output-metadata metadata.tmp
        
        augur curate normalize-strings \
            --metadata metadata.tmp \
            --id-column {params.strain_id_field} \
        | augur curate format-dates \
            --id-column {params.strain_id_field} \
            --no-mask-failure \
            --expected-date-formats {params.date_format} \
            --date-fields {params.date_fields} \
            --output-metadata {output.metadata}

        rm metadata.tmp
        """

rule add_reference_to_include:
    """
    Create an include file for augur filter
    """
    input:
        "resources/include.txt",
    output:
        "results/include.txt",
    shell:
        """
        cat {input} >> {output}
        echo "{REFERENCE_ACCESSION}" >> {output}
        echo ancestral_sequence >> {output}
        """

if STATIC_ANCESTRAL_INFERRENCE and INFERRENCE_RERUN:
    rule static_inferrence:
        message:
            """
            Running "inferred-root" snakefile for inference of the ancestral root. 
            This reference will be included in the Nextclade reference tree.
            WARNING: This will overwrite your {output.seq} & {output.meta} files!
            """
        input:
            dir = "inferred-root",
            dataset_path = "dataset",
            meta = rules.curate.output.metadata,
            seq = SEQUENCES,
            meta_ancestral = "resources/static_inferred_metadata.tsv",
            include = "results/include.txt"
        params:
            strain_id_field = ID_FIELD,
        output:
            inref = INFERRED_ANCESTOR,
            seq = INFERRED_SEQ_PATH,
            meta = INFERRED_META_PATH,
        threads: workflow.cores
        shell:
            r"""
            set -euo pipefail

            echo "Cleaning previous results..."
            rm -rf {input.dir}/results/* {input.dir}/resources/inferred-root.fasta

            echo "Running inferred-root workflow..."
            cd {input.dir}
            snakemake --cores {threads} all_sub
            cd - > /dev/null

            echo "Combining sequences with ancestral root..."
            cat {input.seq} {output.inref} > {output.seq}

            echo "Merging metadata..."
            augur merge \
                --metadata metadata={input.meta} ancestral={input.meta_ancestral} \
                --metadata-id-columns {params.strain_id_field} \
                --output-metadata {output.meta}

            echo "Static ancestral inference completed successfully!"
            """

if STATIC_ANCESTRAL_INFERRENCE and not INFERRENCE_RERUN:
    rule add_ancestral:
        input:
            meta = rules.curate.output.metadata,
            seq = SEQUENCES,
            meta_ancestral = "resources/static_inferred_metadata.tsv",
            inref = INFERRED_ANCESTOR,
            RIVM = "resources/subgenotypes_rivm.csv",
        output:
            seq = INFERRED_SEQ_PATH,
            meta = INFERRED_META_PATH,
        params:
            strain_id_field="accession",
        shell:
            """
            echo "Combining sequences with ancestral root..."
            cat {input.seq} {input.inref} > {output.seq}

            echo "Merging metadata..."
            augur merge \
                --metadata metadata={input.meta} ancestral={input.meta_ancestral} rivm={input.RIVM} \
                --metadata-id-columns {params.strain_id_field} \
                --output-metadata {output.meta}

            echo "Static ancestral sequence imported successfully!"
            """

rule index_sequences:
    message:
        """
        Creating an index of sequence composition for filtering
        """
    input:
        sequences = INFERRED_SEQ_PATH,
    output:
        sequence_index = "results/sequence_index.tsv"
    shell:
        """
        augur index \
            --sequences {input.sequences} \
            --output {output.sequence_index}
        """

rule filter:
    """
    Exclude sequences from before {MIN_DATE} and subsample to {MAX_SEQS} sequences.
    Only take sequences longer than {MIN_LENGTH}
    """
    input:
        sequences = INFERRED_SEQ_PATH,
        sequence_index = rules.index_sequences.output.sequence_index,
        metadata = INFERRED_META_PATH,
        include = rules.add_reference_to_include.output,
    output:
        filtered_sequences = "results/filtered_sequences_raw.fasta",
        filtered_metadata = "results/filtered_metadata_raw.tsv",
    params: 
        min_date="" if MIN_DATE == "" else "--min-date " + MIN_DATE,
        min_length="" if MIN_LENGTH == "" else "--min-length " + MIN_LENGTH,
        max_seqs=MAX_SEQS,
        categories = "country year", #TODO: add subsampling per category?
        strain_id_field = ID_FIELD,
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --sequence-index {input.sequence_index} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id_field} \
            {params.min_length} \
            {params.min_date} \
            --include {input.include} \
            --group-by {params.categories} \
            --subsample-max-sequences {params.max_seqs} \
            --output-sequences {output.filtered_sequences} \
            --output-metadata {output.filtered_metadata}
        """
        
rule align:
    message:
        """
        Aligning sequences to {input.reference} using Nextclade3.
        """
    input:
        sequences = rules.filter.output.filtered_sequences,
        reference = REFERENCE_PATH,
        annotation = GFF_PATH,
    output:
        alignment = "results/aligned.fasta",
        tsv = "results/nextclade.tsv",
    params:
        translation_template = lambda w: "results/translations/cds_{cds}.translation.fasta",
        penalty_gap_extend = config["alignmentParams"]["penaltyGapExtend"],
        penalty_gap_open = config["alignmentParams"]["penaltyGapOpen"],
        penalty_gap_open_in_frame = config["alignmentParams"]["penaltyGapOpenInFrame"],
        penalty_gap_open_out_of_frame = config["alignmentParams"]["penaltyGapOpenOutOfFrame"],
        kmer_length = config["alignmentParams"]["kmerLength"],
        kmer_distance = config["alignmentParams"]["kmerDistance"],
        min_match_length = config["alignmentParams"]["minMatchLength"],
        allowed_mismatches = config["alignmentParams"]["allowedMismatches"],
        min_length = config["alignmentParams"]["minLength"],
        gap_alignment_side = config["alignmentParams"]["gapAlignmentSide"],  
        min_seed_cover = config["alignmentParams"]["minSeedCover"],
    shell:
        """
        nextclade3 run \
        -j {threads} \
        {input.sequences} \
        --input-ref {input.reference} \
        --input-annotation {input.annotation} \
        --alignment-preset high-diversity \
        --penalty-gap-open {params.penalty_gap_open} \
        --penalty-gap-extend {params.penalty_gap_extend} \
        --penalty-gap-open-in-frame {params.penalty_gap_open_in_frame} \
        --penalty-gap-open-out-of-frame {params.penalty_gap_open_out_of_frame} \
        --gap-alignment-side {params.gap_alignment_side} \
        --kmer-length {params.kmer_length} \
        --kmer-distance {params.kmer_distance} \
        --min-match-length {params.min_match_length} \
        --allowed-mismatches {params.allowed_mismatches} \
        --min-seed-cover {params.min_seed_cover} \
        --min-length {params.min_length} \
        --max-alignment-attempts 5 \
        --include-reference false \
        --output-tsv {output.tsv} \
        --output-translations {params.translation_template} \
        --output-fasta {output.alignment} 
        """

rule get_outliers:
    """
    Automatically identify sequences with >{ALLOWED_DIVERGENCE} substitutions
    (likely to be sequencing errors or low quality/misannotated sequences) and put them in outliers.txt
    """
    input:
        nextclade = rules.align.output.tsv,
    output:
        outliers = "results/outliers.txt",
        tmp = "tmp/outliers.txt",
    params:
        allowed_divergence = lambda w: ALLOWED_DIVERGENCE,
    shell:
        """
        tsv-filter -H -v --is-numeric totalSubstitutions {input.nextclade} \
        > {output.tmp}
        tsv-filter -H \
            --is-numeric totalSubstitutions \
            --gt totalSubstitutions:{params.allowed_divergence} \
            {input.nextclade} \
        | tail -n +2 >> {output.tmp}
        cat {output.tmp} \
        | tsv-select -H -f seqName \
        | tail -n +2 > {output.outliers}
        """


rule exclude:
    """
    Rule to allow for manual and automatic exclusion of sequences
    without triggering a new subsampling that could
    surface new bad sequences resulting in an infinite loop
    """
    input:
        sequences = rules.align.output.alignment,
        sequence_index = rules.index_sequences.output.sequence_index,
        metadata = INFERRED_META_PATH,
        exclude = EXCLUDE,
        outliers = rules.get_outliers.output.outliers,
        example = INCLUDE_EXAMPLES,
    params:
        strain_id_field = ID_FIELD,
    output:
        filtered_sequences = "results/filtered_aligned.fasta",
        filtered_metadata = "results/filtered_metadata.tsv",
        strains = "results/tree_strains.txt",
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --sequence-index {input.sequence_index} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id_field} \
            --exclude {input.exclude} {input.outliers} {input.example} \
            --output-sequences {output.filtered_sequences} \
            --output-metadata {output.filtered_metadata} \
            --output-strains {output.strains}
        """

rule tree:
    message:
        """
        Creating a maximum likelihood tree
        """
    input:
        alignment = rules.exclude.output.filtered_sequences,
    output:
        tree = "results/tree_raw.nwk",
    threads: workflow.cores
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --nthreads {threads}\
            --output {output.tree} \
        """

rule refine:
    input:
        tree=rules.tree.output.tree,
        alignment=rules.exclude.output.filtered_sequences,
    output:
        tree="results/tree.nwk",
        node_data="results/branch_lengths.json",
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --root {ROOTING} \
            --keep-polytomies \
            --divergence-unit mutations-per-site \
            --output-node-data {output.node_data} \
            --output-tree {output.tree}
        """

if STAR_ROOT==True:
    """
    Rule to root the tree using a star-like rooting method.
    """
    rule star_like_rooting:
        input:
            tree=ancient(rules.refine.output.tree),
            clades="resources/clade_map.tsv", # TODO: the "clade_map" has to be downloaded from Nextstrain auspice metadata: 
                                                # "accession\tclade" format - use a tree with all sequences if possible `df = df.loc[:,["accession", "clade_membership"]]`
            recombinant_accessions="resources/recombinants.tsv",  # Format: "accession" (one per line, with header)
            alignment=rules.exclude.output.filtered_sequences,
        output:
            tree="results/star_tree.nwk",
            node_data="results/branch_lengths.json",
        params:
            strain_id_field=ID_FIELD,
            recombinant_clades = ["C2r", "C1-like","C2-like", "E", "F", "A", "C6"],
            root_name="NODE_0000000"
        log:
            "logs/star_like_rooting.log"
        shell:
            """
            python scripts/star_like_rooting.py \
                --input_tree {input.tree} \
                --input_clades {input.clades} \
                --recombinant_accessions {input.recombinant_accessions} \
                --output_tree {output.tree} \
                --strain_id_field {params.strain_id_field} \
                --recombinant_clades {params.recombinant_clades} \
                --root_name {params.root_name} \
                2>&1 | tee {log}

            augur refine \
            --tree {output.tree} \
            --alignment {input.alignment} \
            --root {ROOTING} \
            --keep-polytomies \
            --divergence-unit mutations-per-site \
            --output-node-data {output.node_data} \
            --output-tree {output.tree}
            """


rule ancestral:
    input:
        tree=TREE,
        alignment=rules.exclude.output.filtered_sequences,
        annotation=GENBANK_PATH,
    output:
        node_data="results/muts.json",
        ancestral_sequences="results/ancestral_sequences.fasta",
    params:
        translation_template=r"results/translations/cds_%GENE.translation.fasta",
        output_translation_template=r"results/translations/cds_%GENE.ancestral.fasta",
        genes=" ".join(GENES),
        root = "results/ancestral_sequences_star.fasta" if "{STAR_ROOT}"==True else "results/ancestral_sequences.fasta",
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --annotation {input.annotation} \
            --root-sequence {input.annotation} \
            --genes {params.genes} \
            --translations {params.translation_template} \
            --output-node-data {output.node_data} \
            --output-translations {params.output_translation_template}\
            --output-sequences {output.ancestral_sequences}
        """


rule clades:
    input:
        tree=TREE,
        mutations = rules.ancestral.output.node_data,
        clades = CLADES
    output:
        json = "results/clades.json",
    shell:
        """
        augur clades --tree {input.tree} \
            --mutations {input.mutations} \
            --clades {input.clades} \
            --output-node-data {output}
        """

rule recombinant_clades:
    """
    Mark accessions listed in resources/recombinants.tsv as 'recombinant' in the
    augur clades JSON produced by the clades rule. Writes a new clades JSON
    for downstream export.
    """
    input:
        clades=rules.clades.output,
        recombinants="resources/recombinants.tsv",  # one-column TSV or header+column
    output:
        node_data="results/clades_recombinant.json",
    run:
        import json
        import pandas as pd
        import os

        # load recombinants (assume first column is accession)
        df = pd.read_csv(input.recombinants, sep="\t", dtype=str, header=0)
        recomb = set(df["accession"].to_list())

        with open(input.clades[0], "r") as fh:
            clades = json.load(fh)

        if not recomb:
            print("No recombinants loaded; writing input clades file unchanged.")
            with open(output.node_data, "w") as fh:
                json.dump(clades, fh, indent=2)
            return

        # clades JSON is expected to have structure {"nodes": {node_name: {...}, ...}, ...}
        nodes = clades.get("nodes", {})
        changed = 0
        missing = []
        for acc in sorted(recomb):
            if acc in nodes:
                node = nodes[acc]
                node["clade_membership"] = "recombinant"
                node.pop("clade_annotation", None)
                changed += 1
            else:
                missing.append(acc)

        print(f"Marked {changed} nodes as recombinant (from {len(recomb)} accessions provided).")
        if missing:
            print(f"{len(missing)} accessions not found in clades JSON; first few: {missing[:10]}")

        with open(output.node_data, "w") as fh:
            json.dump(clades, fh, indent=2)


rule get_dates:
    """Create ordering for color assignment"""
    input:
        metadata = rules.exclude.output.filtered_metadata
    output:
        ordering = "results/color_ordering.tsv"
    run:
        import pandas as pd
        column = "date"
        meta = pd.read_csv(input.metadata, delimiter='\t')

        if column not in meta.columns:
            print(f"The column '{column}' does not exist in the file.")
            sys.exit(1)

        deflist = meta[column].dropna().tolist()
        # Store unique values (ordered)
        deflist = sorted(set(deflist))
        if "XXXX-XX-XX" in deflist:
            deflist.remove("XXXX-XX-XX")

        result_df = pd.DataFrame({
            'column': ['date'] * len(deflist),
            'value': deflist
        })

        result_df.to_csv(output.ordering, sep='\t', index=False, header=False)
        
rule epitopes:
    input:
        anc_seqs = rules.ancestral.output.node_data,
        tree=TREE,
    output:
        node_data = "results/epitopes.json"
    params:
        translation = "results/translations/cds_VP1.ancestral.fasta",
        epitopes = {
        'BC':     list(range(95, 107)),         # Huang et al., 2015; Foo et al., 2008; structural mapping of neutralizing antibodies
        'DE':     list(range(142, 152)),        # Liu et al., 2011; Zaini et al., 2012.
        'EF':     list(range(165, 173)),        # Lyu et al., 2014; Wang et al., 2010.
        'CTERM':  list(range(281, 291)),        # Chang et al., 2012 (monoclonal antibody studies), structural models.
        'GH':     list(range(209, 224)),        # mutations at S215, K218 have been noted to impact neutralization. Often used in vaccine design (e.g., in VLPs or epitope grafting studies).      
        'Esc_CHN':[283, 293]},                  # Escape Mutations in C-Terminal in Chinese Samples
        min_count = 6 # number of sequences?
    run:
        import json
        from collections import defaultdict
        from Bio import Phylo, SeqIO

        manyXList = ["XXXXXXXXXXXX", "KEXXXXXXXXXX", "KERANXXXXXXX", "KERXXXXXXXXX", "KERAXXXXXXXX"]
        valid_esc_chn = {"SA", "TA", "TX", "TS", "SS", "XX"}  # Set of valid values for Esc_CHN
        # with open(input.anc_seqs) as fh:
        #     anc = json.load(fh)["nodes"]

        # Read translation files
        vp1_anc = SeqIO.to_dict(SeqIO.parse(params.translation, "fasta"))

        T = Phylo.read(input.tree, 'newick')
        for node in T.find_clades(order='preorder'):
            for child in node:
                child.parent = node

        nodes = {}
        epitope_counts = {epi: defaultdict(int) for epi in params.epitopes}

        for node in T.find_clades(order='preorder'):
            n = node.name
            aa = vp1_anc[n].seq
            nodes[n] = {}
            for epi,pos in params.epitopes.items():
                pos = [p - 1 for p in pos]  # Convert to 0-based indexing
                nodes[n][epi] = "".join([aa[p] for p in pos])
                if epi == 'CTERM':
                    if nodes[n]['CTERM'] in manyXList:
                        nodes[n]['CTERM'] = "many x"
                    elif 'X' in nodes[n]['CTERM']:
                        nodes[n]['CTERM'] = nodes[node.parent.name]['CTERM']
                if epi == 'Esc_CHN':
                    if nodes[n]['Esc_CHN'] not in valid_esc_chn:
                        nodes[n]['Esc_CHN'] = "other"
                if not n.startswith('NODE_'):
                    epitope_counts[epi][nodes[n][epi]] += 1

        for node in nodes:
            for epi,seq in nodes[node].items():
                min_count2 = params.min_count if epi != "CTERM" else 6
                if epi == "CTERM" and seq in manyXList:
                    nodes[node][epi]='many X'
                elif epitope_counts[epi][seq]<min_count2:#params.min_count:
                    nodes[node][epi]='other'

        with open(output.node_data, 'w') as fh:
            json.dump({"epitopes": params.epitopes, "nodes":nodes}, fh)


rule colors:
    """Assign colors based on ordering"""
    input:
        ordering=rules.get_dates.output.ordering,
        color_schemes=COLORS_SCHEMES,
        colors=COLORS,
    output:
        colors="results/colors_dates.tsv",
        final_colors="results/final_colors.tsv"
    shell:
        """
        python3 scripts/assign-colors.py \
            --ordering {input.ordering} \
            --color-schemes {input.color_schemes} \
            --output {output.colors}

        echo -e '\ndate\tXXXX-XX-XX\t#a6acaf' >> {output.colors}

        cat {output.colors} {input.colors} >> {output.final_colors}
        """

rule export: 
    input:
        tree=TREE,
        metadata = rules.exclude.output.filtered_metadata,
        mutations = rules.ancestral.output.node_data,
        branch_lengths = rules.refine.output.node_data,
        clades = rules.clades.output.json, # dummy_clades if not set yet
        auspice_config = AUSPICE_CONFIG,
        colors = rules.colors.output.final_colors,
        epitopes = rules.epitopes.output.node_data,
        lat_long = "resources/lat_longs.tsv",
        recombinants = rules.recombinant_clades.output.node_data
    params:
        strain_id_field = ID_FIELD,
    output:
        auspice = "results/auspice.json",
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id_field} \
            --auspice-config {input.auspice_config} \
            --lat-longs {input.lat_long} \
            --node-data {input.mutations} {input.branch_lengths} {input.clades} {input.recombinants}  \
            --colors {input.colors} \
            --output {output.auspice}
        """
        # {input.epitopes}

rule extract_clades_tsv:
    input:
        json=rules.clades.output.json,
    output:
        tsv = "results/clades_metadata.tsv"
    run:
        import json
        import csv

        with open(input.json) as f:
            data = json.load(f)

        nodes = data.get("nodes", {})

        with open(output.tsv, "w", newline="") as out_f:
            writer = csv.writer(out_f, delimiter="\t")
            writer.writerow(["accession", "clade"])

            for accession, values in nodes.items():
                clade = values.get("clade_membership", None)
                if clade:
                    writer.writerow([accession, clade])


rule subsample_example_sequences:
    input:
        all_sequences = INFERRED_SEQ_PATH,
        metadata = INFERRED_META_PATH,
        exclude = EXCLUDE,
        outliers = rules.get_outliers.output.outliers,
        incl_examples = INCLUDE_EXAMPLES,
        clades =  rules.extract_clades_tsv.output.tsv,
        tree_strains = "results/tree_strains.txt",  # strains in the tree
    output:
        example_sequences = "results/example_sequences.fasta",
    params:
        strain_id_field = ID_FIELD,
    shell:
        """
        augur merge \
            --metadata metadata={input.metadata} clades={input.clades} \
            --metadata-id-columns {params.strain_id_field} \
            --output-metadata metadata.tmp
        augur filter \
            --sequences {input.all_sequences} \
            --metadata metadata.tmp \
            --metadata-id-columns {params.strain_id_field} \
            --min-date 2010 --group-by clade \
            --subsample-max-sequences 30  \
            --min-length 4000 \
            --include {input.incl_examples} \
            --exclude {input.exclude} {input.outliers} \
            --exclude-where "clade=E" "clade=F" "clade=A" \
            --exclude-ambiguous-dates-by year \
            --probabilistic-sampling \
            --output-sequences {output.example_sequences}
        rm metadata.tmp
        """
        # seqkit grep -v -f {input.tree_strains} {input.all_sequences} \
        # | seqkit sample -n 100 -s 41 > {output.example_sequences}

rule assemble_dataset:
    input:
        tree = rules.export.output.auspice,
        reference = REFERENCE_PATH,
        annotation = GFF_PATH,
        sequences = rules.subsample_example_sequences.output.example_sequences,
        pathogen = PATHOGEN_JSON,
        readme = README_PATH,
        changelog = CHANGELOG_PATH,
    params:
        pathogen = "out-dataset/pathogen.json",
    output:
        tree = "out-dataset/tree.json",
        reference = "out-dataset/reference.fasta",
        annotation = "out-dataset/genome_annotation.gff3",
        sequences = "out-dataset/sequences.fasta",
        readme = "out-dataset/README.md",
        changelog = "out-dataset/CHANGELOG.md",
        dataset_zip = "dataset.zip",
    shell:
        """
        cp {input.tree} {output.tree}
        cp {input.reference} {output.reference}
        cp {input.annotation} {output.annotation}
        cp {input.sequences} {output.sequences}
        cp {input.pathogen} {params.pathogen}
        cp {input.readme} {output.readme}
        cp {input.changelog} {output.changelog}
        zip -rj dataset.zip  out-dataset/*
        """


rule test:
    input:
        dataset = rules.assemble_dataset.output.dataset_zip,
        sequences = rules.assemble_dataset.output.sequences,
    output:
        output = directory("test_out"),
    shell:
        """
        nextclade3 run \
            --input-dataset {input.dataset} \
            --output-all {output.output} \
            {input.sequences}
        """

rule mutLabels:
    input:
        table = "results/nextclade.tsv",
        clade = "results/clades_metadata.tsv",
        json = PATHOGEN_JSON,
    params:
        min_proportion = 0.2,
        high_threshold_proportion = 0.70,
        clades_high_threshold = ["A","B1","B2","B3","C1", "C2","C3","C4"],
        clades_to_drop = ["unassigned"],
    output:
        clade_meta = "results/clades_mut_metadata.tsv",
        properties = "results/virus_properties.json",
        json = "out-dataset/pathogen.json"
    shell:
        """
        augur merge \
            --metadata meta={input.table} clade={input.clade} \
            --metadata-id-columns meta=seqName clade=accession \
            --output-metadata {output.clade_meta}

        python3 scripts/generate_virus_properties.py \
            --clade_meta {output.clade_meta} \
            --properties {output.properties} \
            --min-prop {params.min_proportion} \
            --high-min-prop {params.high_threshold_proportion} \
            --high-prop-clades "{params.clades_high_threshold}" \
            --exclude-clades "{params.clades_to_drop}"


        jq --slurpfile v {output.properties} \
           '.mutLabels.nucMutLabelMap = $v[0].nucMutLabelMap |
            .mutLabels.nucMutLabelMapReverse = $v[0].nucMutLabelMapReverse' \
           {input.json} > {output.json}

        zip -rj dataset.zip  out-dataset/*
        """

rule fragment_testing:
    input:
        nextstrain = "testing/nextstrain_a71_vp1.tsv",
        sequences = "results/aligned.fasta",
    output:
        fragments = "testing/EV-A71_fragments.fasta"
    params:
        length = range(100, 3000, 100),  # lengths from 200 to 3000
        gene = ["VP1", "3D"]  # genes to sample from; atm only VP1 and 3D supported
    run:
        import os
        import random
        from Bio import SeqIO
        import pandas as pd

        # Read all sequences from the input file
        records = list(SeqIO.parse(input.sequences, "fasta"))
        os.makedirs(os.path.dirname(output.fragments), exist_ok=True)

        # filter records in nextstrain file
        ns_ids = list(pd.read_csv(input.nextstrain).accession)
        records = [r for r in records if r.id in ns_ids]

        with open(output.fragments, "w") as out_handle:
            for length in params.length:
                record = random.choice(records)
                seq_len = len(record.seq)
                if "VP1" in params.gene or "3D" in params.gene:
                    if "VP1" in params.gene: 
                        seq1 = record.seq[2389:3315]
                        l = len(seq1) - seq1.count("-") - seq1.count("N")
                        if l > length:
                            s = random.randint(0, l - length)
                            seq1 = seq1[s:s+length]
                            header = f"{record.id}_partial_{length}_VP1"
                            out_handle.write(f">{header}\n{seq1}\n")
                    if "3D" in params.gene:
                        seq2 = record.seq[5926:7296]
                        l = len(seq2)
                        if l > length:
                            s = random.randint(0, l - length)
                            seq2 = seq2[s:s+length]
                            header = f"{record.id}_partial_{length}_3D"
                            out_handle.write(f">{header}\n{seq2}\n")
                else: 
                    print(f"Gene {params.gene} not recognized.")
                        
                while seq_len < length:
                    record = random.choice(records)
                    seq_len = len(record.seq)
                start = random.randint(0, seq_len - length)
                fragment_seq = record.seq[start:start+length]
                header = f"{record.id}_partial_{length}"
                out_handle.write(f">{header}\n{fragment_seq}\n")


rule recombinant_testing:
    input:
        sequences = SEQUENCES,
        nextstrain = "testing/nextstrain_a71_vp1.tsv",
        clades = "results/clades_metadata.tsv",
        evD_seq = "testing/EV-D_sequence.fasta"
    output:
        recombinants = "testing/EV-A71_recombinants.fasta"
    params:
        inter_recombinants = 10,
        intra_recombinants = 10,
        min_length = 3500,
    run:
        import random
        from Bio import SeqIO
        import pandas as pd

        def eligible(records, ml):
            return [r for r in records if len(r) >= ml]

        # Load sequences and filter by Nextstrain IDs & min_length
        seqs = list(SeqIO.parse(input.sequences, "fasta"))
        ns_ids = list(pd.read_csv(input.nextstrain, sep="\t").accession)
        
        seqs = eligible([r for r in seqs if r.id in ns_ids], params.min_length)

        # Map clade assignments
        clade_map = pd.read_csv(input.clades, sep="\t").set_index("accession")["clade"].to_dict()
        clade2seqs = {}
        for r in seqs:
            clade = clade_map.get(r.id, "NA")
            clade2seqs.setdefault(clade, []).append(r)
        clades = [c for c in clade2seqs if c != "NA" and len(clade2seqs[c]) > 0]

        # EV-D sequences for intertypic recombination
        evd = eligible(list(SeqIO.parse(input.evD_seq, "fasta")), params.min_length)

        with open(output.recombinants, "w") as out:
            # Intra-typic: between clades
            for i in range(params.intra_recombinants):
                c1, c2 = random.sample(clades, 2)
                p1, p2 = random.choice(clade2seqs[c1]), random.choice(clade2seqs[c2])
                minlen = min(len(p1.seq), len(p2.seq))
                if minlen < params.min_length: continue
                x = random.randint(1, minlen-1)
                out.write(f">intra_{p1.id}_{c1}_{x}_{p2.id}_{c2}\n{p1.seq[:x]}{p2.seq[x:]}\n")

            # Inter-typic: A71 x EV-D
            for i in range(params.inter_recombinants):
                p1 = random.choice(seqs)
                p2 = random.choice(evd)
                minlen = min(len(p1.seq), len(p2.seq))
                if minlen < params.min_length: continue
                x = random.randint(1, minlen-1)
                out.write(f">inter_{p1.id}_A71_{x}_{p2.id}_D\n{p1.seq[:x]}{p2.seq[x:]}\n")

rule clean:
    shell:
        """
        rm ingest/data/* data/*
        rm -r results out-dataset test_out dataset.zip tmp
        """