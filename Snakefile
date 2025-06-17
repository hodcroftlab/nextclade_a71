# Set the parameters
REFERENCE_ACCESSION =   "U22521"
TAXON_ID =              39054
GENES =                 ["VP4", "VP2", "VP3", "VP1", "2A", "2B", "2C", "3A", "3B", "3C", "3D"]
ALLOWED_DIVERGENCE =    "3000" # was 
MIN_DATE =              "1960-01-01"
MIN_LENGTH =            "6000" # was 6000 for whole genome build on Nextstrain
MAX_SEQS =              "2000" # tree will be subsampled
ROOTING =               "AB575911 KF501389" #"mid_point" or alternative root using outgroup, e.g. the reference U22521
ID_FIELD=               "accession" # either accession or strain, used for meta-id-column in augur

# Set the paths
GFF_PATH =              "dataset/genome_annotation.gff3"
PATHOGEN_JSON =         "dataset/pathogen.json"
GENBANK_PATH =          "resources/reference.gbk"
REFERENCE_PATH =        "dataset/reference.fasta"
README_PATH =           "dataset/README.md"
CHANGELOG_PATH =        "dataset/CHANGELOG.md"
AUSPICE_CONFIG =        "resources/auspice_config.json"
EXCLUDE =               "resources/exclude.txt"
SEQUENCES =             "data/sequences.fasta"
METADATA =              "data/metadata.tsv"
CLADES =                "resources/clades.tsv"
ACCESSION_STRAIN =      "resources/accession_strain.tsv"
EXTRA_META =            "resources/meta_public.tsv"
INCLUDE_EXAMPLES =      "resources/include_examples.txt"
REFINE_DROP =          "resources/dropped_refine.txt"

FETCH_SEQUENCES = True
STAR_ROOT = True

rule all:
    input:
        auspice = "results/auspice.json",
        augur_jsons = "test_out/",
        data = "dataset.zip",
        seqs = "results/example_sequences.fasta",


if FETCH_SEQUENCES == True:
    rule fetch:
        input:
            dir = "ingest"
        output:
            sequences=SEQUENCES,
            metadata=METADATA
        shell:
            """
            cd {input.dir} 
            snakemake --cores 9 all
            cd ../
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
        echo "{REFERENCE_ACCESSION}" >> results/include.txt
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
        date_format = ['%Y-%m-%d','%Y','XX-%m-%Y','%Y-%m-%dT%H:%M:%SZ','201X-XX-XX', 'XX-XX-%Y', 'XX-XX-XXXX','%m.%Y','%m-%Y', '%d.%m.%Y', "%b-%Y", "%d-%b-%Y"],  # Date format for metadata
        date_fields = ["collection_date","date"]
    output:
        metadata = "results/metadata.tsv",  # Final output file for publications metadata
    shell:
        """
        augur merge --metadata metadata={input.meta} strains={input.strains}\
            --metadata-id-columns {params.strain_id_field} \
            --output-metadata metadata.tmp
        augur merge --metadata metadata=metadata.tmp public={input.public}\
            --metadata-id-columns {params.strain_id_field} \
            --output-metadata metadata2.tmp
        
        augur curate normalize-strings \
            --metadata metadata2.tmp \
            --id-column {params.strain_id_field} \
            --output-metadata metadata3.tmp

        augur curate format-dates \
            --metadata metadata3.tmp \
            --id-column {params.strain_id_field} \
            --expected-date-formats {params.date_format} \
            --date-fields {params.date_fields} \
            --output-metadata {output.metadata}

        rm metadata*.tmp
        """


rule index_sequences:
    message:
        """
        Creating an index of sequence composition for filtering
        """
    input:
        sequences = SEQUENCES,
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
        sequences = SEQUENCES,
        sequence_index = rules.index_sequences.output.sequence_index,
        metadata = rules.curate.output.metadata,
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
        #high-diversity 
        penalty_gap_extend = 1, #make longer gaps more costly - default is 0
        penalty_gap_open = 13,  #make gaps more expensive relative to mismatches - default is 13
        penalty_gap_open_in_frame = 18, #make gaps more expensive relative to mismatches - default is 7
        penalty_gap_open_out_of_frame = 23, #make out of frame gaps more expensive - default is 8 # prev was 19
        kmer_length = 6, #reduce to find more matches - default is 10
        kmer_distance = 25, #reduce to try more seeds - default is 50
        min_match_length = 30, #reduce to keep more seeds - default is 40
        allowed_mismatches = 15, #increase to keep more seeds - default is 8
        min_length = 30, # min_length - default is 100
        #cost of a mutation is 4
    shell:
        """
        nextclade3 run \
        -j {threads} \
        {input.sequences} \
        --input-ref {input.reference} \
        --input-annotation {input.annotation} \
        --penalty-gap-open {params.penalty_gap_open} \
        --penalty-gap-extend {params.penalty_gap_extend} \
        --penalty-gap-open-in-frame {params.penalty_gap_open_in_frame} \
        --penalty-gap-open-out-of-frame {params.penalty_gap_open_out_of_frame} \
        --kmer-length {params.kmer_length} \
        --kmer-distance {params.kmer_distance} \
        --min-match-length {params.min_match_length} \
        --allowed-mismatches {params.allowed_mismatches} \
        --min-length {params.min_length} \
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
        metadata = rules.curate.output.metadata,
        exclude = EXCLUDE,
        outliers = rules.get_outliers.output.outliers,
        examples = INCLUDE_EXAMPLES,
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
            --exclude {input.exclude} {input.outliers}  \
            --output-sequences {output.filtered_sequences} \
            --output-metadata {output.filtered_metadata} \
            --output-strains {output.strains}
        """
        #{input.examples}


rule tree:
    message:
        """
        Creating a maximum likelihood tree
        """
    input:
        alignment = rules.exclude.output.filtered_sequences,
    output:
        tree = "results/tree_raw.nwk",
    threads: 9
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

if STAR_ROOT == True:
    rule star_like_rooting:
        input:
            tree=rules.refine.output.tree,
            clades="resources/clade_map.tsv"
        output:
            tree="results/star_tree.nwk"
        params:
            strain_id_field=ID_FIELD
        script:
            "scripts/star_like_rooting.py"


rule ancestral:
    input:
        tree = rules.star_like_rooting.output.tree if STAR_ROOT == True else rules.refine.output.tree,
        alignment=rules.exclude.output.filtered_sequences,
        annotation=GENBANK_PATH,
    output:
        node_data="results/muts.json",
    params:
        translation_template=r"results/translations/cds_%GENE.translation.fasta",
        output_translation_template=r"results/translations/cds_%GENE.ancestral.fasta",
        genes=" ".join(GENES),
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
            --output-translations {params.output_translation_template}
        """

rule clades:
    input:
        tree = rules.star_like_rooting.output.tree if STAR_ROOT == True else rules.refine.output.tree,
        mutations = rules.ancestral.output.node_data,
        clades = CLADES
    output:
        json = "results/clades.json",
    shell:
        """
        augur clades --tree {input.tree} \
            --mutations {input.mutations} \
            --clades {input.clades} \
            --output-node-data {output.json}
        """

rule export:
    input:
        tree = rules.star_like_rooting.output.tree if STAR_ROOT == True else rules.refine.output.tree,
        metadata = rules.exclude.output.filtered_metadata,
        mutations = rules.ancestral.output.node_data,
        branch_lengths = rules.refine.output.node_data,
        clades = rules.clades.output.json,
        auspice_config = AUSPICE_CONFIG,
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
            --node-data {input.mutations} {input.branch_lengths} {input.clades} \
            --output {output.auspice}
        """


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
        all_sequences = SEQUENCES,
        metadata = rules.curate.output.metadata,
        exclude = EXCLUDE,
        refine = REFINE_DROP,
        outliers = rules.get_outliers.output.outliers,
        incl_examples = INCLUDE_EXAMPLES,
        clades =  rules.extract_clades_tsv.output.tsv,
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
            --min-date 2000 --group-by year clade \
            --subsample-max-sequences 25  \
            --min-length 4000 \
            --include {input.incl_examples} \
            --exclude {input.exclude} {input.outliers} {input.refine} \
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
    output:
        tree = "out-dataset/tree.json",
        reference = "out-dataset/reference.fasta",
        annotation = "out-dataset/genome_annotation.gff3",
        sequences = "out-dataset/sequences.fasta",
        pathogen = "out-dataset/pathogen.json",
        readme = "out-dataset/README.md",
        changelog = "out-dataset/CHANGELOG.md",
        dataset_zip = "dataset.zip",
    shell:
        """
        cp {input.tree} {output.tree}
        cp {input.reference} {output.reference}
        cp {input.annotation} {output.annotation}
        cp {input.sequences} {output.sequences}
        cp {input.pathogen} {output.pathogen}
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

rule clean:
    shell:
        """
        rm ingest/data/*
        rm data/*
        rm -r results
        rm -r out-dataset/
        rm -r test_out/
        rm -r dataset.zip
        """