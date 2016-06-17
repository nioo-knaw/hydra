from snakemake.utils import R

configfile: "config.json"

PROJECT = config["project"] + "/"

rule final:
    input: expand("{project}/fastqc_raw/{data}_R1_fastqc.zip \
                   {project}/fastqc_pandaseq/{data}_fastqc.zip \
                   {project}/{prog}/clst/{ds}.minsize{minsize}.usearch_smallmem.fasta \
                   {project}/{prog}/sina/{ds}.minsize{minsize}.{clmethod}.sina.taxonomy \
                   {project}/{prog}/{ds}.minsize{minsize}.{clmethod}.taxonomy.sina.biom \
                   {project}/{prog}/{ds}.minsize{minsize}.{clmethod}.tre ".split(),data=config["data"],project=config['project'],prog=["vsearch"],ds=config['project'],minsize=2,clmethod="usearch_smallmem") 

rule unpack_and_rename:
    input:
        forward = lambda wildcards: config["data"][wildcards.data]['forward'],
        reverse = lambda wildcards: config["data"][wildcards.data]['reverse']
    output:
        forward="{project}/gunzip/{data}_R1.fastq",
        reverse="{project}/gunzip/{data}_R2.fastq"
    params:
        prefix="{data}"
    threads: 2
    run: 
        shell("zcat {input.forward} | awk '{{print (NR%4 == 1) ? \"@{params.prefix}_\" substr($0,2) : $0}}' > {output.forward}")
        shell("zcat {input.reverse} | awk '{{print (NR%4 == 1) ? \"@{params.prefix}_\" substr($0,2) : $0}}' > {output.reverse}")

rule fastqc:
    input:
        forward="{project}/gunzip/{data}_R1.fastq",
        reverse="{project}/gunzip/{data}_R2.fastq"
    output:
        zip="{project}/fastqc_raw/{data}_R1_fastqc.zip"
    params:
        dir="{project}/fastqc_raw",
        adapters = config["adapters_fasta"]
    log: "fastqc_raw.log"
    threads: 2
    run:
        shell("fastqc -q -t {threads} --contaminants {params.adapters} --outdir {params.dir} {input.forward} > {params.dir}/{log}")
        shell("fastqc -q -t {threads} --contaminants {params.adapters} --outdir {params.dir} {input.reverse} > {params.dir}/{log}")

rule pandaseq:
    input:      
        forward="{project}/gunzip/{data}_R1.fastq",
        reverse="{project}/gunzip/{data}_R2.fastq"
    output:
        fastq = "{project}/pandaseq/{data}.fastq"
    params:
        overlap = config['pandaseq_overlap'],
        quality = config['pandaseq_quality'],
        minlength = config['pandaseq_minlength'],
        maxlength = config['pandaseq_maxlength'],
        forward_primer = config['forward_primer'],
        reverse_primer = config['reverse_primer']
    log: "{project}/pandaseq/{data}_pandaseq.stdout"
    threads: 1
    #shell: "source /data/tools/RDP_Assembler/1.0.3/env.sh; pandaseq -N -o {params.overlap} -e {params.quality} -F -d rbfkms -l {params.minlength} -L {params.maxlength} -T {threads} -f {input.forward} -r {input.reverse}  1> {output.fastq} 2> {log}"
    shell: "/data/tools/pandaseq/2.9/bin/pandaseq -N -f {input.forward} -r {input.reverse} -p {params.forward_primer} -q {params.reverse_primer} -A rdp_mle -T {threads} -w {output.fastq} -g {log}"

rule fastqc_pandaseq:
    input:
        fastq = "{project}/pandaseq/{data}.fastq"
    output: 
        dir="{project}/pandaseq/{data}_fastqc/",zip="{project}/fastqc_pandaseq/{data}_fastqc.zip"
    params:
        dir="{project}/fastqc_pandaseq",
        adapters = config["adapters_fasta"]
    log: "fastqc.log"
    threads: 8
    shell: "fastqc -q -t {threads} --contaminants {params.adapters} --outdir {params.dir} {input.fastq} > {params.dir}/{log}"


rule primer_matching:
    input:
        "{project}/pandaseq/{data}.fastq"
    output:
        "{project}/flexbar/{data}.fastq"
    params:
        prefix_forward="flexbar_{data}_forward",
        prefix_reverse="flexbar_{data}_reverse"
    log: "{project}/trim/flexbar_{data}.log"
    run:
        shell("/data/tools/flexbar/2.5/flexbar -t {params.prefix_forward} -b primers.fasta -r {input} --barcode-min-overlap 10 --barcode-threshold 3 --min-read-length 50 >> {log}")
        shell("cat {params.prefix_forward}* | fastx_reverse_complement | /data/tools/flexbar/2.5/flexbar -t {params.prefix_reverse} -b primers.fasta --reads - --barcode-trim-end LEFT --barcode-min-overlap 10 --barcode-threshold 3 --min-read-length 50 --barcode-unassigned >> {log}")
        shell("cat {params.prefix_reverse}* > {output}")


rule fastq2fasta:
    input:
        fastq = "{project}/pandaseq/{data}.fastq"
    output:
        fasta = "{project}/pandaseq/{data}.fasta"
    shell: "fastq_to_fasta -i {input.fastq} -o {output.fasta}"

# Combine per sample files to a single project file
rule mergefiles:
    input:
        fasta = expand(PROJECT + "pandaseq/{data}.fasta", data=config["data"]),
    output: 
        fasta="{project}/mergefiles/{project}.fasta"
    params:
        samples=config["data"]
    shell: """cat {input}  > {output}"""

########
# VSEARCH PIPELINE
########
USEARCH = "/data/tools/usearch/7.0.1090/usearch"
VSEARCH = "/data/tools/vsearch/1.0.10/vsearch-1.0.10-linux-x86_64"
# VSEARCH commands are compatible with USEARCH commands
# depending on the output file that is requested the required program is chosen:
# Example: snakemake -j 8  vsearch/B2R.results.txt -p

# Dereplication
rule derep:
    input:
        "{project}/mergefiles/{ds}.fasta",
    output:
        temp("{project}/{prog}/{ds}.derep.fasta")
    threads: 8
    run:
        cmd = ""
        if wildcards.prog == "vsearch":
            cmd = VSEARCH
        elif wildcards.prog == "usearch":
            cmd = USEARCH
        shell("{cmd} -derep_fulllength {input} -output {output} -sizeout -threads {threads}")

# Abundance sort and discard singletons
rule sortbysize:
    input:
        "{project}/{prog}/{ds}.derep.fasta"
    output:
        temp("{project}/{prog}/{ds}.sorted.minsize{minsize}.fasta")
    params:
        minsize="{minsize}"
    threads: 8
    run:
        cmd = ""
        if wildcards.prog == "vsearch":
            cmd = VSEARCH
            shell("{cmd} -sortbysize {input} -fasta_width 0 -output {output} -threads {threads} -minsize {params.minsize}")      
        elif wildcards.prog == "usearch":
            cmd = USEARCH
            shell("{cmd} -sortbysize {input} -output {output} -minsize {params.minsize}")

# Uclust clustering
rule smallmem:
    input:
        "{project}/{prog}/{ds}.sorted.minsize{minsize}.fasta"
    output:
        otus=protected("{project}/{prog}/clst/{ds}.minsize{minsize}.usearch_smallmem.fasta")
    threads: 8
    run:
        cmd = ""
        if wildcards.prog == "vsearch":
            cmd = VSEARCH      
        elif wildcards.prog == "usearch":
            cmd = USEARCH
        shell("{cmd} --cluster_smallmem {input} --usersort -centroids {output.otus} --id 0.97 -sizeout")

rule cluster_fast:
    input:   
        "{project}/{prog}/{ds}.sorted.minsize{minsize}.fasta"
    output:
        otus=protected("{project}/{prog}/{ds}.minsize{minsize}.usearch_cluster_fast.fasta")
    threads: 8
    run:
        cmd = ""
        if wildcards.prog == "vsearch":
            cmd = VSEARCH
        elif wildcards.prog == "usearch":
            cmd = USEARCH
        shell("{cmd} --cluster_fast {input} --usersort -centroids {output.otus} --id 0.97 -sizeout")

rule uparse:
    input:
        "{project}/{prog}/{ds}.sorted.minsize{minsize}.fasta"
    output:
        otus=protected("{project}/{prog}/{ds}.minsize{minsize}.uparse.fasta")
    shell: "{USEARCH} -cluster_otus {input} -otus {output.otus} -relabel OTU_ -sizeout"


#
# Chimera checking
#

rule uchime:
    input:
        "{project}/{prog}/clst/{ds}.minsize{minsize}.{clmethod}.fasta"
    output:
        chimeras="{project}/{prog}/uchime/{ds}.minsize{minsize}.{clmethod}.chimeras",
        nonchimeras="{project}/{prog}/uchime/{ds}.minsize{minsize}.{clmethod}.fasta"
    log: "{project}/{prog}/uchime/{ds}.minsize{minsize}.{clmethod}.uchime.log"
    run:
        cmd = ""
        if wildcards.prog == "vsearch":
            cmd = VSEARCH
        elif wildcards.prog == "usearch":
            cmd = USEARCH
        shell("{cmd} --uchime_denovo {input} --nonchimeras {output.nonchimeras} --chimeras {output.chimeras} > {log}")

# 
# Mapping
#
#TODO Check mapping accuracy!!!
rule make_otu_names:
    input:
        "{project}/{prog}/uchime/{ds}.minsize{minsize}.{clmethod}.fasta"
    output:
        "{project}/{prog}/otus/{ds}.minsize{minsize}.{clmethod}.fasta"
    shell: "python2.7 /data/tools/usearch/uparse_scripts/fasta_number.py {input} OTU_ > {output}"

rule mapping:
    input:
        otus="{project}/{prog}/otus/{ds}.minsize{minsize}.{clmethod}.fasta",
        reads="{project}/mergefiles/{ds}.fasta"
    output:
        "{project}/{prog}/otus/{ds}.minsize{minsize}.{clmethod}.uc"
    run:
        cmd = ""
        if wildcards.prog == "vsearch":
            cmd = VSEARCH
        elif wildcards.prog == "usearch":
            cmd = USEARCH
        shell("{cmd} -usearch_global {input.reads} -db {input.otus} -strand plus -id 0.97 -uc {output}")

rule create_otutable:
    input:
        "{project}/{prog}/otus/{ds}.minsize{minsize}.{clmethod}.uc"
    output:
        "{project}/{prog}/otus/{ds}.minsize{minsize}.{clmethod}.otutable.txt"
    shell: "python2.7 /data/tools/usearch/uparse_scripts/uc2otutab.py {input} > {output}"

# convert to biom file
rule biom_otu:
    input: 
        "{project}/{prog}/otus/{ds}.minsize{minsize}.{clmethod}.otutable.txt"
    output:
        "{project}/{prog}/otus/{ds}.minsize{minsize}.{clmethod}.biom" 
    shell: "/data/tools/qiime/1.9/qiime1.9/bin/biom convert -i {input} --to-json -o {output} --table-type='OTU table'"

#
# Taxonomy
#
SINA_VERSION=config['sina_version']
SILVADB = config['silva_db']
SINA_MIN_SIM = config['sina_min_sim']
SILVA_ARB = config['silva_arb']



rule sina_parallel_edgar:
    input:
        "{project}/{prog}/otus/{ds}.minsize{minsize}.{clmethod}.fasta"
    output:
        #align="{project}/{prog}/sina/{ds}.minsize{minsize}.{clmethod}.sina.{chunk}.align",
        #aligncsv="{project}/{prog}/sina/{ds}.minsize{minsize}.{clmethod}.sina.{chunk}.align.csv",
        log="{project}/{prog}/sina/{ds}.minsize{minsize}.{clmethod}.sina.log",
        align="{project}/{prog}/sina/{ds}.minsize{minsize}.{clmethod}.sina.align",
    log: "{project}/{prog}/sina/{ds}.minsize{minsize}.{clmethod}.sina.log"
    priority: -1
    threads: 8
    # TODO: turn is set to all to get classification. Reverse the reads in earlier stage!
    shell: "cat {input} | parallel --block 1000K -j{threads} --recstart '>' --pipe /data/tools/sina/{SINA_VERSION}/sina --log-file {log} -i /dev/stdin -o {output.align} --outtype fasta --meta-fmt csv --ptdb {SILVA_ARB} --overhang remove --turn all --search --search-db {SILVA_ARB} --search-min-sim 0.95 --search-no-fast --search-kmer-len 10 --lca-fields tax_slv"

rule sina_get_taxonomy_from_logfile_edgar:
    input:
        log="{project}/{prog}/sina/{ds}.minsize{minsize}.{clmethod}.sina.log"
    output:
        taxonomy="{project}/{prog}/sina/{ds}.minsize{minsize}.{clmethod}.sina.taxonomy"
    # Parse the log file from Sina to get the taxonomic classification
    # The csv output does not contain the sequence identifier, thats way this approach is better
    # The first space needs to be replaced in order to keep the space in the taxonomy string (would be splitted otherwise)
    # Brackets are escaped by an extra bracket, because they are internaly recognised by Snakemake
    shell: "cat {input.log} | sed 's/ /|/1' | awk -F '|'  '/^sequence_identifier:/ {{id=$2}} /^lca_tax_slv:/{{split(id,a,\" \"); print a[1] \"\t\" $2}}' | tr ' ' '_' > {output.taxonomy}"

# Tree
rule filter_alignment:
    input:
        align="{project}/{prog}/sina/{ds}.minsize{minsize}.{clmethod}.sina.align"
    output:
        filtered="{project}/{prog}/sina/{ds}.minsize{minsize}.{clmethod}.sina_pfiltered.fasta"
    params:
        outdir="{project}/{prog}/sina/"
    shell: "set +u; source /data/tools/qiime/1.9/env.sh; set -u; filter_alignment.py -i {input.align} -o {params.outdir} --suppress_lane_mask_filter --entropy_threshold 0.10"

rule make_tree:
    input:
        align="{project}/{prog}/sina/{ds}.minsize{minsize}.{clmethod}.sina_pfiltered.fasta"
    output:
        tree="{project}/{prog}/{ds}.minsize{minsize}.{clmethod}.tre"
    shell: "set +u; source /data/tools/qiime/1.9/env.sh; source /data/tools/arb/6.0.1/env.sh; set -u; make_phylogeny.py -i {input.align} -t fasttree -o {output.tree}"

rule biom_tax_sina:
    input:
        taxonomy="{project}/{prog}/sina/{ds}.minsize{minsize}.{clmethod}.sina.taxonomy",
        biom="{project}/{prog}/otus/{ds}.minsize{minsize}.{clmethod}.biom",
        meta=config["metadata"]
    output:
        biom=protected("{project}/{prog}/{ds}.minsize{minsize}.{clmethod}.taxonomy.sina.biom"),
        otutable=protected("{project}/{prog}/{ds}.minsize{minsize}.{clmethod}.taxonomy.sina.otutable.txt")
    run:
        shell("/data/tools/qiime/1.9/qiime1.9/bin/biom add-metadata -i {input.biom} -o {output.biom} --output-as-json --observation-metadata-fp {input.taxonomy} --observation-header OTUID,taxonomy --sc-separated taxonomy --float-fields confidence --sample-metadata-fp {input.meta}")
        shell("/data/tools/qiime/1.9/qiime1.9/bin/biom convert --to-tsv --header-key=taxonomy -i {output.biom} -o {output.otutable}")
