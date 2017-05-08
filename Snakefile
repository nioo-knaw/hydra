from snakemake.utils import R
import os

if os.path.isfile("config.yaml"):
    configfile: "config.yaml"

PROJECT = config["project"] + "/"

rule final:
    input: expand("{project}/stats/readstat.{data}.csv \
                   {project}/{prog}/clst/{ds}.minsize{minsize}.usearch_smallmem.fasta \
                   {project}/{prog}/sina/{ds}.minsize{minsize}.{clmethod}.sina.taxonomy \
                   {project}/{prog}/{ds}.minsize{minsize}.{clmethod}.taxonomy.sina.biom \
                   {project}/{prog}/{ds}.minsize{minsize}.{clmethod}.tre".split(),data=config["data"],project=config['project'],prog=["vsearch"],ds=config['project'],minsize=2,clmethod="usearch_smallmem") 

from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider()

rule unpack_and_rename:
    input:
       forward = lambda wildcards: FTP.remote(config["data"][wildcards.data]["path"][0], keep_local=False) if config["remote"] \
else lambda wildcards: config["data"][wildcards.data]["path"][0],
       reverse = lambda wildcards: FTP.remote(config["data"][wildcards.data]["path"][1], keep_local=False) if config["remote"] \
else lambda wildcards: config["data"][wildcards.data]["path"][1]
    output:
        forward="{project}/gunzip/{data}_R1.fastq",
        reverse="{project}/gunzip/{data}_R2.fastq"
    params:
        prefix="{data}"
    threads: 1
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
    conda: "envs/fastqc.yaml"
    shell: "fastqc -q -t {threads} --contaminants {params.adapters} --outdir {params.dir} {input.forward} > {params.dir}/{log} && fastqc -q -t {threads} --contaminants {params.adapters} --outdir {params.dir} {input.reverse} > {params.dir}/{log}"

if config['mergepairs'] == 'pandaseq':
    rule pandaseq:
        input:      
            forward="{project}/gunzip/{data}_R1.fastq",
            reverse="{project}/gunzip/{data}_R2.fastq"
        output:
            fasta = "{project}/mergepairs/{data}.fasta"
        params:
            overlap = config['pandaseq_overlap'],
            quality = config['pandaseq_quality'],
            minlength = config['pandaseq_minlength'],
            maxlength = config['pandaseq_maxlength'],
            forward_primer = config['forward_primer'],
            reverse_primer = config['reverse_primer']
        log: "{project}/mergepairs/{data}_pandaseq.stdout"
        threads: 1
        conda: "envs/pandaseq.yaml"
        shell: "pandaseq -N -A rdp_mle -o {params.overlap} -l {params.minlength} -L {params.maxlength} -f {input.forward} -r {input.reverse} -T {threads} -w {output.fasta} -g {log}"

rule fastqc_pandaseq:
    input:
        fastq = "{project}/mergepairs/{data}.fastq"
    output: 
        dir="{project}/mergepairs/{data}_fastqc/",zip="{project}/fastqc_pandaseq/{data}_fastqc.zip"
    params:
        dir="{project}/fastqc_pandaseq",
        adapters = config["adapters_fasta"]
    log: "fastqc.log"
    threads: 8
    conda: "envs/fastqc.yaml"
    shell: "fastqc -q -t {threads} --contaminants {params.adapters} --outdir {params.dir} {input.fastq} > {params.dir}/{log}"

rule readstat_mergepairs:
    input:
        fasta = "{project}/mergepairs/{data}.fasta"
    output:
        temporary("{project}/stats/readstat.{data}.csv")
    log:
        "{project}/stats/readstat.{data}.log"
    conda:
        "envs/khmer.yaml"
    threads: 1
    shell: "readstats.py {input} --csv -o {output} 2> {log}"

rule readstat_all:
    input:
        expand("{project}/stats/readstat.{data}.csv", project=config['project'], data=config["data"])
    output:
        protected("{project}/stats/readstat.csv")
    shell: "cat {input[0]} | head -n 1 > {output} && for file in {input}; do tail -n +2 $file >> {output}; done;"

if config['mergepairs'] == 'vsearch':
    rule mergepairs:
        input:
            forward="{project}/gunzip/{data}_R1.fastq",
            reverse="{project}/gunzip/{data}_R2.fastq"
        output:
            fasta = "{project}/mergepairs/{data}.fasta"
        threads: 1
        conda: "envs/vsearch.yaml"
        shell: "vsearch --threads {threads} --fastq_mergepairs {input.forward} --reverse {input.reverse} --fastq_allowmergestagger --fastq_minmergelen 200 --fastaout {output}"


# Combine per sample files to a single project file
rule mergefiles:
    input:
        fasta = expand(PROJECT + "mergepairs/{data}.fasta", data=config["data"]),
    output: 
        fasta="{project}/mergefiles/{project}.fasta"
    params:
        samples=config["data"]
    shell: """cat {input}  > {output}"""

rule length:
    input:
        "{project}/mergefiles/{project}.fasta"
    output:
        protected("{project}/stats/readlength.csv")
    shell: "awk -f ../src/hydra/seqlen.awk {input} > {output}"

# Dereplication
rule derep:
    input:
        "{project}/mergefiles/{ds}.fasta",
    output:
        temp("{project}/{prog}/{ds}.derep.fasta")
    threads: 8
    conda: "envs/vsearch.yaml"
    shell: "vsearch -derep_fulllength {input} -output {output} -sizeout -threads {threads}"

# Abundance sort and discard singletons
rule sortbysize:
    input:
        "{project}/{prog}/{ds}.derep.fasta"
    output:
        temp("{project}/{prog}/{ds}.sorted.minsize{minsize}.fasta")
    params:
        minsize="{minsize}"
    threads: 8
    conda: "envs/vsearch.yaml"
    shell: "vsearch -sortbysize {input} -output {output} -minsize {params.minsize}"

# Uclust clustering
rule smallmem:
    input:
        "{project}/{prog}/{ds}.sorted.minsize{minsize}.fasta"
    output:
        otus=protected("{project}/{prog}/clst/{ds}.minsize{minsize}.usearch_smallmem.fasta")
    threads: 8
    conda: "envs/vsearch.yaml"
    shell: "vsearch --cluster_smallmem {input} --usersort -centroids {output.otus} --id 0.97 -sizeout"

rule cluster_fast:
    input:   
        "{project}/{prog}/{ds}.sorted.minsize{minsize}.fasta"
    output:
        otus=protected("{project}/{prog}/{ds}.minsize{minsize}.usearch_cluster_fast.fasta")
    threads: 8
    conda: "envs/vsearch.yaml"
    shell: "vsearch --cluster_fast {input} --usersort -centroids {output.otus} --id 0.97 -sizeout"


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
    conda: "envs/vsearch.yaml"
    shell: "vsearch --uchime_denovo {input} --nonchimeras {output.nonchimeras} --chimeras {output.chimeras} > {log}"

# 
# Mapping
#

rule make_otu_names:
    input:
        "{project}/{prog}/uchime/{ds}.minsize{minsize}.{clmethod}.fasta"
    output:
        "{project}/{prog}/otus/{ds}.minsize{minsize}.{clmethod}.fasta"
    shell: "python2.7 uparse_scripts/fasta_number.py {input} OTU_ > {output}"

rule mapping:
    input:
        otus="{project}/{prog}/otus/{ds}.minsize{minsize}.{clmethod}.fasta",
        reads="{project}/mergefiles/{ds}.fasta"
    output:
        "{project}/{prog}/otus/{ds}.minsize{minsize}.{clmethod}.uc"
    conda: "envs/vsearch.yaml"
    shell: "vsearch -usearch_global {input.reads} -db {input.otus} -strand plus -id 0.97 -uc {output}"

rule create_otutable:
    input:
        "{project}/{prog}/otus/{ds}.minsize{minsize}.{clmethod}.uc"
    output:
        "{project}/{prog}/otus/{ds}.minsize{minsize}.{clmethod}.otutable.txt"
    shell: "python2.7 uparse_scripts/uc2otutab.py {input} > {output}"

# convert to biom file
rule biom_otu:
    input: 
        "{project}/{prog}/otus/{ds}.minsize{minsize}.{clmethod}.otutable.txt"
    output:
        "{project}/{prog}/otus/{ds}.minsize{minsize}.{clmethod}.biom" 
    conda: "envs/biom-format.yaml"
    shell: "biom convert -i {input} --to-json -o {output} --table-type='OTU table'"

#
# Taxonomy
#
rule sina_parallel_edgar:
    input:
        "{project}/{prog}/otus/{ds}.minsize{minsize}.{clmethod}.fasta"
    output:
        #align="{project}/{prog}/sina/{ds}.minsize{minsize}.{clmethod}.sina.{chunk}.align",
        #aligncsv="{project}/{prog}/sina/{ds}.minsize{minsize}.{clmethod}.sina.{chunk}.align.csv",
        log="{project}/{prog}/sina/{ds}.minsize{minsize}.{clmethod}.sina.log",
        align="{project}/{prog}/sina/{ds}.minsize{minsize}.{clmethod}.sina.align",
    params:
        silva_arb = config["silva_arb"]
    log: "{project}/{prog}/sina/{ds}.minsize{minsize}.{clmethod}.sina.log"
    priority: -1
    threads: 24
    # TODO: turn is set to all to get classification. Reverse the reads in earlier stage!
#    conda: "envs/sina.yaml"
    shell: "cat {input} | parallel --block 1000K -j{threads} --recstart '>' --pipe sina --log-file {log} -i /dev/stdin --intype fasta -o {output.align} --outtype fasta --meta-fmt csv --ptdb {params.silva_arb} --overhang remove --turn all --search --search-db {params.silva_arb} --search-min-sim 0.95 --search-no-fast --search-kmer-len 10 --lca-fields tax_slv"

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
    conda: "envs/qiime.yaml"
    shell: "filter_alignment.py -i {input.align} -o {params.outdir} --suppress_lane_mask_filter --entropy_threshold 0.10"

rule make_tree:
    input:
        align="{project}/{prog}/sina/{ds}.minsize{minsize}.{clmethod}.sina_pfiltered.fasta"
    output:
        tree="{project}/{prog}/{ds}.minsize{minsize}.{clmethod}.tre"
    conda: "envs/qiime.yaml"
    shell: "make_phylogeny.py -i {input.align} -t fasttree -o {output.tree}"



rule biom_tax_sina:
    input:
        taxonomy="{project}/{prog}/sina/{ds}.minsize{minsize}.{clmethod}.sina.taxonomy",
        biom="{project}/{prog}/otus/{ds}.minsize{minsize}.{clmethod}.biom",
        meta=config["metadata"]
    output:
        taxonomy="{project}/{prog}/sina/{ds}.minsize{minsize}.{clmethod}.sina.qiimeformat.taxonomy",
        biom=protected("{project}/{prog}/{ds}.minsize{minsize}.{clmethod}.taxonomy.sina.biom"),
        otutable=protected("{project}/{prog}/{ds}.minsize{minsize}.{clmethod}.taxonomy.sina.otutable.txt")
    conda: "envs/biom-format.yaml"
    shell: """cat {input.taxonomy} | awk -F"[;\t]" 'BEGIN{{print "OTUs,Domain,Phylum,Class,Order,Family,Genus"}}{{print $1"\\tk__"$2"; p__"$3"; c__"$4"; o__"$5"; f__"$6"; g__"$7"; s__"$8}}' > {output.taxonomy} && \
              biom add-metadata -i {input.biom} -o {output.biom} --observation-metadata-fp {output.taxonomy} --observation-header OTUID,taxonomy --sc-separated taxonomy --float-fields confidence --sample-metadata-fp {input.meta} && \
              biom convert --to-tsv --header-key=taxonomy -i {output.biom} -o {output.otutable}
           """

rule report:
    input:
        readstat = "{project}/stats/readstat.csv",
        biom = expand("{{project}}/{prog}/{ds}.minsize{minsize}.{clmethod}.taxonomy.sina.biom", prog=["vsearch"],ds=config['project'],minsize=2,clmethod="usearch_smallmem")
    output:
        "{project}/stats/report.html"
    params:
        prefix="{project}/stats/report",
        mergemethod = config['mergepairs']
    conda: "envs/report.yaml"
    script:
        "report.Rmd"

