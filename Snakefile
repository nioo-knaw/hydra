from snakemake.utils import R
from snakemake.utils import min_version
import os

min_version("3.12") # R Markdown reports have been introduced in Snakemake 3.12

if os.path.isfile("config.yaml"):
    configfile: "config.yaml"

PROJECT = config["project"] + "/"

rule final:
    input: expand("{project}/stats/{data}_phix_stats.txt \
                   {project}/barcode/{data}/barcodes.fastq \
                   {project}/{prog}/clst/{ds}.minsize{minsize}.{clmethod}.fasta \
                   {project}/{prog}/{ds}.minsize{minsize}.{clmethod}.taxonomy.biom \
                   {project}/stats/readstat.{data}.csv \
                   {project}/stats/report.html".split(),data=config["data"],project=config['project'],prog=["vsearch"],ds=config['project'],minsize=2,clmethod=config['clustering']) 


rule unpack_and_rename:
    input:
        forward = lambda wildcards: config["data"][wildcards.data]["path"][0],
        reverse = lambda wildcards: config["data"][wildcards.data]["path"][1]
    output:
        forward="{project}/gunzip/{data}_R1.fastq",
        reverse="{project}/gunzip/{data}_R2.fastq"
    params:
        prefix="{data}"
    threads: 2
    run: 
        shell("zcat {input.forward} | awk '{{print (NR%4 == 1) ? \"@{params.prefix}_\" substr($0,2) : $0}}' > {output.forward}")
        shell("zcat {input.reverse} | awk '{{print (NR%4 == 1) ? \"@{params.prefix}_\" substr($0,2) : $0}}' > {output.reverse}")

rule filter_phix:
     input:
        forward="{project}/gunzip/{data}_R1.fastq",
        reverse="{project}/gunzip/{data}_R2.fastq"
     output:
        forward="{project}/filter/{data}_R1.fastq",
        reverse="{project}/filter/{data}_R2.fastq",
        stats="{project}/stats/{data}_phix_stats.txt"
     params:
         phix="/data/db/contaminants/phix/phix.fasta"
     log: "{project}/filter/{data}.log"
     conda: "envs/bbmap.yaml"
     threads: 16
     shell:"""bbduk2.sh -Xmx8g in={input.forward} in2={input.reverse} out={output.forward} out2={output.reverse} \
              fref={params.phix} qtrim="rl" trimq=30 threads={threads} stats={output.stats} 2> {log}"""

rule remove_barcodes:
    input:
        forward="{project}/filter/{data}_R1.fastq",
        reverse="{project}/filter/{data}_R2.fastq",
    output:
        barcodes=temp("{project}/barcode/{data}/barcodes.fastq"),
        barcodes_fasta=temp("{project}/barcode/{data}/barcodes.fasta"),
        forward="{project}/barcode/{data}_R1.fastq",
        reverse="{project}/barcode/{data}_R2.fastq",
        forward_unpaired="{project}/barcode/{data}_R1_unpaired.fastq",
        reverse_unpaired="{project}/barcode/{data}_R2_unpaired.fastq",
    params:
        outdir="{project}/barcode/{data}/",
        threshold=5
    log: "{project}/barcode/{data}.log"
    conda: "envs/barcode.yaml"
    threads: 8
    shell: """extract_barcodes.py -f {input.forward}  -s# -l 8 -o {params.outdir} -c barcode_in_label && fastq_to_fasta < {output.barcodes} > {output.barcodes_fasta} && \
              trimmomatic PE -threads {threads} -phred33 {input.forward} {input.reverse} {output.forward} {output.forward_unpaired} {output.reverse} {output.reverse_unpaired} ILLUMINACLIP:{output.barcodes_fasta}:0:0:{params.threshold} 2> {log}
"""

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
            forward="{project}/barcode/{data}_R1.fastq",
            reverse="{project}/barcode/{data}_R2.fastq"
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
            forward="{project}/barcode/{data}_R1.fastq",
            reverse="{project}/barcode/{data}_R2.fastq"
        output:
            fasta = "{project}/mergepairs/{data}.fasta"
        threads: 1
        conda: "envs/vsearch.yaml"
        shell: "vsearch --threads {threads} --fastq_mergepairs {input.forward} --reverse {input.reverse} --fastq_allowmergestagger --fastq_minmergelen 200 --fastaout {output}"

if config['its'] == True:
    rule extract_its:
            input:
                    fasta="{project}/mergepairs/{data}.fasta"
            output:
                    fasta="{project}/itsx/{data}.ITS2.fasta"
            params:
                    basename="{project}/itsx/{data}",
                    dir="{project}/itsx"
            log: "itsx.log"
            threads: 4
            # TODO: Filter on specific list of organisms? 
            # Only ITS2 region?
            shell: "source /data/tools/hmmer/3.0/env.sh; /data/tools/ITSx/1.0.10/ITSx --cpu {threads} --preserve TRUE -i {input.fasta} -o {params.basename} > {params.dir}/{log}"

# Combine per sample files to a single project file
rule mergefiles:
    input:
        fasta = expand(PROJECT + "itsx/{data}.ITS2.fasta", data=config["data"]) if config['its'] \
                else expand(PROJECT + "mergepairs/{data}.fasta", data=config["data"]),
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
        "{project}/mergefiles/{ds}.fasta"
    output:
        temporary("{project}/{prog}/{ds}.derep.fasta")
    threads: 8
    conda: "envs/vsearch.yaml"
    shell: "vsearch -derep_fulllength {input} -output {output} -sizeout -threads {threads}"

# Abundance sort and discard singletons
rule sortbysize:
    input:
        "{project}/{prog}/{ds}.derep.fasta"
    output:
        temporary("{project}/{prog}/{ds}.sorted.minsize{minsize}.fasta")
    params:
        minsize="{minsize}"
    threads: 8
    conda: "envs/vsearch.yaml"
    shell: "vsearch -sortbysize {input} -output {output} -minsize {params.minsize}"

if config['clustering'] == "usearch_smallmem":
    # Uclust clustering
    rule smallmem:
        input:
	        "{project}/{prog}/{ds}.sorted.minsize{minsize}.fasta"
        output:
	        otus=protected("{project}/{prog}/clst/{ds}.minsize{minsize}.usearch_smallmem.fasta")
        threads: 8
        conda: "envs/vsearch.yaml"
        shell: "vsearch --cluster_smallmem {input} --usersort -centroids {output.otus} --id 0.97 -sizeout"

#
# Swarm
#
if config['clustering'] == "swarm":

    # Swarm
    rule swarm:
        input: 
            "{project}/{prog}/{ds}.sorted.minsize{minsize}.fasta"
        output:
            swarms="{project}/{prog}/{ds}.minsize{minsize}.swarm.swarms",
            stats="{project}/{prog}/{ds}.minsize{minsize}.swarm.stats"
        params: d="1" 
        threads: 16
        conda: "envs/swarm.yaml"
        shell: "swarm -d {params.d} -t {threads} -z -u uclust.out -s {output.stats} < {input}  > {output.swarms}"

    rule swarm_get_seed:
        input: 
            swarms="{project}/{prog}/{ds}.minsize{minsize}.swarm.swarms",
            amplicons="{project}/{prog}/clst/{ds}.sorted.minsize{minsize}.fasta"
        output:
            seeds="{project}/{prog}/clst/{ds}.minsize{minsize}.swarm.fasta"
        shell: "SEEDS=$(mktemp); cut -d ' ' -f 1 {input.swarms} | sed -e 's/^/>/' > '${{SEEDS}}'; grep -A 1 -F -f '${{SEEDS}}' {input.amplicons} | sed -e '/^--$/d' > {output.seeds}"

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
if config["classification"] == "silva":
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
            taxonomy="{project}/{prog}/sina/{ds}.minsize{minsize}.{clmethod}.sina.taxonomy.txt"
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
            taxonomy="{project}/{prog}/sina/{ds}.minsize{minsize}.{clmethod}.sina.taxonomy.txt",
            biom="{project}/{prog}/otus/{ds}.minsize{minsize}.{clmethod}.biom",
            meta=config["metadata"]
        output:
            taxonomy="{project}/{prog}/sina/{ds}.minsize{minsize}.{clmethod}.sina.qiimeformat.taxonomy",
            biom=protected("{project}/{prog}/{ds}.minsize{minsize}.{clmethod}.taxonomy.biom"),
            otutable=protected("{project}/{prog}/{ds}.minsize{minsize}.{clmethod}.taxonomy.otutable.txt")
        conda: "envs/biom-format.yaml"
        shell: """cat {input.taxonomy} | awk -F"[;\t]" 'BEGIN{{print "OTUs,Domain,Phylum,Class,Order,Family,Genus"}}{{print $1"\\tk__"$2"; p__"$3"; c__"$4"; o__"$5"; f__"$6"; g__"$7"; s__"$8}}' > {output.taxonomy} && \
                  biom add-metadata -i {input.biom} -o {output.biom} --observation-metadata-fp {output.taxonomy} --observation-header OTUID,taxonomy --sc-separated taxonomy --float-fields confidence --sample-metadata-fp {input.meta} --output-as-json && \
                  biom convert --to-tsv --header-key=taxonomy -i {output.biom} -o {output.otutable}
               """

if config["classification"] == "stampa":
    rule stampa:
        input:
            "{project}/{prog}/otus/{ds}.minsize{minsize}.{clmethod}.fasta"
        output:
            swarm="{project}/{prog}/stampa/{ds}.minsize{minsize}.{clmethod}.fasta",
            hits="{project}/{prog}/stampa/hits.{ds}.minsize{minsize}.{clmethod}_usearch_global",
            results="{project}/{prog}/stampa/results.{ds}.minsize{minsize}.{clmethod}_usearch_global",
            taxonomy="{project}/{prog}/stampa/{ds}.minsize{minsize}.{clmethod}.taxonomy.txt",
        params:
             stampadir="{project}/{prog}/stampa/",
             db = config['stampa_db']
        conda: "envs/vsearch.yaml"
        threads: 16
        # Create STAMPA compatible input
        # Replace underscore in otu names and add fake abundance information
        shell:"""
            sed 's/_/:/' {input} | awk '/^>/ {{$0=\">\" substr($0,2) \"_1\"}}1' > {output.swarm} && \
            vsearch --usearch_global {output.swarm}    --threads {threads}     --dbmask none     --qmask none     --rowlen 0     --notrunclabels     --userfields query+id1+target     --maxaccepts 0     --maxrejects 32  --top_hits_only  --output_no_hits     --db {params.db}     --id 0.5     --iddef 1     --userout {output.hits} && \
            python2.7 stampa_merge.py {params.stampadir} && \
            sed 's/:/_/' {output.results} | sed 's/|/;/g' | cut -f 1,4 > {output.taxonomy}
            """

    rule biom_tax_stampa:
        input:
            taxonomy="{project}/{prog}/stampa/{ds}.minsize{minsize}.{clmethod}.taxonomy.txt",
            biom="{project}/{prog}/otus/{ds}.minsize{minsize}.{clmethod}.biom",
        output:
            biom="{project}/{prog}/{ds}.minsize{minsize}.{clmethod}.taxonomy.biom",
            otutable="{project}/{prog}/{ds}.minsize{minsize}.{clmethod}.taxonomy.otutable.txt"
        conda: "envs/biom-format.yaml"
        shell: """biom add-metadata -i {input.biom} -o {output.biom} --observation-metadata-fp {input.taxonomy} --observation-header OTUID,taxonomy --sc-separated taxonomy --float-fields confidence --output-as-json && \
                  biom convert --to-tsv --header-key=taxonomy -i {output.biom} -o {output.otutable}
               """

rule report:
    input:
        readstat = "{project}/stats/readstat.csv",
        biom = expand("{{project}}/{prog}/{ds}.minsize{minsize}.{clmethod}.taxonomy.biom", prog=["vsearch"],ds=config['project'],minsize=2,clmethod=config['clustering']),
        otutable = expand("{{project}}/{prog}/{ds}.minsize{minsize}.{clmethod}.taxonomy.otutable.txt", prog=["vsearch"],ds=config['project'],minsize=2,clmethod=config['clustering'])        
    output:
        "{project}/stats/report.html"
    params:
        prefix="{project}/stats/report",
        mergemethod = config['mergepairs']
    conda: "envs/report.yaml"
    script:
        "report.Rmd"

