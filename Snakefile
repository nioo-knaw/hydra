from snakemake.utils import R
from snakemake.utils import min_version
import os

#min_version("3.12") # R Markdown reports have been introduced in Snakemake 3.12

if os.path.isfile("config.yaml"):
    configfile: "config.yaml"

# For a bash shell, needed on docker images to activate conda environments with source
# http://snakemake.readthedocs.io/en/stable/project_info/faq.html#i-want-to-configure-the-behavior-of-my-shell-for-all-rules-how-can-that-be-achieved-with-snakemake
shell.executable("/bin/bash")

PROJECT = config["project"] + "/"

rule final:
    input: expand("{project}/stats/contaminants.txt \
                   {project}/{prog}/clst/{ds}.minsize{minsize}.{clmethod}.fasta \
                   {project}/{prog}/{ds}.minsize{minsize}.{clmethod}.taxonomy.biom \
                   {project}/report/report.html".split(),data=config["data"],project=config['project'],prog=["vsearch"],ds=config['project'],minsize=config['minsize'],clmethod=config['clustering'])


from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider()

rule rename:
    input:
      forward = lambda wildcards: FTP.remote(config["data"][wildcards.data]["path"][0], keep_local=True, immediate_close=True) if config["remote"] else \
                                  config["data"][wildcards.data]["path"][0],
      reverse = lambda wildcards: FTP.remote(config["data"][wildcards.data]["path"][1], keep_local=True, immediate_close=True) if config["remote"] else \
                                  config["data"][wildcards.data]["path"][1]
    output:
        forward=protected("{project}/gunzip/{data}_R1.fastq.gz"),
        forward_md5=protected("{project}/gunzip/{data}_R1.fastq.gz.md5"),
        reverse=protected("{project}/gunzip/{data}_R2.fastq.gz"),
        reverse_md5=protected("{project}/gunzip/{data}_R2.fastq.gz.md5")

    params:
        prefix="{data}"
    threads: 8
    run: 
        if config["convert_to_casava1.8"]:
            # BUGFIX: For baseclear data, convert ti casava 1.8 format and add 0 as tag
            shell("zcat {input.forward} | awk '{{print (NR%4 == 1) ? \"@{params.prefix}_\" gsub(\"/1$\",\" 1:N:0:0\") substr($0,2) : $0}}' | gzip -c > {output.forward}")
            shell("md5sum {output.forward} > {output.forward_md5}")
            shell("zcat {input.reverse} | awk '{{print (NR%4 == 1) ? \"@{params.prefix}_\" gsub(\"/2$\",\" 2:N:0:0\") substr($0,2) : $0}}' | gzip -c > {output.reverse}")
            shell("md5sum {output.reverse} > {output.reverse_md5}")
        else:
            shell("zcat {input.forward} | awk '{{print (NR%4 == 1) ? \"@{params.prefix}_\" substr($0,2) : $0}}' | gzip -c > {output.forward}")
            shell("md5sum {output.forward} > {output.forward_md5}")
            shell("zcat {input.reverse} | awk '{{print (NR%4 == 1) ? \"@{params.prefix}_\" substr($0,2) : $0}}' | gzip -c > {output.reverse}")
            shell("md5sum {output.reverse} > {output.reverse_md5}")

rule filter_contaminants:
     input:
        forward="{project}/gunzip/{data}_R1.fastq.gz",
        reverse="{project}/gunzip/{data}_R2.fastq.gz"
     output:
        forward=temporary("{project}/filter/{data}_R1.fastq"),
        reverse=temporary("{project}/filter/{data}_R2.fastq"),
        stats="{project}/stats/{data}_contaminants_stats.txt"
     params:
         phix="refs/phix.fasta",
         adapters="refs/illumina_scriptseq_and_truseq_adapters.fa",
         quality=config["quality_control"]["trimming"]["quality"]
     log: "{project}/filter/{data}.log"
     conda: "envs/bbmap.yaml"
     threads: 16
     shell:"""bbduk.sh -Xmx8g in={input.forward} in2={input.reverse} out={output.forward} out2={output.reverse} \
              ref={params.adapters},{params.phix} qtrim="rl" trimq={params.quality} threads={threads} stats={output.stats} 2> {log}"""

rule contaminants_stats:
    input: expand("{project}/stats/{data}_contaminants_stats.txt",  project=config['project'], data=config["data"])
    output: 
        "{project}/stats/contaminants.txt"
    shell: "grep '#' -v {input} | tr ':' '\t' > {output}"

if config["barcode_in_header"]: 
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
            threshold=config['quality_control']['barcode']['threshold'],
            length=config['quality_control']['barcode']['length'],
            sep=config['quality_control']['barcode']['seperator'] 
        log: "{project}/barcode/{data}.log"
        conda: "envs/barcode.yaml"
        threads: 8
        shell: """extract_barcodes.py -f {input.forward}  -s{params.sep} -l {params.length} -o {params.outdir} -c barcode_in_label && fastq_to_fasta < {output.barcodes} > {output.barcodes_fasta} && \
                  trimmomatic PE -threads {threads} -phred33 {input.forward} {input.reverse} {output.forward} {output.forward_unpaired} {output.reverse} {output.reverse_unpaired} ILLUMINACLIP:{output.barcodes_fasta}:0:0:{params.threshold} 2> {log}"""

rule readstat_reverse:
    input:
          "{project}/barcode/{data}_R2.fastq" if config["barcode_in_header"] else\
          "{project}/filter/{data}_R2.fastq",
    output:
        temporary("{project}/stats/reverse/readstat.{data}.csv")
    log:
        "{project}/stats/reverse/readstat.{data}.log"
    conda:
        "envs/khmer.yaml"
    threads: 1
    shell: "readstats.py {input} --csv -o {output} 2> {log}"

rule readstat_reverse_merge:
    input:
        expand("{project}/stats/reverse/readstat.{data}.csv", project=config['project'], data=config["data"])
    output:
        protected("{project}/stats/readstat_R2.csv")
    shell: "cat {input[0]} | head -n 1 > {output} && for file in {input}; do tail -n +2 $file >> {output}; done;"


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

if config['mergepairs'] == 'none':
    rule copy_forward:
        input:
            forward="{project}/barcode/{data}_R1.fastq" if config["barcode_in_header"] else\
                    "{project}/filter/{data}_R1.fastq",
        output:
            fasta = temporary("{project}/mergepairs/{data}.fasta")
        conda: "envs/barcode.yaml"
        shell: "fastq_to_fasta < {input} > {output}"


if config['mergepairs'] == 'pandaseq':
    rule pandaseq:
        input:      
            forward="{project}/barcode/{data}_R1.fastq" if config["barcode_in_header"] else\
                    "{project}/filter/{data}_R1.fastq",
            reverse="{project}/barcode/{data}_R2.fastq" if config["barcode_in_header"] else\
                    "{project}/filter/{data}_R2.fastq",
        output:
            fasta = temporary("{project}/mergepairs/{data}.fasta")
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
        shell: "pandaseq -N -o {params.overlap} -l {params.minlength} -L {params.maxlength} -f {input.forward} -r {input.reverse} -T {threads} -w {output.fasta} -g {log}"

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
    rule vsearch_merge:
        input:
            forward="{project}/barcode/{data}_R1.fastq" if config["barcode_in_header"] else\
                    "{project}/filter/{data}_R1.fastq",
            reverse="{project}/barcode/{data}_R2.fastq" if config["barcode_in_header"] else\
                    "{project}/filter/{data}_R2.fastq",
        output:
            fasta = temporary("{project}/mergepairs/{data}.fasta")
        log: "{project}/mergepairs/{data}.log"
        threads: 1
        conda: "envs/vsearch.yaml"
        shell: "vsearch --threads {threads} --fastq_mergepairs {input.forward} --reverse {input.reverse} --fastq_allowmergestagger --fastq_minmergelen 200 --fastaout {output} > {log}"

if config['its'] == True:
    rule extract_its:
            input:
                    fasta="{project}/mergepairs/{data}.fasta"
            output:
                    fasta=temporary("{project}/itsx/{data}.ITS2.fasta")
            params:
                    basename="{project}/itsx/{data}",
                    dir="{project}/itsx"
            log: "itsx.log"
            threads: 32
            conda: "envs/itsx.yaml"
            # TODO: Filter on specific list of organisms? 
            # Only ITS2 region?
            shell: "ITSx --cpu {threads} --preserve TRUE -i {input.fasta} -o {params.basename} > {params.dir}/{log}"

# Combine per sample files to a single project file
rule mergefiles:
    input:
        fasta = expand(PROJECT + "itsx/{data}.ITS2.fasta", data=config["data"]) if config['its'] \
                else expand(PROJECT + "mergepairs/{data}.fasta", data=config["data"]),
    output: 
        fasta=temporary("{project}/mergefiles/{project}.fasta")
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
            amplicons="{project}/{prog}/{ds}.sorted.minsize{minsize}.fasta"
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
if config["classification"] == "sina":
    rule download_silva_arb:
        output: temporary("SSURef_NR99_128_SILVA_07_09_16_opt.arb")
        shell: """
RELEASE=128
URL="https://www.arb-silva.de/fileadmin/silva_databases/release_${{RELEASE}}/ARB_files"

FILE="SSURef_NR99_${{RELEASE}}_SILVA_07_09_16_opt.arb.gz"

# Download and check
wget -c ${{URL}}/${{FILE}}{{,.md5}} && md5sum -c ${{FILE}}.md5
gunzip ${{FILE}}
"""
    rule create_index_sina:
        input:
            "SSURef_NR99_128_SILVA_07_09_16_opt.arb"
        output:
            temporary("SSURef_NR99_128_SILVA_07_09_16_opt.arb.index.arb"),
            temporary("SSURef_NR99_128_SILVA_07_09_16_opt.arb.index.arb.pt"),
            temporary("SSURef_NR99_128_SILVA_07_09_16_opt.arb.index.ARM"),
            temporary("SSURef_NR99_128_SILVA_07_09_16_opt.arb.index.ARF")

        conda: "envs/sina.yaml"
        shell: """
cp SSURef_NR99_128_SILVA_07_09_16_opt.arb SSURef_NR99_128_SILVA_07_09_16_opt.arb.index.arb
arb_pt_server -build_clean -DSSURef_NR99_128_SILVA_07_09_16_opt.arb.index.arb
arb_pt_server -build -DSSURef_NR99_128_SILVA_07_09_16_opt.arb.index.arb
"""


    rule sina:
        input:
            fasta="{project}/{prog}/otus/{ds}.minsize{minsize}.{clmethod}.fasta",
            arb="SSURef_NR99_128_SILVA_07_09_16_opt.arb",
            index="SSURef_NR99_128_SILVA_07_09_16_opt.arb.index.arb",
            pt="SSURef_NR99_128_SILVA_07_09_16_opt.arb.index.arb.pt",
            arm="SSURef_NR99_128_SILVA_07_09_16_opt.arb.index.ARM",
            arf="SSURef_NR99_128_SILVA_07_09_16_opt.arb.index.ARF"
        output:
            #align="{project}/{prog}/sina/{ds}.minsize{minsize}.{clmethod}.sina.{chunk}.align",
            #aligncsv="{project}/{prog}/sina/{ds}.minsize{minsize}.{clmethod}.sina.{chunk}.align.csv",
            log="{project}/{prog}/sina/{ds}.minsize{minsize}.{clmethod}.sina.log",
            align="{project}/{prog}/sina/{ds}.minsize{minsize}.{clmethod}.sina.align",
        priority: -1
        threads: 24
        # TODO: turn is set to all to get classification. Reverse the reads in earlier stage!
        conda: "envs/sina.yaml"
        shell: "cat {input.fasta} | parallel --block 1000K -j{threads} --recstart '>' --pipe sina --log-file {output.log} -i /dev/stdin --intype fasta -o {output.align} --outtype fasta --meta-fmt csv --ptdb {input.arb} --overhang remove --turn all --search --search-db {input.arb} --search-min-sim 0.95 --search-no-fast --search-kmer-len 10 --lca-fields tax_slv || true"

    rule sina_get_taxonomy:
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

    rule sina_convert_tax:
        input:
            taxonomy="{project}/{prog}/sina/{ds}.minsize{minsize}.{clmethod}.sina.taxonomy.txt",
        output:
            taxonomy="{project}/{prog}/sina/{ds}.minsize{minsize}.{clmethod}.sina.qiimeformat.taxonomy",
        run:
            if config["use_full_lineage"]:
                shell("""awk -F"[;\t]" '{{printf $1"\t"; for(i=2;i<NF;i++){{printf i-1"__%s;", $i}}; printf "\\n"}}' {input.taxonomy} > {output.taxonomy}""")
            else:
                shell("""cat {input.taxonomy} | awk -F"[;\t]" 'BEGIN{{print "OTUs,Domain,Phylum,Class,Order,Family,Genus"}}{{print $1"\\tk__"$2"; p__"$3"; c__"$4"; o__"$5"; f__"$6"; g__"$7"; s__"$8}}' > {output.taxonomy}""")


    rule biom_tax_sina:
        input:
            taxonomy="{project}/{prog}/sina/{ds}.minsize{minsize}.{clmethod}.sina.qiimeformat.taxonomy",
            biom="{project}/{prog}/otus/{ds}.minsize{minsize}.{clmethod}.biom",
            meta=config["metadata"]
        output:
            biom=protected("{project}/{prog}/{ds}.minsize{minsize}.{clmethod}.taxonomy.biom"),
            otutable=protected("{project}/{prog}/{ds}.minsize{minsize}.{clmethod}.taxonomy.otutable.txt")
        conda: "envs/biom-format.yaml"
        shell: """
               biom add-metadata -i {input.biom} -o {output.biom} --observation-metadata-fp {input.taxonomy} --observation-header OTUID,taxonomy --sc-separated taxonomy --float-fields confidence --sample-metadata-fp {input.meta} --output-as-json && \
               biom convert --to-tsv --header-key=taxonomy -i {output.biom} -o {output.otutable}
               """

if config["classification"] == "stampa":
    rule download_silva:
        output:
            "SILVA_128_SSURef_tax_silva_trimmed.fasta"
        params:
            forward_primer=config["forward_primer"],
            reverse_primer=config["reverse_primer"]
        conda: "envs/cutadapt.yaml"
        # Download script adapted from https://github.com/frederic-mahe/stampa
        shell: """
RELEASE=128
URL="https://www.arb-silva.de/fileadmin/silva_databases/release_${{RELEASE}}/Exports"
FILE="SILVA_${{RELEASE}}_SSURef_tax_silva.fasta.gz"

# Download and check
wget -c ${{URL}}/${{FILE}}{{,.md5}} && md5sum -c ${{FILE}}.md5

# Define variables and output files
INPUT="SILVA_${{RELEASE}}_SSURef_tax_silva.fasta.gz"
OUTPUT="${{INPUT/.fasta.gz/_trimmed.fasta}}"
LOG="${{INPUT/.fasta.gz/_trimmed.log}}"
PRIMER_F={params.forward_primer}
PRIMER_R={params.reverse_primer}
MIN_LENGTH=32
MIN_F=$(( ${{#PRIMER_F}} * 2 / 3 ))
MIN_R=$(( ${{#PRIMER_R}} * 2 / 3 ))
CUTADAPT="cutadapt --discard-untrimmed --minimum-length ${{MIN_LENGTH}}"

# Trim forward & reverse primers, format
zcat "${{INPUT}}" | sed '/^>/ ! s/U/T/g' | \
     ${{CUTADAPT}} -g "${{PRIMER_F}}" -O "${{MIN_F}}" - 2> "${{LOG}}" | \
     ${{CUTADAPT}} -a "${{PRIMER_R}}" -O "${{MIN_F}}" - 2>> "${{LOG}}" | \
     sed '/^>/ s/;/|/g ; /^>/ s/ /_/g ; /^>/ s/_/ /1' > "${{OUTPUT}}"
    """

    rule download_unite:
        output:
            "UNITE-7.2-28.06.2017.fasta"
        params:
            forward_primer=config["forward_primer"],
            reverse_primer=config["reverse_primer"]
        conda: "envs/cutadapt.yaml"

        shell:"""
RELEASE=7.2
DATE=28.06.2017
URL="https://unite.ut.ee/sh_files/"
FILE="UNITE_public_${{DATE}}.fasta.zip"

# Download and check
wget -c ${{URL}}/${{FILE}}
unzip ${{FILE}}

# Define variables and output files
PRIMER_F={params.forward_primer}
PRIMER_R={params.reverse_primer}
FNAME=ITS9
RNAME=ITS4
INPUT="UNITE_public_${{DATE}}.fasta"
OUTPUT="UNITE-${{RELEASE}}-${{DATE}}.fasta"
LOG="UNITE-${{RELEASE}}_${{DATE}}.log"
MIN_LENGTH=32
MIN_F=$(( ${{#PRIMER_F}} * 2 / 3 ))
MIN_R=$(( ${{#PRIMER_R}} * 2 / 3 ))
CUTADAPT="cutadapt --discard-untrimmed --minimum-length ${{MIN_LENGTH}}"

# Trim forward & reverse primers, format
cat "${{INPUT}}" | sed '/^>/ ! s/U/T/g' | \
     ${{CUTADAPT}} -g "${{PRIMER_F}}" -O "${{MIN_F}}" - 2> "${{LOG}}" | \
     ${{CUTADAPT}} -a "${{PRIMER_R}}" -O "${{MIN_F}}" - 2>> "${{LOG}}" | \
     sed 's/\ /+/'g | sed 's/|/\ /'g | sed 's/;/|/g' | \
     awk '/^>/ {{env=index($0,"s__Fungi_sp");}} {{if (!env) print;}}' | \
     awk '/^>/ {{env=index($0,"k__Fungi");}} {{if (env) print;}}' > "${{OUTPUT}}"
"""

    rule stampa:
        input:
            fasta="{project}/{prog}/otus/{ds}.minsize{minsize}.{clmethod}.fasta",
            db = "SILVA_128_SSURef_tax_silva_trimmed.fasta" if config["reference_db"] == "silva" else "UNITE_public_28.06.2017.fasta" if config["reference_db"] == "unite" else None,
        output:
            swarm="{project}/{prog}/stampa/{ds}.minsize{minsize}.{clmethod}.fasta",
            hits="{project}/{prog}/stampa/hits.{ds}.minsize{minsize}.{clmethod}_usearch_global",
            results="{project}/{prog}/stampa/results.{ds}.minsize{minsize}.{clmethod}_usearch_global",
            taxonomy="{project}/{prog}/stampa/{ds}.minsize{minsize}.{clmethod}.taxonomy.txt",
        params:
             stampadir="{project}/{prog}/stampa/",
        conda: "envs/vsearch.yaml"
        threads: 32
        # Create STAMPA compatible input
        # Replace underscore in otu names and add fake abundance information
        shell:"""
            sed 's/_/:/' {input.fasta} | awk '/^>/ {{$0=\">\" substr($0,2) \"_1\"}}1' > {output.swarm} && \
            vsearch --usearch_global {output.swarm}    --threads {threads}     --dbmask none     --qmask none     --rowlen 0     --notrunclabels     --userfields query+id1+target     --maxaccepts 0     --maxrejects 32  --top_hits_only  --output_no_hits     --db {input.db}     --id 0.5     --iddef 1     --userout {output.hits} && \
            python2.7 stampa_merge.py {params.stampadir} && \
            sed 's/:/_/' {output.results} | sed 's/|/;/g' | cut -f 1,4 > {output.taxonomy}
            """

    rule biom_tax_stampa:
        input:
            taxonomy="{project}/{prog}/stampa/{ds}.minsize{minsize}.{clmethod}.taxonomy.txt",
            biom="{project}/{prog}/otus/{ds}.minsize{minsize}.{clmethod}.biom",
            meta=config["metadata"]
        output:
            taxonomy="{project}/{prog}/stampa/{ds}.minsize{minsize}.{clmethod}.taxonomy.qiimeformat.txt",
            biom="{project}/{prog}/{ds}.minsize{minsize}.{clmethod}.taxonomy.biom",
            otutable="{project}/{prog}/{ds}.minsize{minsize}.{clmethod}.taxonomy.otutable.txt"
        conda: "envs/biom-format.yaml"
        shell: """cat {input.taxonomy} | awk -F"[;\t]" 'BEGIN{{print "OTUs,Domain,Phylum,Class,Order,Family,Genus"}}{{print $1"\\tk__"$2"; p__"$3"; c__"$4"; o__"$5"; f__"$6"; g__"$7"; s__"$8}}' > {output.taxonomy} && \
                  biom add-metadata -i {input.biom} -o {output.biom} --observation-metadata-fp {output.taxonomy} --observation-header OTUID,taxonomy --sc-separated taxonomy --float-fields confidence --sample-metadata-fp {input.meta} --output-as-json && \
                  biom convert --to-tsv --header-key=taxonomy -i {output.biom} -o {output.otutable}
               """

if config["classification"] == "blast":
    rule blastn:
        input:
            "{project}/{prog}/otus/{ds}.minsize{minsize}.{clmethod}.fasta"
        output:
            hits="{project}/{prog}/blast/{ds}.minsize{minsize}.{clmethod}.blastout.txt",
        params:
            db = config["blast_db"],
            max_hits = config["blast_max_hits"]
        threads: 32
        conda: "envs/blast.yaml"
        shell: """blastn -query {input} -db {params.db} -evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"  -out {output.hits} -num_threads {threads} -max_target_seqs {params.max_hits}"""

    rule lca:
        input:
            hits="{project}/{prog}/blast/{ds}.minsize{minsize}.{clmethod}.blastout.txt",
            fasta = "{project}/{prog}/otus/{ds}.minsize{minsize}.{clmethod}.fasta"
        output:
            lca="{project}/{prog}/blast/{ds}.minsize{minsize}.{clmethod}.lca.txt",
            taxonomy="{project}/{prog}/blast/{ds}.minsize{minsize}.{clmethod}.taxonomy.txt"
        params:
            taxref = "/data/db/pr2/gb203/lotus_lca.tax",
        run:
            shell("./LCA -i {input.hits} -r {params.taxref} -o {output.lca} -id '99,96,93,91,88,78,0'")
            shell("""cat {output.lca} | awk -F"\t" 'BEGIN{{}}{{gsub(" ","_",$0);gsub("\\"","",$0);print $1"\\td__"$2";p__"$3";c__"$4";o__"$5";f_"$6";g__"$7";s__"$8}}' > {output.taxonomy}""")
            # Detect OTUs without blast hit
            hits = []
            with open(output.lca) as fh:
                for line in fh:
                    line = line.strip().split("\t")
                    hits.append(line[0])
            with open(input.fasta) as fasta_file, open(output.taxonomy, "a") as outfile:
                for line in fasta_file:
                    line = line.strip()
                    if line.startswith(">"):
                        if line[1:] not in hits:
                            outfile.write("%s\tk__?;p__?;c__?;o__?;f__?;g__?;s__?\n" % line[1:])

    rule biom_tax_blast:
        input:
            taxonomy="{project}/{prog}/blast/{ds}.minsize{minsize}.{clmethod}.taxonomy.txt",
            biom="{project}/{prog}/otus/{ds}.minsize{minsize}.{clmethod}.biom",
            meta=config["metadata"]
        output:
            biom="{project}/{prog}/{ds}.minsize{minsize}.{clmethod}.taxonomy.biom",
            otutable="{project}/{prog}/{ds}.minsize{minsize}.{clmethod}.taxonomy.otutable.txt"
        conda: "envs/biom-format.yaml"
        shell: """biom add-metadata -i {input.biom} -o {output.biom} --observation-metadata-fp {input.taxonomy} --observation-header OTUID,taxonomy --sc-separated taxonomy --float-fields confidence --sample-metadata-fp {input.meta} --output-as-json && \
                  biom convert --to-tsv --header-key=taxonomy -i {output.biom} -o {output.otutable}
               """

if config["classification"] == "rdp":
    rule rdp:
        input:
            fasta = "{project}/{prog}/otus/{ds}.minsize{minsize}.{clmethod}.fasta"
        output:
             "{project}/{prog}/rdp/{ds}.minsize{minsize}.{clmethod}.rdp"
        params:
            traindir="/data/db/unite/UNITE_retrained/" 
        conda: "envs/rdp.yaml"
        shell: "classifier classify -Xms512M -Xmx8g -t {params.traindir}/rRNAClassifier.properties -o {output} {input}"

    rule rdp_filter:
        input:
            rdpout = "{project}/{prog}/rdp/{ds}.minsize{minsize}.{clmethod}.rdp",
            fasta = "{project}/{prog}/otus/{ds}.minsize{minsize}.{clmethod}.fasta"
        output:
             "{project}/{prog}/rdp/{ds}.minsize{minsize}.{clmethod}.filtered.rdp"
        params:
            cutoff = config["rdp_confidence_cutoff"]
        run:
            hits = []
            with open(input.rdpout) as rdpout, open(output[0], "w") as outfile:
                for line in rdpout:
                    parts = line.strip().split("\t")
                    tax_list = []
                   
                    for i in range(5, len(parts), 3):
                         confidence = float(parts[i+2])
                         taxonomy = parts[i].split("|")[-1].replace(" ", "_")
                         tax_level = parts[i+1]

                         if tax_level == "domain" or tax_level == "kingdom":
                             tax_prefix = "k__" 
                         elif tax_level == "phylum":
                             tax_prefix = "p__"
                         elif tax_level == "class":
                             tax_prefix = "c__"
                         elif tax_level == "order":
                             tax_prefix = "o__"
                         elif tax_level == "family":
                             tax_prefix = "f__"
                         elif tax_level == "genus":
                             tax_prefix = "g__"
                         elif tax_level == "species":
                             tax_prefix = "s__"

                         if confidence > float(params.cutoff):
                             tax_list.append(tax_prefix + taxonomy)
                         else:
                             tax_list.append(tax_prefix + "unclassified")

                    tax_string = ";".join(tax_list)
                    outfile.write("%s\t%s\n" % (parts[0],tax_string))
                    # Keep track of OTUs with classification
                    hits.append(parts[0])
            
            # Add OTUs without RDP classification
            with open(input.fasta) as fasta_file, open(output[0], "a") as outfile:
                for line in fasta_file:
                    line = line.strip()
                    if line.startswith(">"):
                        if line[1:] not in hits:
                            outfile.write("%s\td__unclassified;p__unclassified;c__unclassified;o__unclassified;f__unclassified;g__unclassified;s__unclassified\n" % line[1:])

             

    rule biom_tax_rdp:
        input:
            biom="{project}/{prog}/otus/{ds}.minsize{minsize}.{clmethod}.biom",
            taxonomy="{project}/{prog}/rdp/{ds}.minsize{minsize}.{clmethod}.filtered.rdp",
            meta=config["metadata"]
        output:
            biom="{project}/{prog}/{ds}.minsize{minsize}.{clmethod}.taxonomy.biom",
            otutable="{project}/{prog}/{ds}.minsize{minsize}.{clmethod}.taxonomy.otutable.txt"
        conda: "envs/biom-format.yaml"
        shell: """biom add-metadata -i {input.biom} -o {output.biom} --observation-metadata-fp {input.taxonomy} --observation-header OTUID,taxonomy --sc-separated taxonomy --float-fields confidence --sample-metadata-fp {input.meta} --output-as-json && \
                  biom convert --to-tsv --header-key=taxonomy -i {output.biom} -o {output.otutable}
               """

rule workflow_graph:
    output: temporary("{project}/report/workflow.svg")
    conda: "envs/rulegraph.yaml"
    shell: "snakemake --rulegraph | dot -Tsvg > {output}"


rule report:
    input:
        workflow =  "{project}/report/workflow.svg",
        readstat = "{project}/stats/readstat.csv",
        readstat_reverse = "{project}/stats/readstat_R2.csv",
        biom = expand("{{project}}/{prog}/{ds}.minsize{minsize}.{clmethod}.taxonomy.biom", prog=["vsearch"],ds=config['project'],minsize=2,clmethod=config['clustering']),
        otutable = expand("{{project}}/{prog}/{ds}.minsize{minsize}.{clmethod}.taxonomy.otutable.txt", prog=["vsearch"],ds=config['project'],minsize=2,clmethod=config['clustering']),
        otus= expand("{{project}}/{prog}/otus/{ds}.minsize{minsize}.{clmethod}.fasta", prog=["vsearch"],ds=config['project'],minsize=2,clmethod=config['clustering']),
    output:
        "{project}/report/report.html"
    params:
        prefix="{project}/report/report",
        mergemethod = config['mergepairs']
    conda: "envs/report.yaml"
    script:
        "report.Rmd"

