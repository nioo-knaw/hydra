derep_swarm:
    input:
        "{project}/trim/{data}.fasta"
    output:
        "{project}/swarm/{data}.fasta"
    params:
        minsize="1"
    shell: """
grep -v "^>" {input} | \
grep -v [^ACGTacgt] | sort -d | uniq -c | \
while read abundance sequence ; do
    if [[ ${{abundance}} -ge {params.minsize} ]]; then
        hash=$(printf "${{sequence}}" | sha1sum)
        hash=${{hash:0:40}}
        printf ">%s_%d_%s\\n" "${{hash}}" "${{abundance}}" "${{sequence}}"
    fi
done | sort -t "_" -k2,2nr -k1.2,1d | \
sed -e 's/\_/\\n/2' > {output}
"""

rule swarm_amplicon_contigency:
    input:
        fasta = expand(PROJECT + "swarm/{data}.fasta", data=config["data"])
    output:
        count="{project}/swarm/{project}_amplicon_contingency_table.csv",
        fasta="{project}/swarm/{project}.fasta"
    run:
        shell("python2.7 /data/tools/swarm/1.2.20/amplicon_contingency_table.py {input} > {output.count}")
        shell("cat {input} > {output.fasta}")

rule swarm:
    input: 
        "{project}/swarm/{project}.fasta"
    output:
        swarms="{project}/swarm/{project}.swarm_d{d}.swarms",
        stats="{project}/swarm/{project}.swarm_d{d}.stats"
    params: d="{d}" 
    threads: 16
    # The alternative algorithm only works with d=1!
    shell: "/data/tools/swarm/1.2.20/swarm -d {params.d} -t {threads} -s {output.stats} < {input}  > {output.swarms}"

rule swarm_get_seed:
    input: 
        swarms="{project}/swarm/{project}.swarm_d{d}.swarms",
        amplicons="{project}/swarm/{project}.fasta"
    output:
        seeds="{project}/swarm/{project}.swarm_d{d}.seeds.fasta"
    shell: "SEEDS=$(mktemp); cut -d ' ' -f 1 {input.swarms} | sed -e 's/^/>/' > '${{SEEDS}}'; grep -A 1 -F -f '${{SEEDS}}' {input.amplicons} | sed -e '/^--$/d' > {output.seeds}"

rule swarm_otu_contigency:
     input: 
        stats="{project}/swarm/{project}.swarm_d{d}.stats",
        swarms="{project}/swarm/{project}.swarm_d{d}.swarms",
        amplicon_count="{project}/swarm/{project}_amplicon_contingency_table.csv"
     output:
        otu_count="{project}/swarm/{project}.swarm_d{d}.otu_contingency_table.csv"
     priority: 50
     shell:"""
echo -e "OTU\\t$(head -n 1 "{input.amplicon_count}")" > "{output.otu_count}"

awk -v SWARM="{input.swarms}"     -v TABLE="{input.amplicon_count}"     'BEGIN {{FS = " "
            while ((getline < SWARM) > 0) {{
                swarms[$1] = $0
            }}
            FS = "\\t"
            while ((getline < TABLE) > 0) {{
                table[$1] = $0
            }}
           }}

     {{# Parse the stat file (OTUs sorted by decreasing abundance)
      seed = $3 "_" $4
      n = split(swarms[seed], OTU, "[ _]")
      for (i = 1; i < n; i = i + 2) {{
          s = split(table[OTU[i]], abundances, "\\t")
          for (j = 1; j < s; j++) {{
              samples[j] += abundances[j+1]
          }}
      }}
      printf "%s\\t%s", NR, $3
      for (j = 1; j < s; j++) {{
          printf "\\t%s", samples[j]
      }}
     printf "\\n"
     delete samples
     }}' "{input.stats}" >> "{output.otu_count}"
    """

rule swarm_biom:
    input:
        otu_count="{project}/swarm/{project}.swarm_d{d}.otu_contingency_table.csv"
    output:
        otutable="{project}/swarm/{project}.swarm_d{d}.otutable.txt",
        biom="{project}/swarm/{project}.swarm_d{d}.biom"
    run:
        #shell("cut -f 1,3,4 {input} > {output.otutable}")
        shell("""awk  '{{delim = ""; for (i=1;i<=NF-1;i++) if (i!=2) {{printf delim "%s", $i; delim = "\\t"}}; printf "\\n"}}' < {input} > {output.otutable}""")
        shell("/data/tools/qiime/1.9/qiime1.9/bin/biom convert -i {output.otutable} -o {output.biom} --table-type='otu table'")

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
        otus=protected("{project}/{prog}/{ds}.minsize{minsize}.usearch_smallmem.fasta")
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

# Swarm
rule swarm_vsearch:
    input: 
        "{project}/{prog}/{ds}.sorted.minsize{minsize}.fasta"
    output:
        swarms="{project}/{prog}/{ds}.minsize{minsize}.swarm_d{d}.swarms",
        stats="{project}/{prog}/{ds}.minsize{minsize}.swarm_d{d}.stats"
    params: d="{d}" 
    threads: 16
    # The alternative algorithm only works with d=1!
    shell: "/data/tools/swarm/1.2.20/swarm -d {params.d} -t {threads} -z -u uclust.out -s {output.stats} < {input}  > {output.swarms}"

rule swarm_get_seed_vsearch:
    input: 
        swarms="{project}/{prog}/{ds}.minsize{minsize}.swarm_d{d}.swarms",
        amplicons="{project}/{prog}/{ds}.sorted.minsize{minsize}.fasta"
    output:
        seeds="{project}/{prog}/{ds}.minsize{minsize}.swarm_d{d,\d+}.fasta"
    shell: "SEEDS=$(mktemp); cut -d ' ' -f 1 {input.swarms} | sed -e 's/^/>/' > '${{SEEDS}}'; grep -A 1 -F -f '${{SEEDS}}' {input.amplicons} | sed -e '/^--$/d' > {output.seeds}"


# HPC-Clust
rule hpc_clust:
    input: "{project}/filterseqs/{project}.unique.good.filter.fasta"
    output: 
        log="{project}/hpc-clust/{project}.{cl}",
        otus="{project}/hpc-clust/{project}.{cl}.otus"
    params:
        prefix="{project}/hpc-clust/{project}",
        method="{cl}"
    log: "{project}/hpc-clust/{project}.log"
    threads: 16
    run: 
        shell("source /data/tools/hpc-clust/1.2.0/env.sh; hpc-clust -nthreads {threads} -{params.method} true -t 0.95 {input} -ofile {params.prefix} -ignoreMemThres true > {log}")
        shell("/data/tools/hpc-clust/1.2.0/bin/make-otus-mothur.sh {input} {output.log} 0.97 > {output.otus}")

#
# Chimera checking
#

rule uchime:
    input:
        "{project}/{prog}/{ds}.minsize{minsize}.{clmethod}.fasta"
    output:
        chimeras="{project}/{prog}/{ds}.minsize{minsize}.{clmethod}.chimeras.fasta",
        nonchimeras="{project}/{prog}/{ds}.minsize{minsize}.{clmethod}.nonchimeras.fasta"
    log: "{project}/{prog}/{ds}.minsize{minsize}.{clmethod}.uchime.log"
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
        "{project}/{prog}/{ds}.minsize{minsize}.{clmethod}.nonchimeras.fasta"
    output:
        "{project}/{prog}/{ds}.minsize{minsize}.{clmethod}.otus.fasta"
    shell: "python2.7 /data/tools/usearch/uparse_scripts/fasta_number.py {input} OTU_ > {output}"

rule mapping:
    input:
        otus="{project}/{prog}/{ds}.minsize{minsize}.{clmethod}.otus.fasta",
        reads="{project}/mergefiles/{ds}.fasta"
    output:
        "{project}/{prog}/{ds}.minsize{minsize}.{clmethod}.otus.uc"
    run:
        cmd = ""
        if wildcards.prog == "vsearch":
            cmd = VSEARCH
        elif wildcards.prog == "usearch":
            cmd = USEARCH
        shell("{cmd} -usearch_global {input.reads} -db {input.otus} -strand plus -id 0.97 -uc {output}")

#http://rpackages.ianhowson.com/bioc/phyloseq/man/import_usearch_uc.html
rule create_otutable_phyloseq:
    input:
        "{project}/{prog}/{ds}.minsize{minsize}.{clmethod}.otus.uc"
    output:
        "{project}/{prog}/{ds}.minsize{minsize}.{clmethod}.otus.phyloseq.otutable"
    run:
        R("""
          require('phyloseq')
          OTU <- import_usearch_uc('{input}', colRead = 9, colOTU = 10,readDelimiter = "_", verbose = TRUE)
          write.table(t(OTU[,2:ncol(OTU)]),file="{output}",sep=\"\t\")
          """)

rule create_otutable:
    input:
        "{project}/{prog}/{ds}.minsize{minsize}.{clmethod}.otus.uc"
    output:
        "{project}/{prog}/{ds}.minsize{minsize}.{clmethod}.otus.otutable.txt"
    shell: "python2.7 /data/tools/usearch/uparse_scripts/uc2otutab.py {input} > {output}"

# convert to biom file
rule biom_otu:
    input: 
        "{project}/{prog}/{ds}.minsize{minsize}.{clmethod}.otus.otutable.txt"
    output:
        "{project}/{prog}/{ds}.minsize{minsize}.{clmethod}.biom" 
    shell: "/data/tools/qiime/1.9/qiime1.9/bin/biom convert -i {input} --to-json -o {output} --table-type='OTU table'"

#
# Taxonomy
#
SINA_VERSION=config['sina_version']
SILVADB = config['silva_db']
SINA_MIN_SIM = config['sina_min_sim']

rule sina_parallel:
    input:
        "{project}/swarm/{project}.swarm_d{d}.seeds.fasta"
    output:
        align="{project}/sina/{project}.swarm_d{d}.sina.align",
        log="{project}/sina/{project}.swarm_d{d}.sina.log"
    log: "{project}/sina/{project}.swarm_d{d}.sina.log"
    priority: -1
    threads: 8
    # TODO: turn is set to all to get classification. Reverse the reads in earlier stage!
    shell: "cat {input} | parallel --block 1000K -j{threads} --recstart '>' --pipe /data/tools/sina/{SINA_VERSION}/sina --log-file {log} -i /dev/stdin -o {output.align} --outtype fasta --meta-fmt csv --ptdb /scratch/silva/SSURef_NR99_119_SILVA_14_07_14_opt.arb --overhang remove --turn all --search --search-db /scratch/silva/SSURef_NR99_119_SILVA_14_07_14_opt.arb --search-min-sim 0.95 --search-no-fast --search-kmer-len 10 --lca-fields tax_slv"

rule sina_get_taxonomy_from_logfile:
    input: log="{project}/sina/{project}.swarm_d{d}.sina.log"
    output: taxonomy="{project}/sina/{project}.swarm_d{d}.sina.taxonomy"
    # Parse the log file from Sina to get the taxonomic classification
    # The csv output does not contain the sequence identifier, thats way this approach is better
    # The first space needs to be replaced in order to keep the space in the taxonomy string (would be splitted otherwise)
    # Brackets are escaped by an extra bracket, because they are internaly recognised by Snakemake
    shell: "cat {input.log} | sed 's/ /|/1' | awk -F '|'  '/^sequence_identifier:/ {{id=$2}} /^lca_tax_slv:/{{split(id,a,\" \"); print a[1] \"\t\" $2}}' | tr ' ' '_' > {output.taxonomy}"

rule rdp:
    input:
        "{project}/swarm/{project}.swarm_d{d}.seeds.fasta"
    output:
        otus="{project}/rdp/{project}.swarm_d{d}.seeds.otus.fasta",
        rdp="{project}/rdp/{project}.swarm_d{d}.seeds.rdp",
        taxonomy="{project}/rdp/{project}.swarm_d{d}.seeds.rdp.taxonomy.txt"
    run:
        shell("python2.7 /data/tools/usearch/uparse_scripts/fasta_number.py {input} > {output.otus}") 
        shell("java -Xmx1g -jar /data/tools/rdp-classifier/2.10/classifier.jar classify -c 0.8 {output.otus} -f filterbyconf -o {output.rdp}")
        shell("""cat {output.rdp} | awk -F"\\t" 'BEGIN{{print "OTUs,Domain,Phylum,Class,Order,Family,Genus"}}{{gsub(" ","_",$0);gsub("\\"","",$0);print $1"\\t"$2";"$3";"$4";"$5";"$6";"$7}}' > {output.taxonomy}""")

rule rdp_edgar:
    input:
        "{project}/{prog}/{ds}.minsize{minsize}.{clmethod}.otus.fasta"
    output:
        rdp="{project}/{prog}/rdp/{ds}.minsize{minsize}.{clmethod}.otus.rdp",
        taxonomy="{project}/{prog}/rdp/{ds}.minsize{minsize}.{clmethod}.otus.rdp.taxonomy"
    run: 
        shell("java -Xmx1g -jar /data/tools/rdp-classifier/2.10/classifier.jar classify -c 0.8 {input} -f filterbyconf -o {output.rdp}")
        shell("""cat {output.rdp} | awk -F"\\t" 'BEGIN{{print "OTUs,Domain,Phylum,Class,Order,Family,Genus"}}{{gsub(" ","_",$0);gsub("\\"","",$0);print $1"\\t"$2";"$3";"$4";"$5";"$6";"$7}}' > {output.taxonomy}""")

rule sina_parallel_edgar:
    input:
        "{project}/{prog}/{ds}.minsize{minsize}.{clmethod}.otus.fasta"
    output:
        align="{project}/{prog}/sina/{ds}.minsize{minsize}.{clmethod}.sina.align",
        log="{project}/{prog}/sina/{ds}.minsize{minsize}.{clmethod}.sina.log"
    log: "{project}/{prog}/sina/{ds}.minsize{minsize}.{clmethod}.sina.log"
    priority: -1
    threads: 8
    # TODO: turn is set to all to get classification. Reverse the reads in earlier stage!
    shell: "cat {input} | parallel --block 100K -j{threads} --recstart '>' --pipe /data/tools/sina/{SINA_VERSION}/sina --log-file {log} -i /dev/stdin -o {output.align} --outtype fasta --meta-fmt csv --ptdb /scratch/silva/SSURef_NR99_119_SILVA_14_07_14_opt.arb --overhang remove --turn all --search --search-db /scratch/silva/SSURef_NR99_119_SILVA_14_07_14_opt.arb --search-min-sim 0.95 --search-no-fast --search-kmer-len 10 --lca-fields tax_slv"

rule sina_get_taxonomy_from_logfile_edgar:
    input: log="{project}/{prog}/sina/{ds}.minsize{minsize}.{clmethod}.sina.log"
    output: taxonomy="{project}/{prog}/sina/{ds}.minsize{minsize}.{clmethod}.sina.taxonomy"
    # Parse the log file from Sina to get the taxonomic classification
    # The csv output does not contain the sequence identifier, thats way this approach is better
    # The first space needs to be replaced in order to keep the space in the taxonomy string (would be splitted otherwise)
    # Brackets are escaped by an extra bracket, because they are internaly recognised by Snakemake
    shell: "cat {input.log} | sed 's/ /|/1' | awk -F '|'  '/^sequence_identifier:/ {{id=$2}} /^lca_tax_slv:/{{split(id,a,\" \"); print a[1] \"\t\" $2}}' | tr ' ' '_' > {output.taxonomy}"

rule biom_tax_rdp:
    input:
        taxonomy="{project}/{prog}/rdp/{ds}.minsize{minsize}.{clmethod}.otus.rdp.taxonomy",
        biom="{project}/{prog}/{ds}.minsize{minsize}.{clmethod}.biom"
    output:
        biom=protected("{project}/{prog}/{ds}.minsize{minsize}.{clmethod}.taxonomy.biom"),
        otutable=protected("{project}/{prog}/{ds}.minsize{minsize}.{clmethod}.taxonomy.otutable.txt")
    run:
        shell("/data/tools/qiime/1.9/qiime1.9/bin/biom add-metadata -i {input.biom} -o {output.biom} --observation-metadata-fp {input.taxonomy} --observation-header OTUID,taxonomy --sc-separated taxonomy --float-fields confidence")
        shell("/data/tools/qiime/1.9/qiime1.9/bin/biom convert --to-tsv --header-key=taxonomy -i {output.biom} -o {output.otutable}")

rule biom_tax_sina:
    input:
        taxonomy="{project}/{prog}/sina/{ds}.minsize{minsize}.{clmethod}.sina.taxonomy",
        biom="{project}/{prog}/{ds}.minsize{minsize}.{clmethod}.biom"
    output:
        biom=protected("{project}/{prog}/sina/{ds}.minsize{minsize}.{clmethod}.taxonomy.biom"),
        otutable=protected("{project}/{prog}/sina/{ds}.minsize{minsize}.{clmethod}.taxonomy.otutable.txt")
    run:
        shell("/data/tools/qiime/1.9/qiime1.9/bin/biom add-metadata -i {input.biom} -o {output.biom} --observation-metadata-fp {input.taxonomy} --observation-header OTUID,taxonomy --sc-separated taxonomy --float-fields confidence")
        shell("/data/tools/qiime/1.9/qiime1.9/bin/biom convert --to-tsv --header-key=taxonomy -i {output.biom} -o {output.otutable}")


rule biom_tax_swarm:
    input:
        taxonomy="{project}/rdp/{project}.swarm_d{d}.seeds.rdp.taxonomy.txt",
        biom="{project}/swarm/{project}.swarm_d{d}.biom"
    output:
        biom="{project}/rdp/{project}.swarm_d{d}.taxonomy.biom",
        otutable="{project}/rdp/{project}.swarm_d{d}.taxonomy.otutable.txt"
    run:
        shell("/data/tools/qiime/1.9/qiime1.9/bin/biom add-metadata -i {input.biom} -o {output.biom} --observation-metadata-fp {input.taxonomy} --observation-header OTUID,taxonomy --sc-separated taxonomy --float-fields confidence")
        shell("/data/tools/qiime/1.9/qiime1.9/bin/biom convert -b --header-key=taxonomy -i {output.biom} -o {output.otutable}")


           
# Filter OTUs
#filter_otus_from_otu_table.py  -i swarm.biom -o filter.biom -n 5 
