import logging
import multiprocessing
import re
import os
import tempfile
import yaml
import sys
from collections import OrderedDict
import click
import urllib

# Adapted from: https://github.com/pnnl/atlas/blob/master/atlas/conf.py

logging.basicConfig(level=logging.INFO, datefmt="%Y-%m-%d %H:%M", format="[%(asctime)s %(levelname)s] %(message)s")

host = "ftp.sra.ebi.ac.uk"
project = "PRJEB14409"
#project = "PRJNA319605"

# http://stackoverflow.com/a/3675423
def replace_last(source_string, replace_what, replace_with):
    head, _sep, tail =  source_string.rpartition(replace_what)
    if _sep == '':
        return tail
    else: 
        return head + replace_with + tail

def get_ena(project):
    from urllib import request
    samples = ""
    try:
        samples = request.urlopen("http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=%s&result=read_run&fields=fastq_ftp" % project).readlines()[1:]
    except urllib.error.HTTPError:
        print("Not a valid ENA project")
    for sample in samples:
        for fastq in sample.strip().split(b';'):
            dirpath = os.path.dirname(fastq).decode("utf-8")
            filename = os.path.basename(fastq).decode("utf-8")
            yield (dirpath,"",[filename])

def get_sample_files(path, remote):
    samples = OrderedDict()
    seen = set()
    walker = ""

    if remote != None:
        walker = get_ena(remote)
    else:
        walker = os.walk(path, followlinks=True)
    for dir_name, sub_dirs, files in walker:
        for fname in files:
            if ".fastq" in fname or ".fq" in fname:
                sample_id = fname.partition(".fastq")[0]

                if ".fq" in sample_id:
                    sample_id = fname.partition(".fq")[0].replace("_","-")

                sample_id = sample_id.replace("_R1", "").replace("_r1", "").replace("_R2", "").replace("_r2", "")
                sample_id = re.sub("_1$", "", sample_id)
                sample_id = re.sub("_2$", "", sample_id)
                sample_id = sample_id.replace("_", "-").replace(" ", "-")

                fq_path = os.path.join(dir_name, fname)
                fastq_paths = [fq_path]
                if fq_path in seen: continue

                if "_R1" in fname or "_r1" in fname or "_1" in fname:
                    fname = replace_last(fname,"_1.","_2.")
                    r2_path = os.path.join(dir_name, fname.replace("_R1", "_R2").replace("_r1", "_r2"))
                    if not r2_path == fq_path:
                        seen.add(r2_path)
                        fastq_paths.append(r2_path)

                if "_R2" in fname or "_r2" in fname or "_2" in fname:
                    fname = replace_last(fname,"_2.","_1.")
                    r1_path = os.path.join(dir_name, fname.replace("_R2", "_R1").replace("_r2", "_r1"))
                    if not r1_path == fq_path:
                        seen.add(r1_path)
                        fastq_paths.insert(0, r1_path)

                if sample_id in samples:
                    logging.warn("Duplicate sample %s was found after renaming; skipping..." % sample_id)
                    continue

                samples[sample_id] = {'path': fastq_paths }
    return samples

def create_metadata_template(outfile, samples):
   with open(outfile, "w") as f:
       print("#SampleID\tAlias", file=f)
       for sample in samples:
           print("%s\t%s" %  (sample,sample), file=f)

@click.command()
@click.option('--project', prompt="Give your project a unique name", required=True, help='Give your project a nice name')
@click.option('--config', default="config.yaml", show_default=True, help='File to write the configuration to')
@click.option('--remote', help='Specify a ENA project to use as remote data (for example PRJEB14409')
@click.option('--path', default="../data", show_default=True, help='path to data folder')
@click.option('--rename', required=False, help='provide a file for renaming samples')

@click.option('--forward_primer', prompt="Which forward primer did you use?", required=True, default="CCTACGGGNGGCWGCAG", help="Which forward primer did you use?")
@click.option('--reverse_primer', prompt="Which reverse primer did you use?", required=True, default="GACTACHVGGGTATCTAATCC", help="Which reverse primer did you use?")
@click.option('--mergepairs', prompt="Choose wich method to use for stitching paired reads (vsearch, pandaseq)", required=True, default="vsearch", type=click.Choice(['pandaseq', 'vsearch', 'none']), help="Choose wich method to use for stitching paired reads")
@click.option('--classification', prompt="Choose wich classification option you want to use (sina, stampa, rdp, blast)", required=True, type=click.Choice(['sina', 'stampa', 'rdp', 'blast']), help="Choose wich classification option you want to use")
@click.option('--reference_db', prompt="Choose wich reference database to use (silva, unite)", required=True, type=click.Choice(['silva', 'unite']), help="Choose wich reference database to use")
@click.option('--clustering', prompt="Choose wich clustering method you want to use (usearch_smallmem, swarm)", required=True, default="usearch_smallmem", type=click.Choice(['usearch_smallmem', 'swarm']), help="Choose wich clustering method you want to use")
def make_config(project,config,path,remote, rename, forward_primer, reverse_primer, mergepairs, classification, reference_db, clustering):
    """Write the file `config` and complete the sample names and paths for all files in `path`."""
    represent_dict_order = lambda self, data:  self.represent_mapping('tag:yaml.org,2002:map', data.items())
    yaml.add_representer(OrderedDict, represent_dict_order)
    path = os.path.realpath(path)

    conf = OrderedDict()
    samples = get_sample_files(path, remote)
    if rename:
        renamed = 0
        for line in open(rename):
            sample, newname = line.split()
            if sample in samples:
                newname = newname.replace("_","-")
                samples[newname] = samples.pop(sample)
                renamed += 1

    create_metadata_template("metadata.txt", samples.keys())

    logging.info("Found %d samples under %s" % (len(samples), path if remote == None else "remote project %s " % remote))
    if rename:
        logging.info("Renamed %d samples" % renamed)

    conf["project"] = project
    conf["minsize"] = 2
    conf["adapters_fasta"] = "/data/ngs/adapters/contaminant_list.txt"
    conf["pandaseq_overlap"] = "10"
    conf["pandaseq_quality"] = "25"
    conf["pandaseq_minlength"] = "100"
    conf["pandaseq_maxlength"] = "700"

    conf["quality_control"] = OrderedDict()
    conf["quality_control"]["barcode"] = OrderedDict()
    conf["quality_control"]["barcode"]["threshold"] = 5
    conf["quality_control"]["barcode"]["length"] = 8
    conf["quality_control"]["barcode"]["seperator"] = "#"

    conf["quality_control"]["trimming"] = OrderedDict()
    conf["quality_control"]["trimming"]["quality"] = 25

    conf["forward_primer"] = forward_primer
    conf["reverse_primer"] = reverse_primer

    conf["mergepairs"] = mergepairs
    conf["vsearch_minmergelen"] = "200"
    conf["metadata"] = "metadata.txt"
    if remote != None:
        conf["remote"] = True
    else:
        conf["remote"] = False
    conf["barcode_in_header"]  = False

    conf["its"] = False
    conf["its_region"] = "ITS2"
    conf["clustering"] =  clustering
    conf["classification"] = classification
    conf["use_full_lineage"] = False
    conf["rdp_confidence_cutoff"] = 0.80
    conf["reference_db"] = reference_db

    conf["convert_to_casava1.8"] = False
    conf["data"] = samples

    with open(config, "w") as f:
        print(yaml.dump(conf, default_flow_style=False), file=f)
    logging.info("Configuration file written to %s" % config)

if __name__ == "__main__": 
    make_config()
