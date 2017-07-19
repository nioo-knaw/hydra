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
        walker = os.walk(path)
    for dir_name, sub_dirs, files in walker:
        for fname in files:
            if ".fastq" in fname or ".fq" in fname:
                sample_id = fname.partition(".fastq")[0]
                if ".fq" in sample_id:
                    sample_id = fname.partition(".fq")[0]

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

@click.command()
@click.option('--project', prompt=True, required=True, help='Give your project a nice name')
@click.option('--config', default="config.yaml", show_default=True, help='File to write the configuration to')
@click.option('--remote', help='Specify a ENA project to use as remote data (for example PRJEB14409')
@click.option('--path', default="../data", show_default=True, help='path to data folder')
def make_config(project,config,path,remote):
    """Write the file `config` and complete the sample names and paths for all files in `path`."""
    represent_dict_order = lambda self, data:  self.represent_mapping('tag:yaml.org,2002:map', data.items())
    yaml.add_representer(OrderedDict, represent_dict_order)
    path = os.path.realpath(path)

    conf = OrderedDict()
    samples = get_sample_files(path, remote)

    logging.info("Found %d samples under %s" % (len(samples), path if remote == None else "remote project %s " % remote))

    conf["project"] = project
    conf["minsize"] 2
    conf["adapters_fasta"] = "/data/ngs/adapters/contaminant_list.txt"
    conf["pandaseq_overlap"] = "10"
    conf["pandaseq_quality"] = "25"
    conf["pandaseq_minlength"] = "100"
    conf["pandaseq_maxlength"] = "700"
    
    
    conf["forward_primer"] = "CCTACGGGNGGCWGCAG"
    conf["reverse_primer"] = "GACTACHVGGGTATCTAATCC"

    conf["silva_arb"] = "/data/db/Silva/128/SSURef_NR99_128_SILVA_07_09_16_opt.arb"
    conf["mergepairs"] = "vsearch"  
    conf["metadata"] = "metadata.txt"
    if remote != None:
        conf["remote"] = True
    else:
        conf["remote"] = False

    conf["its"] = False 
    conf["clustering"] =  "usearch_smallmem"
    conf["classification"] = "stampa"

    conf["stampa_db"] = "/data/db/unite/itsx.ITS2.stampa.fasta"

    conf["convert_to_casava1.8"] = False
    conf["data"] = samples

    with open(config, "w") as f:
        print(yaml.dump(conf, default_flow_style=False), file=f)
    logging.info("Configuration file written to %s" % config)

if __name__ == "__main__": 
    make_config()
