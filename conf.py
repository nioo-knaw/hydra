import logging
import multiprocessing
import re
import os
import tempfile
import yaml
from collections import OrderedDict

# Adapted from: https://github.com/pnnl/atlas/blob/master/atlas/conf.py

# http://stackoverflow.com/a/3675423
def replace_last(source_string, replace_what, replace_with):
    head, _sep, tail =  source_string.rpartition(replace_what)
    if _sep == '':
        return tail
    else: 
        return head + replace_with + tail

def get_sample_files(path):
    samples = OrderedDict()
    seen = set()
    for dir_name, sub_dirs, files in os.walk(path):
        for fname in files:

            if ".fastq" in fname or ".fq" in fname:

                sample_id = fname.partition(".fastq")[0]
                if ".fq" in sample_id:
                    sample_id = fname.partition(".fq")[0]

                sample_id = sample_id.replace("_R1", "").replace("_r1", "").replace("_R2", "").replace("_r2", "")
                sample_id = re.sub("_1$", "", sample_id)
                sample_id = re.sub("_2$", "", sample_id)
#                sample_id = replace_last(sample_id, "_1", "")
#                sample_id = replace_last(sample_id, "_2", "")
                sample_id = sample_id.replace("_", "-").replace(" ", "-")
#                print(sample_id)

                fq_path = os.path.join(dir_name, fname)
                fastq_paths = [fq_path]

                if fq_path in seen: continue

                print(fname)
                if "_R1" in fname or "_r1" in fname or "_1" in fname:
                    fname = replace_last(fname,"_1.","_2.")
                    r2_path = os.path.join(dir_name, fname.replace("_R1", "_R2").replace("_r1", "_r2"))
                    print(r2_path)
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


def make_config(config, path):
    """Write the file `config` and complete the sample names and paths for all files in `path`."""
    represent_dict_order = lambda self, data:  self.represent_mapping('tag:yaml.org,2002:map', data.items())
    yaml.add_representer(OrderedDict, represent_dict_order)
    path = os.path.realpath(path)

    conf = OrderedDict()
    samples = get_sample_files(path)

    logging.info("Found %d samples under %s" % (len(samples), path))
    conf["project"] = "My-Project"
    conf["adapters_fasta"] = "/data/ngs/adapters/contaminant_list.txt"
    conf["pandaseq_overlap"] = "10"
    conf["pandaseq_quality"] = "25"
    conf["pandaseq_minlength"] = "100"
    conf["pandaseq_maxlength"] = "700"

    conf["forward_primer"] = "CCTACGGGNGGCWGCAG"
    conf["reverse_primer"] = "GACTACHVGGGTATCTAATCC"

    conf["silva_arb"] = "/data/db/Silva/128/SSURef_NR99_128_SILVA_07_09_16_opt.arb"
    conf["metadata"] = "data/metadata.txt"

    conf["data"] = samples

    with open(config, "w") as f:
        print(yaml.dump(conf, default_flow_style=False), file=f)
    logging.info("Configuration file written to %s" % config)

if __name__ == "__main__":
    make_config(config="config.yaml", path="data")
  
