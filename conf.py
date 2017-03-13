import logging
import multiprocessing
import os
import tempfile
import yaml
from collections import OrderedDict

# Adapted from: https://github.com/pnnl/atlas/blob/master/atlas/conf.py


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
                sample_id = sample_id.replace("_", "-").replace(" ", "-")

                fq_path = os.path.join(dir_name, fname)
                fastq_paths = [fq_path]

                if fq_path in seen: continue

                if "_R1" in fname or "_r1" in fname:
                    r2_path = os.path.join(dir_name, fname.replace("_R1", "_R2").replace("_r1", "_r2"))
                    if not r2_path == fq_path:
                        seen.add(r2_path)
                        fastq_paths.append(r2_path)

                if "_R2" in fname or "_r2" in fname:
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
    conf["samples"] = samples

    with open(config, "w") as f:
        print(yaml.dump(conf, default_flow_style=False), file=f)
    logging.info("Configuration file written to %s" % config)

if __name__ == "__main__":
    make_config(config="config.yaml", path="data")
  
