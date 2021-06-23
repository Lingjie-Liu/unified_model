#### config file ####
configfile: "config.yml"

comparaReg_dir = config["comparaReg"]
crana_dir = config["comparaReg_analyses"]

#### import packages ####
import os
import pandas as pd
import numpy as np
# for subsetting samples in a multiindex dataframe
idx = pd.IndexSlice

#### metadata ####
metadata_index = ['assay_type', 'biosample_species', 'biosample_type']

def set_metadata_index(metadata):
    metadata.set_index(metadata_index, drop = False, inplace = True)
    metadata.sort_index(inplace = True)
    return metadata

metadata = pd.read_csv("metadata.csv", dtype=str)
metadata = set_metadata_index(metadata)

#### store wildcard for easier access ####
combine_wildcard = "{assay_type}-{biosample_species}-{biosample_type}"
expand_combine_wildcard = "{df.assay_type}-{df.biosample_species}-{df.biosample_type}"

#### rules ####
rule all:
    input:
        expand(os.path.join("results/tidgrng", expand_combine_wildcard + "-{scheme}.RDS"), df = metadata.itertuples(), scheme = 26)

#### load rules ####
include: "rules/unified_model.smk"
