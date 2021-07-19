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

metadata_comparison = pd.read_csv("metadata_comparison.csv", dtype=str)

#### store wildcard for easier access ####
combine_wildcard = "{assay_type}-{biosample_species}-{biosample_type}"
expand_combine_wildcard = "{df.assay_type}-{df.biosample_species}-{df.biosample_type}"
ind1_wildcard = "{assay_type_1}-{biosample_species_1}-{biosample_type_1}"
ind2_wildcard = "{assay_type_2}-{biosample_species_2}-{biosample_type_2}"

#### rules ####
rule all:
    input:
        expand(os.path.join("results/tidgrng", expand_combine_wildcard + "-{scheme}.RDS"), df = metadata.itertuples(), scheme = 26),
        expand(os.path.join("results/between_samples", "{df.assay_type_1}-{df.biosample_species_1}-{df.biosample_type_1}" + "_vs_" + "{df.assay_type_2}-{df.biosample_species_2}-{df.biosample_type_2}", "S{scheme}-{normalization}", "alpha.csv"), df = metadata_comparison, scheme = 26, normalization = ["identity", "qnorm"])

#### load rules ####
include: "rules/unified_model.smk"
