rule find_dominant_promoter_and_TSS:
    input:
        tq = os.path.join(crana_dir + "results/denr/combine/tq/human_rhesus", "template-{scheme}.RDS"),
        texp = os.path.join(crana_dir + "results/denr/combine/texp/human_rhesus", combine_wildcard + "-{scheme}.csv")
    params:
        texp_cutoff = 10, # filter out genes with dominant promoter with expressions lower than this cutoff
        tid_cutoff = 1000 # filter out TSSs region from dominant promoter larger than this distance
    threads:1
    log:
        os.path.join("logs/find_dominant_promoter_and_TSS", combine_wildcard + ".log")
    output:
        tid = os.path.join("results/tidgrng", combine_wildcard + "-{scheme}.RDS")
    script:
        "../scripts/find_dominant_promoter_and_TSS.R"
