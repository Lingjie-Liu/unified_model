rule find_dominant_promoter_and_TSS:
    input:
        tq = os.path.join(crana_dir, "results/denr/combine/tq/human_rhesus", "template-{scheme}.RDS"),
        texp = os.path.join(crana_dir, "results/denr/combine/texp/human_rhesus", combine_wildcard + "-{scheme}.csv"),
        non_olp_gn = os.path.join(comparaReg_dir, "external_resources/human/non_overlapping_coding_genes.csv")
    params:
        texp_cutoff = 1, # filter out genes with dominant promoter with expressions lower than this cutoff
        tid_cutoff = 1000 # filter out TSSs region from dominant promoter larger than this distance
    threads:1
    log:
        os.path.join("logs/find_dominant_promoter_and_TSS", combine_wildcard + "-{scheme}.log")
    output:
        tid = os.path.join("results/tidgrng", combine_wildcard + "-{scheme}.RDS")
    script:
        "../scripts/find_dominant_promoter_and_TSS.R"

rule analyze_two_samples:
    input:
        tq = os.path.join(crana_dir, "results/denr/combine/tq/human_rhesus", "template-{scheme}.RDS"),
        tid1 = os.path.join("results/tidgrng", ind1_wildcard + "-{scheme}.RDS"),
        tid2 = os.path.join("results/tidgrng", ind2_wildcard + "-{scheme}.RDS"),
        bwp1_p5 = os.path.join(crana_dir, "results/denr/combine/bigwig/p5/human_rhesus", ind1_wildcard + "_plus.bw"),
        bwm1_p5 = os.path.join(crana_dir, "results/denr/combine/bigwig/p5/human_rhesus", ind1_wildcard + "_minus.bw"),
        bwp1_p3 = os.path.join(crana_dir, "results/denr/combine/bigwig/p3/human_rhesus", ind1_wildcard + "_plus.bw"),
        bwm1_p3 = os.path.join(crana_dir, "results/denr/combine/bigwig/p3/human_rhesus", ind1_wildcard + "_minus.bw"),
        bwp2_p5 = os.path.join(crana_dir, "results/denr/combine/bigwig/p5/human_rhesus", ind2_wildcard + "_plus.bw"),
        bwm2_p5 = os.path.join(crana_dir, "results/denr/combine/bigwig/p5/human_rhesus", ind2_wildcard + "_minus.bw"),
        bwp2_p3 = os.path.join(crana_dir, "results/denr/combine/bigwig/p3/human_rhesus", ind2_wildcard + "_plus.bw"),
        bwm2_p3 = os.path.join(crana_dir, "results/denr/combine/bigwig/p3/human_rhesus", ind2_wildcard + "_minus.bw")
    params:
        tid_cutoff = 1000, # length cutoff for how long a TID containing multiple TSSs could span
        tsn_cutoff = 5, # cutoff of minimum number of reads for a TSN
        pause_cutoff = 250,
        gb_min_length = 1e4,
        tts_length = 250, # parameter m
        quantile_normalization = "{normalization}", # "identity", "qnorm"
        result_dir = os.path.join("results/between_samples", ind1_wildcard + "_vs_" + ind2_wildcard)
    log:
        os.path.join("logs/analyze_two_samples", ind1_wildcard + "_vs_" + ind2_wildcard, "S{scheme}-{normalization}.log")
    output:
        alpha = os.path.join("results/between_samples", ind1_wildcard + "_vs_" + ind2_wildcard, "S{scheme}-{normalization}", "alpha.csv"),
        beta = os.path.join("results/between_samples", ind1_wildcard + "_vs_" + ind2_wildcard, "S{scheme}-{normalization}", "beta.csv")
    script:
        "../scripts/analyze_two_samples.R"
