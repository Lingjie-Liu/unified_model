#' @title Plot transcript_quantifier object
#'
#' @description plot transcripts and masks
#'
#' @inheritParams add_data
#' @param gene_name Name of your query gene. Only works if gene ids were given
#' to the \link{transcript_quantifier-class} object
#' @param chrom the chromosome of your query region (single value)
#' @param start the start position of your query region (single value)
#' @param end the end position of your query region (single value)
#' @param strand the strand specificity of your query region (default: ANY)
#' @param ymax_bw maximum value of y-axis for bw files when plotting
#' @param ymax_abundance maximum value of y-axis for abundance when plotting
#' @param all_transcripts whether to show all transcripts in the annotation
#' track or only show one transcript per transcript model. default TRUE, which
#' shows all transcripts.
#'
#' @return plotted tracks of transcripts and masks of your query gene
#'
#' @include transcript_quantifier-class.R
#' @name plot_multiple_samples
#' @export
plot_multiple_samples <- function(tq,
                                   gene_name = NULL,
                                   chrom = NULL,
                                   start = NULL,
                                   end = NULL,
                                   strand = NULL,
                                   bigwig_plus = NULL,
                                   bigwig_minus = NULL,
                                   count_granges = NULL,
                                   ymax_bw = NULL,
                                   ymax_abundance = NULL,
                                   all_transcripts = TRUE,
                                   tx_names_1 = NULL,
                                   tx_names_2 = NULL
) {
    # Check that only gene name or position information is specified
    if (!is.null(gene_name) & (!is.null(chrom) | !is.null(start) |
                               !is.null(end))) {
        stop("Only gene name OR positional information can be specified")
    }

    # Check for gene name specification
    if (!is.null(gene_name)) {
        gid_col <- tq@column_identifiers[2]
        if (is.na(gid_col)) {
            stop("transcript_quantifier object not built with gene ids")
        }
        if (!all(gene_name %in%
                 S4Vectors::elementMetadata(tq@transcripts)[, gid_col])) {
            stop(paste("Gene", gene_name, "does not exist"))
        }
    } else if (!is.null(chrom) &
               is.numeric(start) & is.numeric(end)) {
        if (length(chrom) > 1 | length(start) > 1 | length(end) > 1) {
            stop("Only one value may be specified for each positional option")
        }
        # Check that position overlaps with one or more transcripts
        query_range <-
            GenomicRanges::GRanges(chrom, IRanges::IRanges(start, end),
                                   strand = strand)
        query_inter <-
            GenomicRanges::intersect(query_range, tq@transcripts,
                                     ignore.strand = is.null(strand))
        if (length(query_inter) == 0) {
            stop("Positional query does not intersect any transcripts")
        }
    } else {
        stop("Incorrect positional specification")
    }

    ## ** End checking **
    # Get target transcripts
    target_tx <- DENR:::get_transcripts(tq, gene_name, chrom, start, end, strand)

    # Define bounds of plot range
    if (is.null(start)) {
        chrom <-
            S4Vectors::runValue(droplevels(GenomicRanges::seqnames(target_tx)[1]))
        start <- min(GenomicRanges::start(target_tx))
        end <- max(GenomicRanges::end(target_tx))
    }

    # genome coordination track
    axis_track <- Gviz::GenomeAxisTrack(target_tx)

    tx_col <- tq@column_identifiers[1]
    gid_col <- tq@column_identifiers[2]

    # Transcripts track
    if (!is.na(gid_col)) {
        gene_names <- S4Vectors::elementMetadata(target_tx)[, gid_col]
    } else {
        gene_names <- rep("placeholder", length(target_tx))
    }
    unique_gene_names <- unique(unlist(gene_names))

    # Add columns for gene region track
    target_tx$gene <- unlist(GenomicRanges::values(target_tx)[[gid_col]])
    target_tx$feature <- target_tx$gene
    target_tx$transcript <- unlist(GenomicRanges::values(target_tx)[[tx_col]])

    if (!all_transcripts) {
        # Only show one transcript per transcript model
        tx_key <- tq@transcript_model_key
        tx_key$tx_group_name <-
            paste0("G", tx_key$group, "M", tx_key$model)

        tx_key <- tx_key[!base::duplicated(tx_key$tx_group_name), ]

        target_tx <- target_tx[target_tx$transcript %in% tx_key$tx_name]
        target_tx$tx_group_name <-
            tx_key[base::match(target_tx$transcript,
                               tx_key$tx_name), "tx_group_name"]
    }

    tx_tracks <- list()

    tx_filter_1 <- target_tx$transcript %in% tx_names_1
    tx_filter_2 <- target_tx$transcript %in% tx_names_2
    tx_filter_3 <- !(tx_filter_1 | tx_filter_2)

    for (tx_filter in list(tx_filter_1, tx_filter_2, tx_filter_3)) {
        if (any(tx_filter)) {
            tx_track <- Gviz::GeneRegionTrack(
                target_tx[tx_filter],
                name = "isoforms",
                shape = "arrow"
            )
            tx_tracks <- c(tx_tracks, tx_track)
        }

    }

    # Get masks
    transcripts <- S4Vectors::elementMetadata(target_tx)[, tx_col]
    all_masks <- DENR:::get_masks(tq, transcripts)
    target_masks <- GenomicRanges::reduce(all_masks[[1]])
    add_masks <- all_masks[[2]]
    if (!is.null(add_masks)) add_masks <-
        GenomicRanges::reduce(add_masks, ignore.strand = TRUE)

    # Masks track
    mask_tracks <- list()
    if (length(target_masks) > 0) {
        mask_tracks[["model"]] <- Gviz::AnnotationTrack(target_masks,
                                                        name = "model masks",
                                                        feature = Gviz::strand(target_masks),
                                                        shape = "box",
                                                        chromosome = chrom)
    }
    if (length(add_masks) > 0) {
        mask_tracks[["add"]] <- Gviz::AnnotationTrack(add_masks,
                                                      name = "additional masks",
                                                      shape = "box",
                                                      chromosome = chrom,
                                                      stacking = "dense")
    } else if (!is.null(tq@add_mask) & tq@add_mask_scale) {
        sub_masks <- subsetByOverlaps(tq@add_mask, target_tx)
        mask_tracks[["add"]] <- Gviz::AnnotationTrack(sub_masks,
                                                      name = "additional masks",
                                                      shape = "box",
                                                      chromosome = chrom,
                                                      stacking = "dense")
    }

    counting_track <- Gviz::AnnotationTrack(count_granges,
                                            name = "gene regions",
                                            feature = Gviz::strand(count_granges),
                                            shape = "box",
                                            chromosome = chrom)

    # Some objects for data retrieval
    strand_col <- c(`+` = "blue", `-` = "red")
    bw_max <- 0

    # make DataTrack for bigwig files
    make_data_track_bw <- function(bw_files, bsize = tq@bin_size) {
        data_tracks <- list()
        for (s in names(bw_files)) {
            file <- bw_files[s]
            if (!file.exists(file)) {
                stop(paste("File", file, "does not exist"))
            }
            bw <-
                rtracklayer::import(file,
                                    which = GenomicRanges::GRanges(
                                        seqnames = chrom,
                                        ranges = IRanges::IRanges(start = start,
                                                                  end = end)
                                    ))
            bw$score <- abs(bw$score)

            # Get potential ymax
            bw_max <-
                max(bw_max, max(stats::quantile(abs(bw$score), 0.99))) * 1.05
            if (length(bw) > 0) {
                data_tracks[[s]] <- Gviz::DataTrack(
                    range = bw,
                    type = "h",
                    window = -1,
                    windowSize = bsize,
                    name = paste0("PRO-seq (", s, ")"),
                    col = strand_col[s],
                    strand = s,
                    chromosome = chrom)
            }
        }
        return(list(data_tracks, bw_max))
    }

    data_tracks <- list()

    for (i in 1:length(bigwig_plus)) {
        bw_files <- c(`+` = bigwig_plus[i], `-` = bigwig_minus[i])
        bw_tracks <- make_data_track_bw(bw_files)
        data_tracks <- c(data_tracks, bw_tracks[[1]])
        if (!is.na(bw_tracks[[2]]) & (bw_tracks[[2]] > bw_max)) bw_max <- bw_tracks[[2]]
    }

    # Plot tracks
    args <- list(
        trackList = c(
            list(axis_track),
            tx_tracks,
            mask_tracks,
            list(counting_track),
            data_tracks),
        from = start,
        to = end,
        chromosome = chrom,
        transcriptAnnotation = "transcript"
    )

    # Override default ymax with user specification if given
    if (!is.null(ymax_bw)) {
        bw_max <- ymax_bw
    }

    # if (!is.null(ymax_abundance)) {
    #     abundance_max <- ymax_abundance
    # }

    datatrack_names <- c("Summarized read counts (+)",
                         "Summarized read counts (-)",
                         "PRO-seq (+)", "PRO-seq (-)")

    # Set max for all valid data tracks
    for (track in seq_along(args$trackList)) {
        if (class(args$trackList[[track]]) == "DataTrack") {
            if (args$trackList[[track]]@name %in% datatrack_names) {
                args$trackList[[track]] <-
                    DENR:::set_datatrack_ylim(args$trackList[[track]], c(0, bw_max))
            }
            if (args$trackList[[track]]@name %in%
                c("Predicted abundance (+)", "Predicted abundance (-)")) {
                args$trackList[[track]] <-
                    DENR:::set_datatrack_ylim(args$trackList[[track]],
                                              c(0, abundance_max))
            }
        }
    }

    args$trackList <- #nolint
        args$trackList[!unlist(lapply(args$trackList, is.null))] #nolint

    # Add colors for genes and masks
    gene_colors <- viridisLite::viridis(length(unique_gene_names))
    names(gene_colors) <- unique_gene_names

    mask_colors <- c("blue", "red")
    names(mask_colors) <- c("+", "-")
    extension <- c("extend.left" = 0.15, "extend.right" = 0.02)
    # # Add legend back to the last DataTrack
    # if (length(colnames(S4Vectors::mcols(abundance_1))) > 1) {
    #     dt_idx <- max(which(sapply(args$trackList, class) == "DataTrack"))
    #     suppressMessages(
    #         args$trackList[[dt_idx]] <-
    #             Gviz::setPar(args$trackList[[dt_idx]], "legend", TRUE)
    #     )
    # }

    args <- c(args, gene_colors, mask_colors, extension)

    return(do.call(Gviz::plotTracks, args))
}
