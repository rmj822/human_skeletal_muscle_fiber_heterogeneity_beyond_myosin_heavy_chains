################################################################################################################################################
#################################################     Panel A  ##############################################################
################################################################################################################################################


# Load raw data -----------------------------------------------------------

data_raw <- vroom::vroom(here::here("data-raw/muscle_disease_Report_ProteinQuant.tsv"))

metadata <- vroom::vroom(here::here("data/metadata_MD_w_fiber_type_w_anonim.csv")) |>
    dplyr::rename("fiberID" = 1)

contaminants <- data_raw |>
    dplyr::filter(grepl("KRT", PG.Genes)) |>
    dplyr::pull(PG.Genes)

data_raw <- data_raw |>
    dplyr::filter(!PG.Genes %in% contaminants)

old_column_names <- colnames(data_raw)

new_column_names <- gsub(
    pattern = ".*DIA_S",
    replacement = "",
    old_column_names
)

data_wrangled <- data_raw

new_column_names <- gsub(
    pattern = "_.*",
    replacement = "",
    new_column_names
)

BiocGenerics::colnames(data_wrangled) <- new_column_names

data_wrangled <- data_wrangled |>
    dplyr::select(!PG.ProteinAccessions) |>
    dplyr::select(!PG.Genes) |>
    dplyr::select(!PG.ProteinDescriptions) |>
    tibble::column_to_rownames("PG.ProteinGroups") |>
    t() |>
    as.data.frame() |>
    tibble::rownames_to_column("fiber_number") |>
    dplyr::mutate(fiber_number = as.numeric(fiber_number)) |>
    dplyr::inner_join(metadata |>
                          dplyr::select(
                              fiberID, fiber_number
                          )) |>
    dplyr::select(!fiber_number) |>
    tibble::column_to_rownames("fiberID") |>
    t() |>
    as.data.frame()

# Filtering genes with less than 70% valid --------------------------------

#' Filtering missing values from rows
#'
#' @param .data dataset to filter
#' @param percentage_accepted_missing % of accepted missing values
#'
#' @return a dataframe
#' @export
#'
#' @examples
filtering_rows_Na <- function(.data, percentage_accepted_missing) {
    row_keep_vector <- .data |>
        is.na() |>
        rowSums()

    row_keep_vector <- row_keep_vector / ncol(.data)

    row_keep_vector <- row_keep_vector <= percentage_accepted_missing

    data_filtered <- .data |>
        tibble::add_column(row_keep_vector) |>
        dplyr::filter(row_keep_vector == T) |>
        dplyr::select(!starts_with("row"))

    return(data_filtered)
}
data_row_filtered <- data_wrangled |>
    filtering_rows_Na(0.3) |>
    tibble::rownames_to_column("protein_groups")

write.csv(data_row_filtered,
          here::here("data/lnc_quant/muscle_disease/data_row_filtered.csv"))

data <- data_wrangled |>
    log2() |>
    as.data.frame() |>
    tibble::rownames_to_column("protein_groups") |>
    dplyr::filter(grepl("ENS",
                        protein_groups))

muscle_disease_data <- data_wrangled |>
    log2() |>
    as.data.frame()

metadata_MD <- metadata

# Doing pseudobulk on fast and slow fibers --------------------------------

pseudobulk_maker_fiber_type <- function(.data, metadata, subject_id, grouping, grouping_2, colname) {
    selection_vector <- metadata |>
        dplyr::filter(condition == grouping) |>
        dplyr::filter(fiber_type == grouping_2) |>
        dplyr::filter(subject == subject_id) |>
        dplyr::pull(fiberID)

    pseudobulk_median_calculation <- .data |>
        as.data.frame() |>
        dplyr::select(selection_vector) |>
        t() |>
        as.data.frame() |>
        dplyr::mutate(dplyr::across(.cols = everything(), median, na.rm = T)) |>
        unique() |>
        t()

    colnames(pseudobulk_median_calculation) <- colname
    return(pseudobulk_median_calculation)
}

data_pseudobulk_ft <- data.frame(
    pseudobulk_maker_fiber_type(
        .data = muscle_disease_data,
        metadata = metadata_MD,
        grouping = "control",
        grouping_2 = "Type 1",
        subject_id = "C1",
        colname = "C1_T1"
    ),
    pseudobulk_maker_fiber_type(
        .data = muscle_disease_data,
        metadata = metadata_MD,
        grouping = "control",
        grouping_2 = "Type 1",
        subject_id = "C2",
        colname = "C2_T1"
    ),
    pseudobulk_maker_fiber_type(
        .data = muscle_disease_data,
        metadata = metadata_MD,
        grouping = "control",
        grouping_2 = "Type 1",
        subject_id = "C3",
        colname = "C3_T1"
    ),
    pseudobulk_maker_fiber_type(
        .data = muscle_disease_data,
        metadata = metadata_MD,
        grouping = "control",
        grouping_2 = "Type 2A",
        subject_id = "C1",
        colname = "C1_T2A"
    ),
    pseudobulk_maker_fiber_type(
        .data = muscle_disease_data,
        metadata = metadata_MD,
        grouping = "control",
        grouping_2 = "Type 2A",
        subject_id = "C2",
        colname = "C2_T2A"
    ),
    pseudobulk_maker_fiber_type(
        .data = muscle_disease_data,
        metadata = metadata_MD,
        grouping = "control",
        grouping_2 = "Type 2A",
        subject_id = "C3",
        colname = "C3_T2A"
    ),
    pseudobulk_maker_fiber_type(
        .data = muscle_disease_data,
        metadata = metadata_MD,
        grouping = "ACTA1_nemaline_myopaty",
        grouping_2 = "Type 1",
        subject_id = "A1",
        colname = "A1_T1"
    ),
    pseudobulk_maker_fiber_type(
        .data = muscle_disease_data,
        metadata = metadata_MD,
        grouping = "ACTA1_nemaline_myopaty",
        grouping_2 = "Type 1",
        subject_id = "A2",
        colname = "A2_T1"
    ),
    pseudobulk_maker_fiber_type(
        .data = muscle_disease_data,
        metadata = metadata_MD,
        grouping = "ACTA1_nemaline_myopaty",
        grouping_2 = "Type 1",
        subject_id = "A3",
        colname = "A3_T1"
    ),
    pseudobulk_maker_fiber_type(
        .data = muscle_disease_data,
        metadata = metadata_MD,
        grouping = "ACTA1_nemaline_myopaty",
        grouping_2 = "Type 2A",
        subject_id = "A1",
        colname = "A1_T2A"
    ),
    pseudobulk_maker_fiber_type(
        .data = muscle_disease_data,
        metadata = metadata_MD,
        grouping = "ACTA1_nemaline_myopaty",
        grouping_2 = "Type 2A",
        subject_id = "A2",
        colname = "A2_T2A"
    ),
    pseudobulk_maker_fiber_type(
        .data = muscle_disease_data,
        metadata = metadata_MD,
        grouping = "ACTA1_nemaline_myopaty",
        grouping_2 = "Type 2A",
        subject_id = "A3",
        colname = "A3_T2A"
    ),
    pseudobulk_maker_fiber_type(
        .data = muscle_disease_data,
        metadata = metadata_MD,
        grouping = "TNNT1_nemaline_myopaty",
        grouping_2 = "Type 1",
        subject_id = "T1",
        colname = "T1_T1"
    ),
    pseudobulk_maker_fiber_type(
        .data = muscle_disease_data,
        metadata = metadata_MD,
        grouping = "TNNT1_nemaline_myopaty",
        grouping_2 = "Type 1",
        subject_id = "T2",
        colname = "T2_T1"
    ),
    pseudobulk_maker_fiber_type(
        .data = muscle_disease_data,
        metadata = metadata_MD,
        grouping = "TNNT1_nemaline_myopaty",
        grouping_2 = "Type 1",
        subject_id = "T3",
        colname = "T3_T1"
    ),
    pseudobulk_maker_fiber_type(
        .data = muscle_disease_data,
        metadata = metadata_MD,
        grouping = "TNNT1_nemaline_myopaty",
        grouping_2 = "Type 2A",
        subject_id = "T1",
        colname = "T1_T2A"
    ),
    pseudobulk_maker_fiber_type(
        .data = muscle_disease_data,
        metadata = metadata_MD,
        grouping = "TNNT1_nemaline_myopaty",
        grouping_2 = "Type 2A",
        subject_id = "T2",
        colname = "T2_T2A"
    ),
    pseudobulk_maker_fiber_type(
        .data = muscle_disease_data,
        metadata = metadata_MD,
        grouping = "TNNT1_nemaline_myopaty",
        grouping_2 = "Type 2A",
        subject_id = "T3",
        colname = "T3_T2A"
    )
)

data_limma <- data_pseudobulk_ft

limma::plotDensities(data_limma)

data_limma <- limma::normalizeBetweenArrays(data_limma, method = "quantile")

limma::plotDensities(data_limma)

data_grouping_ft <- data.frame(
    "fiber_ID" = colnames(data_limma),
    "condition" = c(
        "control",
        "control",
        "control",
        "control",
        "control",
        "control",
        "actin",
        "actin",
        "actin",
        "actin",
        "actin",
        "actin",
        "troponin",
        "troponin",
        "troponin",
        "troponin",
        "troponin",
        "troponin"
    ),
    "fiber_type" = c(
        "type_1",
        "type_1",
        "type_1",
        "type_2A",
        "type_2A",
        "type_2A",
        "type_1",
        "type_1",
        "type_1",
        "type_2A",
        "type_2A",
        "type_2A",
        "type_1",
        "type_1",
        "type_1",
        "type_2A",
        "type_2A",
        "type_2A"
    ),
    "subject" = c(
        "C1",
        "C2",
        "C3",
        "C1",
        "C2",
        "C3",
        "A1",
        "A2",
        "A3",
        "A1",
        "A2",
        "A3",
        "T1",
        "T2",
        "T3",
        "T1",
        "T2",
        "T3"
    )
)

data_grouping_ft <- data_grouping_ft |>
    tidyr::unite(
        col = "grouping",
        c(condition, fiber_type),
        sep = "_"
    )

design_matrix <-
    model.matrix(
        ~ 0 + grouping + subject,
        data_grouping_ft
    )

colnames(design_matrix) <- c(
    "actin_type1",
    "actin_type2A",
    "control_type1",
    "control_type2A",
    "troponin_type1",
    "troponin_type2A",
    "A2",
    "A3",
    "C1",
    "C2",
    "C3",
    "T1",
    "T2",
    "T3"

)

fit <- limma::lmFit(
    data_limma,
    design_matrix
)

contrast_matrix <- limma::makeContrasts(
    "control_vs_actin" = (control_type1 + control_type2A)/2 - (actin_type1 + actin_type2A)/2,
    "control_vs_troponin" = (control_type1 + control_type2A)/2 - (troponin_type1 + troponin_type2A)/2,
    "actin_vs_troponin" = (actin_type1 + actin_type2A)/2 - (troponin_type1 + troponin_type2A)/2,
    "T1_control_vs_actin" = control_type1 - actin_type1,
    "T1_control_vs_troponin" = control_type1 - troponin_type1,
    "T1_actin_vs_troponin" = actin_type1 - troponin_type1,
    "T2_control_vs_actin" = control_type2A - actin_type2A,
    "T2_control_vs_troponin" = control_type2A - troponin_type2A,
    "T2_actin_vs_troponin" = actin_type2A - troponin_type2A,
    "main_effect_fiber_type" = (control_type1 + actin_type1 + troponin_type1)/3 - (control_type2A + actin_type2A + troponin_type2A),
    levels = design_matrix
)

tmp <- limma::contrasts.fit(
    fit,
    contrast_matrix
)

tmp <- limma::eBayes(tmp)

# Control - actin ---------------------------------------------------------

DE_results_control_actin <- limma::topTable(tmp,
                                            coef = "control_vs_actin",
                                            sort.by = "P",
                                            n = Inf
) |>
    tibble::rownames_to_column("PG.ProteinGroups") |>
    dplyr::filter(stringr::str_detect(PG.ProteinGroups, "^ENS")) |>
    dplyr::rename(Genes = PG.ProteinGroups) |>
    dplyr::mutate(Genes = gsub(
        pattern = ":.*",
        replacement = "",
        Genes)) |>
    dplyr::mutate("xiao" = 10^-(sqrt(log10(1 / (P.Value^logFC))^2)))


# Control - troponin ------------------------------------------------------

DE_results_control_troponin <- limma::topTable(tmp,
                                               coef = "control_vs_troponin",
                                               sort.by = "P",
                                               n = Inf
) |>
    tibble::rownames_to_column("PG.ProteinGroups") |>
    dplyr::filter(stringr::str_detect(PG.ProteinGroups, "^ENS")) |>
    dplyr::rename(Genes = PG.ProteinGroups) |>
    dplyr::mutate(Genes = gsub(
        pattern = ":.*",
        replacement = "",
        Genes)) |>
    dplyr::mutate("xiao" = 10^-(sqrt(log10(1 / (P.Value^logFC))^2)))


# Actin - troponin --------------------------------------------------------

DE_results_actin_troponin <- limma::topTable(tmp,
                                             coef = "actin_vs_troponin",
                                             sort.by = "P",
                                             n = Inf
) |>
    tibble::rownames_to_column("PG.ProteinGroups") |>
    dplyr::filter(stringr::str_detect(PG.ProteinGroups, "^ENS")) |>
    dplyr::rename(Genes = PG.ProteinGroups) |>
    dplyr::mutate(Genes = gsub(
        pattern = ":.*",
        replacement = "",
        Genes)) |>
    dplyr::mutate("xiao" = 10^-(sqrt(log10(1 / (P.Value^logFC))^2)))


# Fiber type main effect --------------------------------------------------

DE_results_main_effect_ft <- limma::topTable(tmp,
                                             coef = "main_effect_fiber_type",
                                             sort.by = "P",
                                             n = Inf
) |>
    tibble::rownames_to_column("PG.ProteinGroups") |>
    dplyr::filter(stringr::str_detect(PG.ProteinGroups, "^ENS")) |>
    dplyr::rename(Genes = PG.ProteinGroups) |>
    dplyr::mutate(Genes = gsub(
        pattern = ":.*",
        replacement = "",
        Genes)) |>
    dplyr::mutate("xiao" = 10^-(sqrt(log10(1 / (P.Value^logFC))^2)))


# Type 1s control - actin -------------------------------------------------

DE_results_type1_control_actin <- limma::topTable(tmp,
                                                  coef = "T1_control_vs_actin",
                                                  sort.by = "P",
                                                  n = Inf
) |>
    tibble::rownames_to_column("PG.ProteinGroups") |>
    dplyr::filter(stringr::str_detect(PG.ProteinGroups, "^ENS")) |>
    dplyr::rename(Genes = PG.ProteinGroups) |>
    dplyr::mutate(Genes = gsub(
        pattern = ":.*",
        replacement = "",
        Genes)) |>
    dplyr::mutate("xiao" = 10^-(sqrt(log10(1 / (P.Value^logFC))^2)))

# Type 1s control - troponin -------------------------------------------------

DE_results_type1_control_troponin <- limma::topTable(tmp,
                                                     coef = "T1_control_vs_troponin",
                                                     sort.by = "P",
                                                     n = Inf
) |>
    tibble::rownames_to_column("PG.ProteinGroups") |>
    dplyr::filter(stringr::str_detect(PG.ProteinGroups, "^ENS")) |>
    dplyr::rename(Genes = PG.ProteinGroups) |>
    dplyr::mutate(Genes = gsub(
        pattern = ":.*",
        replacement = "",
        Genes)) |>
    dplyr::mutate("xiao" = 10^-(sqrt(log10(1 / (P.Value^logFC))^2)))

# Type 1s actin - troponin -------------------------------------------------

DE_results_type1_actin_troponin <- limma::topTable(tmp,
                                                   coef = "T1_actin_vs_troponin",
                                                   sort.by = "P",
                                                   n = Inf
) |>
    tibble::rownames_to_column("PG.ProteinGroups") |>
    dplyr::filter(stringr::str_detect(PG.ProteinGroups, "^ENS")) |>
    dplyr::rename(Genes = PG.ProteinGroups) |>
    dplyr::mutate(Genes = gsub(
        pattern = ":.*",
        replacement = "",
        Genes)) |>
    dplyr::mutate("xiao" = 10^-(sqrt(log10(1 / (P.Value^logFC))^2)))

# Type 2s control - actin -------------------------------------------------

DE_results_type2_control_actin <- limma::topTable(tmp,
                                                  coef = "T2_control_vs_actin",
                                                  sort.by = "P",
                                                  n = Inf
) |>
    tibble::rownames_to_column("PG.ProteinGroups") |>
    dplyr::filter(stringr::str_detect(PG.ProteinGroups, "^ENS")) |>
    dplyr::rename(Genes = PG.ProteinGroups) |>
    dplyr::mutate(Genes = gsub(
        pattern = ":.*",
        replacement = "",
        Genes)) |>
    dplyr::mutate("xiao" = 10^-(sqrt(log10(1 / (P.Value^logFC))^2)))

# Type 2s control - troponin -------------------------------------------------

DE_results_type2_control_troponin <- limma::topTable(tmp,
                                                     coef = "T2_control_vs_troponin",
                                                     sort.by = "P",
                                                     n = Inf
) |>
    tibble::rownames_to_column("PG.ProteinGroups") |>
    dplyr::filter(stringr::str_detect(PG.ProteinGroups, "^ENS")) |>
    dplyr::rename(Genes = PG.ProteinGroups) |>
    dplyr::mutate(Genes = gsub(
        pattern = ":.*",
        replacement = "",
        Genes)) |>
    dplyr::mutate("xiao" = 10^-(sqrt(log10(1 / (P.Value^logFC))^2)))

# Type 2s actin - troponin -------------------------------------------------

DE_results_type2_actin_troponin <- limma::topTable(tmp,
                                                   coef = "T2_actin_vs_troponin",
                                                   sort.by = "P",
                                                   n = Inf
) |>
    tibble::rownames_to_column("PG.ProteinGroups") |>
    dplyr::filter(stringr::str_detect(PG.ProteinGroups, "^ENS")) |>
    dplyr::rename(Genes = PG.ProteinGroups) |>
    dplyr::mutate(Genes = gsub(
        pattern = ":.*",
        replacement = "",
        Genes)) |>
    dplyr::mutate("xiao" = 10^-(sqrt(log10(1 / (P.Value^logFC))^2)))


output <- DE_results_type1_control_actin |>
    dplyr::mutate("comparison" = "type1_control_vs_actin") |>
    dplyr::bind_rows(
        DE_results_type1_control_troponin |>
            dplyr::mutate("comparison" = "type1_control_vs_troponin"),
        DE_results_type1_actin_troponin |>
            dplyr::mutate("comparison" = "type1_actin_vs_troponin"),
        DE_results_type2_control_actin |>
            dplyr::mutate("comparison" = "type2_control_vs_actin"),
        DE_results_type2_control_troponin |>
            dplyr::mutate("comparison" = "type2_control_vs_troponin"),
        DE_results_type2_actin_troponin |>
            dplyr::mutate("comparison" = "type2_actin_vs_troponin")
    ) |>
    dplyr::mutate(
        significant = xiao < 0.05
    )

# Combined plot for manuscript --------------------------------------------

data_lnc_MD <- data_pseudobulk_ft |>
    as.data.frame() |>
    tibble::rownames_to_column("PG.ProteinGroups") |>
    dplyr::filter(stringr::str_detect(PG.ProteinGroups, "^ENS")) |>
    dplyr::rename(Genes = PG.ProteinGroups) |>
    dplyr::mutate(Genes = gsub(
        pattern = ":.*",
        replacement = "",
        Genes)) |>
    tibble::column_to_rownames("Genes") |>
    t() |>
    as.data.frame() |>
    tibble::rownames_to_column("fiber_ID") |>
    dplyr::inner_join(
        data_grouping_ft) |>
    dplyr::mutate(grouping = dplyr::case_when(
        grouping == "control_type_1" ~ "control_Type 1",
        grouping == "control_type_2A" ~ "control_Type 2A",
        grouping == "actin_type_1" ~ "actin_Type 1",
        grouping == "actin_type_2A" ~ "actin_Type 2A",
        grouping == "troponin_type_1" ~ "troponin_Type 1",
        TRUE ~ "troponin_Type 2A"
    )) |>
    tidyr::separate(col = grouping, into = c("condition", "fiber_type"), sep = "_")

# ENSG00000177822_TR4_ORF9
ENSG00000177822_TR4_ORF9 <- data_lnc_MD |>
    dplyr::select("ENSG00000177822_TR4_ORF9",
                  condition, fiber_type) |>
    dplyr::mutate(condition = factor(condition, levels = c("control", "actin", "troponin"))) |>
    ggplot2::ggplot(ggplot2::aes(x = condition,
                                 y = ENSG00000177822_TR4_ORF9)) +
    ggplot2::geom_boxplot(ggplot2::aes(fill = condition),
                          alpha = 0.85) +
    ggplot2::geom_point(ggplot2::aes(x = condition, y = ENSG00000177822_TR4_ORF9, fill = condition), shape = 21, stroke = 0.5, color = "black") +
    # ggpubr::geom_bracket(xmin = 1, xmax = 2, y.position = 1.25, label = paste("Xiao = ", round(
    #     DE_results_type1_control_actin |>
    #         dplyr::filter(grepl("ENSG00000215483_TR14_ORF67",Genes)) |>
    #         dplyr::pull(xiao),
    #     3
    # )),
    # label.size = 2) +
    # ggpubr::geom_bracket(xmin = 1, xmax = 3, y.position = 17.1, label = paste("Adj.P.Val = ", round(
    #     DE_control_v_troponin |>
    #         dplyr::filter(grepl("ANXA6", Genes)) |>
    #         dplyr::pull(adj.P.Val),
    #     3
    # )),
    # label.size = 2) +
    ggplot2::scale_fill_manual(values=c("#969594",
                                        "#E48C2AFF",
                                        "#d662c4")) +
    ggplot2::ggtitle("ENSG00000177822_TR4_ORF9") +
    ggplot2::labs(
        y = "LFQ intensities (log2)",
        x = ""
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
        legend.position = "none",
        axis.text.x = ggplot2::element_text(size = 6,
                                            vjust = 1,
                                            hjust=1),
        text = ggplot2::element_text(face="bold",
                                     colour="black",
                                     size=7),
        plot.title = ggplot2::element_text(face = "bold",
                                           size = 7,
                                           hjust = 0.5),
        plot.margin = grid::unit(c(1,1,1,1), "mm"),
        strip.text = ggplot2::element_text(size = 6),
        # axis.title.y = ggplot2::element_blank()
    ) +
    ggplot2::facet_grid(~fiber_type)
# ggplot2::ylim(14.8, 17.5)

# ENSG00000215483_TR14_ORF67
ENSG00000215483_TR14_ORF67 <- data_lnc_MD |>
    dplyr::select("ENSG00000215483_TR14_ORF67",
                  condition, fiber_type) |>
    dplyr::mutate(condition = factor(condition, levels = c("control", "actin", "troponin"))) |>
    ggplot2::ggplot(ggplot2::aes(x = condition,
                                 y = ENSG00000215483_TR14_ORF67)) +
    ggplot2::geom_boxplot(ggplot2::aes(fill = condition),
                          alpha = 0.85) +
    ggplot2::geom_point(ggplot2::aes(x = condition, y = ENSG00000215483_TR14_ORF67, fill = condition), shape = 21, stroke = 0.5, color = "black") +
    # ggpubr::geom_bracket(xmin = 1, xmax = 2, y.position = 16.7, label = paste("Adj.P.Val = ", round(
    #     DE_control_v_actin |>
    #         dplyr::filter(grepl("ANXA6",Genes)) |>
    #         dplyr::pull(adj.P.Val),
    #     3
    # )),
    # label.size = 2) +
    # ggpubr::geom_bracket(xmin = 1, xmax = 3, y.position = 17.1, label = paste("Adj.P.Val = ", round(
    #     DE_control_v_troponin |>
    #         dplyr::filter(grepl("ANXA6", Genes)) |>
    #         dplyr::pull(adj.P.Val),
    #     3
    # )),
    # label.size = 2) +
    ggplot2::scale_fill_manual(values=c("#969594",
                                        "#E48C2AFF",
                                        "#d662c4")) +
    ggplot2::ggtitle("ENSG00000215483_TR14_ORF67") +
    ggplot2::labs(
        y = "LFQ intensities (log2)",
        x = ""
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
        legend.position = "none",
        axis.text.x = ggplot2::element_text(size = 6,
                                            vjust = 1,
                                            hjust=1),
        text = ggplot2::element_text(face="bold",
                                     colour="black",
                                     size=7),
        plot.title = ggplot2::element_text(face = "bold",
                                           size = 7,
                                           hjust = 0.5),
        plot.margin = grid::unit(c(1,1,1,1), "mm"),
        strip.text = ggplot2::element_text(size = 6),
        axis.title.y = ggplot2::element_blank()
    ) +
    ggplot2::facet_grid(~fiber_type)
# ggplot2::ylim(14.8, 17.5)

# ENSG00000229425_TR25_ORF40
ENSG00000229425_TR25_ORF40 <- data_lnc_MD |>
    dplyr::select("ENSG00000229425_TR25_ORF40",
                  condition, fiber_type) |>
    dplyr::mutate(condition = factor(condition, levels = c("control", "actin", "troponin"))) |>
    ggplot2::ggplot(ggplot2::aes(x = condition,
                                 y = ENSG00000229425_TR25_ORF40)) +
    ggplot2::geom_boxplot(ggplot2::aes(fill = condition),
                          alpha = 0.85) +
    ggplot2::geom_point(ggplot2::aes(x = condition, y = ENSG00000229425_TR25_ORF40, fill = condition), shape = 21, stroke = 0.5, color = "black") +
    # ggpubr::geom_bracket(xmin = 1, xmax = 2, y.position = 16.7, label = paste("Adj.P.Val = ", round(
    #     DE_control_v_actin |>
    #         dplyr::filter(grepl("ANXA6",Genes)) |>
    #         dplyr::pull(adj.P.Val),
    #     3
    # )),
    # label.size = 2) +
    # ggpubr::geom_bracket(xmin = 1, xmax = 3, y.position = 17.1, label = paste("Adj.P.Val = ", round(
    #     DE_control_v_troponin |>
    #         dplyr::filter(grepl("ANXA6", Genes)) |>
    #         dplyr::pull(adj.P.Val),
    #     3
    # )),
    # label.size = 2) +
    ggplot2::scale_fill_manual(values=c("#969594",
                                        "#E48C2AFF",
                                        "#d662c4")) +
    ggplot2::ggtitle("ENSG00000229425_TR25_ORF40") +
    ggplot2::labs(
        y = "LFQ intensities (log2)",
        x = ""
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
        legend.position = "none",
        axis.text.x = ggplot2::element_text(size = 6,
                                            vjust = 1,
                                            hjust=1),
        text = ggplot2::element_text(face="bold",
                                     colour="black",
                                     size=7),
        plot.title = ggplot2::element_text(face = "bold",
                                           size = 7,
                                           hjust = 0.5),
        plot.margin = grid::unit(c(1,1,1,1), "mm"),
        strip.text = ggplot2::element_text(size = 6),
        axis.title.y = ggplot2::element_blank()
    ) +
    ggplot2::facet_grid(~fiber_type)
# ggplot2::ylim(14.8, 17.5)

# ENSG00000231312_TR1_ORF133
ENSG00000231312_TR1_ORF133 <- data_lnc_MD |>
    dplyr::select("ENSG00000231312_TR1_ORF133",
                  condition, fiber_type) |>
    dplyr::mutate(condition = factor(condition, levels = c("control", "actin", "troponin"))) |>
    ggplot2::ggplot(ggplot2::aes(x = condition,
                                 y = ENSG00000231312_TR1_ORF133)) +
    ggplot2::geom_boxplot(ggplot2::aes(fill = condition),
                          alpha = 0.85) +
    ggplot2::geom_point(ggplot2::aes(x = condition, y = ENSG00000231312_TR1_ORF133, fill = condition), shape = 21, stroke = 0.5, color = "black") +
    # ggpubr::geom_bracket(xmin = 1, xmax = 2, y.position = 16.7, label = paste("Adj.P.Val = ", round(
    #     DE_control_v_actin |>
    #         dplyr::filter(grepl("ANXA6",Genes)) |>
    #         dplyr::pull(adj.P.Val),
    #     3
    # )),
    # label.size = 2) +
    # ggpubr::geom_bracket(xmin = 1, xmax = 3, y.position = 17.1, label = paste("Adj.P.Val = ", round(
    #     DE_control_v_troponin |>
    #         dplyr::filter(grepl("ANXA6", Genes)) |>
    #         dplyr::pull(adj.P.Val),
    #     3
    # )),
    # label.size = 2) +
    ggplot2::scale_fill_manual(values=c("#969594",
                                        "#E48C2AFF",
                                        "#d662c4")) +
    ggplot2::ggtitle("ENSG00000231312_TR1_ORF133") +
    ggplot2::labs(
        y = "LFQ intensities (log2)",
        x = ""
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
        legend.position = "none",
        axis.text.x = ggplot2::element_text(size = 6,
                                            vjust = 1,
                                            hjust=1),
        text = ggplot2::element_text(face="bold",
                                     colour="black",
                                     size=7),
        plot.title = ggplot2::element_text(face = "bold",
                                           size = 7,
                                           hjust = 0.5),
        plot.margin = grid::unit(c(1,1,1,1), "mm"),
        strip.text = ggplot2::element_text(size = 6),
        # axis.title.y = ggplot2::element_blank()
    ) +
    ggplot2::facet_grid(~fiber_type)

# ENSG00000232046_TR1_ORF437
ENSG00000232046_TR1_ORF437 <- data_lnc_MD |>
    dplyr::select("ENSG00000232046_TR1_ORF437",
                  condition, fiber_type) |>
    dplyr::mutate(condition = factor(condition, levels = c("control", "actin", "troponin"))) |>
    ggplot2::ggplot(ggplot2::aes(x = condition,
                                 y = ENSG00000232046_TR1_ORF437)) +
    ggplot2::geom_boxplot(ggplot2::aes(fill = condition),
                          alpha = 0.85) +
    ggplot2::geom_point(ggplot2::aes(x = condition, y = ENSG00000232046_TR1_ORF437, fill = condition), shape = 21, stroke = 0.5, color = "black") +
    # ggpubr::geom_bracket(xmin = 1, xmax = 2, y.position = 16.7, label = paste("Adj.P.Val = ", round(
    #     DE_control_v_actin |>
    #         dplyr::filter(grepl("ANXA6",Genes)) |>
    #         dplyr::pull(adj.P.Val),
    #     3
    # )),
    # label.size = 2) +
    # ggpubr::geom_bracket(xmin = 1, xmax = 3, y.position = 17.1, label = paste("Adj.P.Val = ", round(
    #     DE_control_v_troponin |>
    #         dplyr::filter(grepl("ANXA6", Genes)) |>
    #         dplyr::pull(adj.P.Val),
    #     3
    # )),
    # label.size = 2) +
    ggplot2::scale_fill_manual(values=c("#969594",
                                        "#E48C2AFF",
                                        "#d662c4")) +
    ggplot2::ggtitle("ENSG00000232046_TR1_ORF437") +
    ggplot2::labs(
        y = "LFQ intensities (log2)",
        x = ""
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
        legend.position = "none",
        axis.text.x = ggplot2::element_text(size = 6,
                                            vjust = 1,
                                            hjust=1),
        text = ggplot2::element_text(face="bold",
                                     colour="black",
                                     size=7),
        plot.title = ggplot2::element_text(face = "bold",
                                           size = 7,
                                           hjust = 0.5),
        plot.margin = grid::unit(c(1,1,1,1), "mm"),
        strip.text = ggplot2::element_text(size = 6),
        axis.title.y = ggplot2::element_blank()
    ) +
    ggplot2::facet_grid(~fiber_type)

# Join plots

patchwork::wrap_plots(
    ENSG00000177822_TR4_ORF9,
    ENSG00000215483_TR14_ORF67,
    ENSG00000229425_TR25_ORF40,
    ENSG00000231312_TR1_ORF133,
    ENSG00000232046_TR1_ORF437,
    nrow = 2
)

ggplot2::ggsave(here::here("doc/figures/figure_6_S3/figure_6_S3A.pdf"),
                units = "mm",
                width = 190,
                height = 90)

################################################################################################################################################
#################################################     Panel B  ##############################################################
################################################################################################################################################

data_MD_pseudobulk <- vroom::vroom(here::here("data/data_MD_pseudobulk.csv")) |>
    tibble::column_to_rownames("Genes")

data_grouping <- data.frame(
    "fiber_ID" = colnames(data_MD_pseudobulk),
    "condition" = c(
        "ACTA1-NM",
        "ACTA1-NM",
        "ACTA1-NM",
        "control",
        "control",
        "control",
        "TNNT1-NM",
        "TNNT1-NM",
        "TNNT1-NM"
    )
)

DE_control_v_actin <- vroom::vroom(here::here("data/DE_analysis_MD/actin_v_controls.csv"))

DE_control_v_troponin <- vroom::vroom(here::here("data/DE_analysis_MD/troponin_v_controls.csv"))

DE_actin_v_troponin <- vroom::vroom(here::here("data/DE_analysis_MD/actin_v_troponin.csv"))

metadata_MD <- vroom::vroom(here::here("data/metadata_MD_w_fiber_type_w_anonim.csv")) |>
    dplyr::mutate(condition = dplyr::case_when(
        condition == "actin" ~ "ACTA1-NM",
        condition == "troponin" ~ "TNNT1-NM",
        TRUE ~ "control"
    ))

muscle_disease_data <- vroom::vroom(
    here::here("data/data_muscle_disease.csv")
) |>
    dplyr::select(!...1) |>
    tibble::column_to_rownames("Gene_name") |>
    dplyr::mutate(dplyr::across(
        .cols = everything(),
        log2
    ))

data_wrangled <- muscle_disease_data

# ANXA1 -------------------------------------------------------------------

ANXA1 <- data_MD_pseudobulk |>
    tibble::rownames_to_column("Genes") |>
    dplyr::filter(Genes == "ANXA1") |>
    tibble::column_to_rownames("Genes") |>
    t() |>
    as.data.frame() |>
    tibble::rownames_to_column("sample_ID") |>
    dplyr::mutate(condition = c(
        rep("ACTA1-NM",3),
        rep("control", 3),
        rep("TNNT1-NM", 3)
    ))
# dplyr::inner_join(metadata_MD |>
#                       dplyr::select(fiber_ID, condition))

ANXA1$condition <- factor(ANXA1$condition,
                          levels = c("control",
                                     "ACTA1-NM",
                                     "TNNT1-NM"))

ANXA1_plot <- ggplot2::ggplot(ANXA1,
                              ggplot2::aes(x = condition,
                                           y = ANXA1)) +
    ggplot2::geom_boxplot(ggplot2::aes(fill = condition),
                          alpha = 0.85) +
    ggplot2::geom_point(ggplot2::aes(x = condition, y = ANXA1, fill = condition), shape = 21, stroke = 0.5, color = "black") +
    ggpubr::geom_bracket(xmin = 1, xmax = 2, y.position = 14.75, label = paste("Adj.P.Val = ", round(
        DE_control_v_actin |>
            dplyr::filter(grepl("ANXA1",names)) |>
            dplyr::pull(adj.P.Val),
        3
    )),
    label.size = 2) +
    ggpubr::geom_bracket(xmin = 1, xmax = 3, y.position = 15.5, label = paste("Adj.P.Val = ", round(
        DE_control_v_troponin |>
            dplyr::filter(grepl("ANXA1", Genes)) |>
            dplyr::pull(adj.P.Val),
        3
    )),
    label.size = 2) +
    ggplot2::scale_fill_manual(values=c("#969594",
                                        "#E48C2AFF",
                                        "#d662c4")) +
    ggplot2::ggtitle("ANXA1") +
    ggplot2::labs(
        y = "LFQ intensities (log2)",
        x = ""
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
        legend.position = "none",
        axis.text.x = ggplot2::element_text(size = 8,
                                            angle = 45,
                                            vjust = 1,
                                            hjust=1),
        text = ggplot2::element_text(face="bold",
                                     colour="black",
                                     size=7),
        plot.title = ggplot2::element_text(face = 4,
                                           size = 10,
                                           hjust = 0.5),
        plot.margin = grid::unit(c(0,0,0,0), "mm")
    ) +
    ggplot2::ylim(11, 16)

# ggplot2::ggsave(here::here("doc/figures/figure_6/annexins_expression.png"),
#                 units = "mm",
#                 height = 40,
#                 width = 40)


# ANXA2 -------------------------------------------------------------------

ANXA2 <- data_MD_pseudobulk |>
    tibble::rownames_to_column("Genes") |>
    dplyr::filter(Genes == "ANXA2") |>
    tibble::column_to_rownames("Genes") |>
    t() |>
    as.data.frame() |>
    tibble::rownames_to_column("sample_ID") |>
    dplyr::mutate(condition = c(
        rep("ACTA1-NM",3),
        rep("control", 3),
        rep("TNNT1-NM", 3)
    ))
# dplyr::inner_join(metadata_MD |>
#                       dplyr::select(fiber_ID, condition))

ANXA2$condition <- factor(ANXA2$condition,
                          levels = c("control",
                                     "ACTA1-NM",
                                     "TNNT1-NM"))

ANXA2_plot <- ggplot2::ggplot(ANXA2,
                              ggplot2::aes(x = condition,
                                           y = ANXA2)) +
    ggplot2::geom_boxplot(ggplot2::aes(fill = condition),
                          alpha = 0.85) +
    ggplot2::geom_point(ggplot2::aes(x = condition, y = ANXA2, fill = condition), shape = 21, stroke = 0.5, color = "black") +
    ggpubr::geom_bracket(xmin = 1, xmax = 2, y.position = 16.25, label = paste("Adj.P.Val = ", round(
        DE_control_v_actin |>
            dplyr::filter(grepl("ANXA2",names)) |>
            dplyr::pull(adj.P.Val),
        3
    )),
    label.size = 2) +
    ggpubr::geom_bracket(xmin = 1, xmax = 3, y.position = 16.8, label = paste("Adj.P.Val = ", round(
        DE_control_v_troponin |>
            dplyr::filter(grepl("ANXA2", Genes)) |>
            dplyr::pull(adj.P.Val),
        3
    )),
    label.size = 2) +
    ggplot2::scale_fill_manual(values=c("#969594",
                                        "#E48C2AFF",
                                        "#d662c4")) +
    ggplot2::ggtitle("ANXA2") +
    ggplot2::labs(
        y = "LFQ intensities (log2)",
        x = ""
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
        legend.position = "none",
        axis.text.x = ggplot2::element_text(size = 8,
                                            angle = 45,
                                            vjust = 1,
                                            hjust=1),
        text = ggplot2::element_text(face="bold",
                                     colour="black",
                                     size=7),
        plot.title = ggplot2::element_text(face = 4,
                                           size = 10,
                                           hjust = 0.5),
        plot.margin = grid::unit(c(0,0,0,0), "mm"),
        axis.title.y = ggplot2::element_blank()
    ) +
    ggplot2::ylim(14.25, 17)


# ANXA5 -------------------------------------------------------------------

ANXA5 <- data_MD_pseudobulk |>
    tibble::rownames_to_column("Genes") |>
    dplyr::filter(Genes == "ANXA5") |>
    tibble::column_to_rownames("Genes") |>
    t() |>
    as.data.frame() |>
    tibble::rownames_to_column("sample_ID") |>
    dplyr::mutate(condition = c(
        rep("ACTA1-NM",3),
        rep("control", 3),
        rep("TNNT1-NM", 3)
    ))
# dplyr::inner_join(metadata_MD |>
#                       dplyr::select(fiber_ID, condition))

ANXA5$condition <- factor(ANXA5$condition,
                          levels = c("control",
                                     "ACTA1-NM",
                                     "TNNT1-NM"))

ANXA5_plot <- ggplot2::ggplot(ANXA5,
                              ggplot2::aes(x = condition,
                                           y = ANXA5)) +
    ggplot2::geom_boxplot(ggplot2::aes(fill = condition),
                          alpha = 0.85) +
    ggplot2::geom_point(ggplot2::aes(x = condition, y = ANXA5, fill = condition), shape = 21, stroke = 0.5, color = "black") +
    ggpubr::geom_bracket(xmin = 1, xmax = 2, y.position = 16.75, label = paste("Adj.P.Val = ", round(
        DE_control_v_actin |>
            dplyr::filter(grepl("ANXA5",names)) |>
            dplyr::pull(adj.P.Val),
        3
    )),
    label.size = 2) +
    ggpubr::geom_bracket(xmin = 1, xmax = 3, y.position = 17.35, label = paste("Adj.P.Val = ", round(
        DE_control_v_troponin |>
            dplyr::filter(grepl("ANXA5", Genes)) |>
            dplyr::pull(adj.P.Val),
        3
    )),
    label.size = 2) +
    ggplot2::scale_fill_manual(values=c("#969594",
                                        "#E48C2AFF",
                                        "#d662c4")) +
    ggplot2::ggtitle("ANXA5") +
    ggplot2::labs(
        y = "LFQ intensities (log2)",
        x = ""
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
        legend.position = "none",
        axis.text.x = ggplot2::element_text(size = 8,
                                            angle = 45,
                                            vjust = 1,
                                            hjust=1),
        text = ggplot2::element_text(face="bold",
                                     colour="black",
                                     size=7),
        plot.title = ggplot2::element_text(face = 4,
                                           size = 10,
                                           hjust = 0.5),
        plot.margin = grid::unit(c(0,0,0,0), "mm"),
        axis.title.y = ggplot2::element_blank()
    ) +
    ggplot2::ylim(14, 17.65)


# ANXA6 -------------------------------------------------------------------

ANXA6 <- data_MD_pseudobulk |>
    tibble::rownames_to_column("Genes") |>
    dplyr::filter(Genes == "ANXA6") |>
    tibble::column_to_rownames("Genes") |>
    t() |>
    as.data.frame() |>
    tibble::rownames_to_column("sample_ID") |>
    dplyr::mutate(condition = c(
        rep("ACTA1-NM",3),
        rep("control", 3),
        rep("TNNT1-NM", 3)
    ))
# dplyr::inner_join(metadata_MD |>
#                       dplyr::select(fiber_ID, condition))

ANXA6$condition <- factor(ANXA6$condition,
                          levels = c("control",
                                     "ACTA1-NM",
                                     "TNNT1-NM"))

ANXA6_plot <- ggplot2::ggplot(ANXA6,
                              ggplot2::aes(x = condition,
                                           y = ANXA6)) +
    ggplot2::geom_boxplot(ggplot2::aes(fill = condition),
                          alpha = 0.85) +
    ggplot2::geom_point(ggplot2::aes(x = condition, y = ANXA6, fill = condition), shape = 21, stroke = 0.5, color = "black") +
    ggpubr::geom_bracket(xmin = 1, xmax = 2, y.position = 16.7, label = paste("Adj.P.Val = ", round(
        DE_control_v_actin |>
            dplyr::filter(grepl("ANXA6",Genes)) |>
            dplyr::pull(adj.P.Val),
        3
    )),
    label.size = 2) +
    ggpubr::geom_bracket(xmin = 1, xmax = 3, y.position = 17.1, label = paste("Adj.P.Val = ", round(
        DE_control_v_troponin |>
            dplyr::filter(grepl("ANXA6", Genes)) |>
            dplyr::pull(adj.P.Val),
        3
    )),
    label.size = 2) +
    ggplot2::scale_fill_manual(values=c("#969594",
                                        "#E48C2AFF",
                                        "#d662c4")) +
    ggplot2::ggtitle("ANXA6") +
    ggplot2::labs(
        y = "LFQ intensities (log2)",
        x = ""
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
        legend.position = "none",
        axis.text.x = ggplot2::element_text(size = 8,
                                            angle = 45,
                                            vjust = 1,
                                            hjust=1),
        text = ggplot2::element_text(face="bold",
                                     colour="black",
                                     size=7),
        plot.title = ggplot2::element_text(face = 4,
                                           size = 10,
                                           hjust = 0.5),
        plot.margin = grid::unit(c(0,0,0,0), "mm"),
        axis.title.y = ggplot2::element_blank()
    ) +
    ggplot2::ylim(14.8, 17.5)



# Annexin plot ------------------------------------------------------------

annexin_plot <- patchwork::wrap_plots(ANXA1_plot,
                                      ANXA2_plot,
                                      ANXA5_plot,
                                      ANXA6_plot,
                                      nrow = 1, ncol = 4) +
    patchwork::plot_annotation(
        title = "Annexin expression",
        theme = ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size=12, face="bold"))
    )

ggplot2::ggsave(plot = annexin_plot,
                here::here("doc/figures/figure_6_S3/figure_6_S3B.pdf"),
                units = "mm",
                width = 120,
                height = 60)
