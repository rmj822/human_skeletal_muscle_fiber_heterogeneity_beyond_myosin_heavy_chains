

# muscle disease ----------------------------------------------------------

# Load raw data -----------------------------------------------------------

data_raw <- vroom::vroom(here::here("data-raw/muscle_disease_Report_ProteinQuant_Diana (Pivot).tsv"))

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

Gene.name <- data_wrangled$PG.Genes

selection_vector <- as.character(metadata$random_order)

data_wrangled <- data_wrangled |>
    dplyr::select(PG.ProteinGroups,
                  selection_vector) |>
    tibble::column_to_rownames("PG.ProteinGroups") |>
    t() |>
    as.data.frame() |>
    tibble::rownames_to_column("random_order") |>
    dplyr::mutate(randomized_number = as.numeric(random_order))


data_wrangled <- data_wrangled |>
    dplyr::inner_join(
        metadata |>
            dplyr::mutate(random_order = as.character(random_order)) |>
            dplyr::select(fiberID,
                          random_order)
    ) |>
    tibble::column_to_rownames("fiberID") |>
    dplyr::select(!randomized_number) |>
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
                     here::here("data/lnc_quant/data_row_filtered_muscle_disease.csv"))


# 1000 fiber proteome -----------------------------------------------------

data_raw <- vroom::vroom(here::here("data-raw/20230930_112903_SMF_directDIA_lncRNAs_DS_13092023_Report_proteins.tsv"))

metadata <- vroom::vroom(here::here("data/metadata_proteomics.csv")) |>
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

Gene.name <- data_wrangled$PG.Genes

selection_vector <- as.character(metadata$randomized_number)

data_wrangled <- data_wrangled |>
    dplyr::select(PG.ProteinGroups,
                  selection_vector) |>
    tibble::column_to_rownames("PG.ProteinGroups") |>
    t() |>
    as.data.frame() |>
    tibble::rownames_to_column("randomized_number") |>
    dplyr::mutate(randomized_number = as.numeric(randomized_number))


data_wrangled <- data_wrangled |>
    dplyr::inner_join(
        metadata |>
            # dplyr::mutate(randomized_number = as.numeric(randomized_number)) |>
            dplyr::select(fiberID,
                          randomized_number)
    ) |>
    tibble::column_to_rownames("fiberID") |>
    dplyr::select(!randomized_number) |>
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

data <- data_wrangled |>
    log2() |>
    as.data.frame() |>
    tibble::rownames_to_column("protein_groups") |>
    dplyr::filter(grepl("ENS",
                        protein_groups))

write.csv(data_row_filtered,
          here::here("data/lnc_quant/data_row_filtered_1000_fiber.csv"))

write.csv(data,
          here::here("doc/supplementary_tables/supp_tables_nature_com/table_s16.csv"))
