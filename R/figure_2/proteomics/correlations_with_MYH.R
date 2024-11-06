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


# Load data ---------------------------------------------------------------

data_proteomics <-vroom::vroom(here::here("data/data_proteomics_filtered.csv")) |>
    dplyr::select(!1) |>
    tibble::column_to_rownames("Gene.name") |>
    filtering_rows_Na(percentage_accepted_missing = 0.3) |>
    log2() |>
    as.data.frame() |>
    tibble::rownames_to_column("Genes") |>
    dplyr::mutate(
        Genes = gsub(
            pattern = "-",
            replacement = ".",
            Genes
        )
    )

seurat_clusters <-
    vroom::vroom(here::here("data/proteomics clustering/seurat_clusters_6PC_res04.csv"))

metadata <- vroom::vroom(here::here("data/metadata_proteomics.csv"))|>
    dplyr::rename("fiberID" = 1) |>
    dplyr::inner_join(seurat_clusters) |>
    dplyr::mutate(fiber_type_seurat = dplyr::case_when(
        seurat_clusters == "3" ~ "slow",
        seurat_clusters == "1" ~ "slow",
        TRUE ~ "fast"
    )) |>
    dplyr::mutate("dataset" = dplyr::case_when(
        digestion_batch == 1 ~ "before_Xmas",
        digestion_batch == 2 ~ "before_Xmas",
        digestion_batch == 3 ~ "after_Xmas",
        digestion_batch == 4 ~ "after_Xmas",
        digestion_batch == 5 ~ "after_Xmas",
        digestion_batch == 6 ~ "after_Xmas"
    )) |>
    tibble::column_to_rownames("fiberID")

list_proteomics <- data_proteomics |>
    dplyr::group_split(
        Genes
    )

names(list_proteomics) <- data_proteomics |>
    dplyr::arrange(Genes) |>
    dplyr::pull(Genes)

gene_names <- data_proteomics |>
    dplyr::arrange(Genes) |>
    dplyr::pull(Genes)

# Create correlation matrix -----------------------------------------------

correlation_function <- function(.list, .gene, .myh_to_correlate) {

x <- .list[[.myh_to_correlate]] |>
    tibble::column_to_rownames("Genes") |>
    t() |>
    as.data.frame() |>
    tibble::rownames_to_column("Genes") |>
    dplyr::pull(2, name = Genes)

y <- .list[[.gene]] |>
    tibble::column_to_rownames("Genes") |>
    t() |>
    as.data.frame() |>
    tibble::rownames_to_column("Genes") |>
    dplyr::pull(2, name = Genes)

test <- cor.test(x,
                 y,
                 na.rm = TRUE,
                 method = "pearson")

result <- data.frame(gene = .list[[.gene]]$Genes,
                     statistic = test$statistic,
                     correlation = test$estimate,
                     p_val = test$p.value) |>
    tibble::remove_rownames()

return(result)
}

num_cores <- detectCores() - 1  # Adjust the number of cores as needed
registerDoParallel(cores = num_cores)

correlation_matrix_MYH7 <- purrr::map(
    gene_names,
    ~ correlation_function(
        .list = list_proteomics,
        .myh_to_correlate = "MYH7",
        .gene = .x
    )
) |>
    purrr::list_rbind() |>
    dplyr::mutate(fdr = p.adjust(p_val, method = "BH"))

correlation_matrix_MYH2 <- purrr::map(
    gene_names,
    ~ correlation_function(
        .list = list_proteomics,
        .myh_to_correlate = "MYH2",
        .gene = .x
    )
) |>
    purrr::list_rbind() |>
    dplyr::mutate(fdr = p.adjust(p_val, method = "BH"))

correlation_matrix_MYH1 <- purrr::map(
    gene_names,
    ~ correlation_function(
        .list = list_proteomics,
        .myh_to_correlate = "MYH1",
        .gene = .x
    )
) |>
    purrr::list_rbind() |>
    dplyr::mutate(fdr = p.adjust(p_val, method = "BH"))

# Stop parallel processing
stopImplicitCluster()

labels_MYH7 <- c(correlation_matrix_MYH7 |>
    dplyr::arrange(desc(correlation)) |>
    dplyr::slice_head(n = 10) |>
    dplyr::pull(gene),
    correlation_matrix_MYH7 |>
        dplyr::arrange(correlation) |>
        dplyr::slice_head(n = 10) |>
        dplyr::pull(gene)
)
correlation_matrix_MYH7 |>
    dplyr::mutate(
        color = dplyr::case_when(
            fdr < 0.05 & correlation > 0.7 ~ "positive",
            fdr < 0.05 & correlation < -0.7 ~ "negative",
            TRUE ~ "blank"
        )
    ) |>
    dplyr::mutate(
        color = factor(color, levels = c("positive", "blank", "negative"))
    ) |>
    ggplot2::ggplot(
        ggplot2::aes(
            x = correlation,
            y = -log10(p_val),
            color = color
        )
    ) +
    ggplot2::theme_classic() +
    ggplot2::geom_point(size = 1, alpha = 0.65) +
    ggrepel::geom_label_repel(
        data = correlation_matrix_MYH7 |>
            dplyr::mutate(
            label = dplyr::case_when(
                gene %in% labels_MYH7 ~ gene,
            )) |>
            dplyr::mutate(
                color = dplyr::case_when(
                    fdr < 0.05 & correlation > 0.7 ~ "positive",
                    fdr < 0.05 & correlation < -0.7 ~ "negative",
                    TRUE ~ "blank"
                )
            ) |>
            dplyr::mutate(
                color = factor(color, levels = c("positive", "blank", "negative"))
            ),
        mapping = ggplot2::aes(
            x = correlation,
            y = -log10(p_val),
            fill = color,
            label = label
        ),
        max.overlaps = Inf,
        size = 2,
        color = "black",
        label.padding = 0.1,
        min.segment.length = 0.1,
        segment.size = 0.2,
        force = 20
    ) +
    ggplot2::ggtitle("Proteins correlating with MYH7") +
    ggplot2::ylim(c(0, 350)) +
    ggplot2::scale_color_manual(
        values = c("#990000",
                   "lightgray",
                   "#0570b0")
    ) +
    ggplot2::scale_fill_manual(values = c(
        "#ebcccc",
        "lightgrey",
        "#cde2ef"
    )) +
    ggplot2::xlab("Pearson's r") +
    ggplot2::ylab("-log10(P-value)") +
    ggplot2::theme(
        legend.position = "none",
        text = ggplot2::element_text(size = 6),
        plot.title = ggplot2::element_text(size = 8, face = "bold", hjust = 0.5)
    )

ggplot2::ggsave(here::here("doc/figures/figure_4/proteomics_correlation_MYH7.pdf"),
                units = "mm",
                height = 60,
                width = 60)


# MYH2 --------------------------------------------------------------------

labels_MYH2 <- c(correlation_matrix_MYH2 |>
                     dplyr::arrange(desc(correlation)) |>
                     dplyr::slice_head(n = 10) |>
                     dplyr::pull(gene),
                 correlation_matrix_MYH2 |>
                     dplyr::arrange(correlation) |>
                     dplyr::slice_head(n = 10) |>
                     dplyr::pull(gene))

correlation_matrix_MYH2 |>
    dplyr::mutate(
        color = dplyr::case_when(
            fdr < 0.05 & correlation > 0.7 ~ "positive",
            fdr < 0.05 & correlation < -0.7 ~ "negative",
            TRUE ~ "blank"
        )
    ) |>
    dplyr::mutate(
        color = factor(color, levels = c("positive", "blank", "negative"))
    ) |>
    ggplot2::ggplot(
        ggplot2::aes(
            x = correlation,
            y = -log10(p_val),
            color = color
        )
    ) +
    ggplot2::theme_classic() +
    ggplot2::geom_point(size = 1, alpha = 0.65) +
    ggrepel::geom_label_repel(
        data = correlation_matrix_MYH2 |>
            dplyr::mutate(
                label = dplyr::case_when(
                    gene %in% labels_MYH2 ~ gene,
                )) |>
            dplyr::mutate(
                color = dplyr::case_when(
                    fdr < 0.05 & correlation > 0.7 ~ "positive",
                    fdr < 0.05 & correlation < -0.7 ~ "negative",
                    TRUE ~ "blank"
                )
            ) |>
            dplyr::mutate(
                color = factor(color, levels = c("positive", "blank", "negative"))
            ),
        mapping = ggplot2::aes(
            x = correlation,
            y = -log10(p_val),
            fill = color,
            label = label
        ),
        max.overlaps = Inf,
        size = 2,
        color = "black",
        label.padding = 0.1,
        min.segment.length = 0.1,
        segment.size = 0.2,
        force = 20
    ) +
    ggplot2::ggtitle("Proteins correlating with MYH2") +
    ggplot2::ylim(c(0, 350)) +
    ggplot2::scale_color_manual(
        values = c("#990000",
                   "lightgray",
                   "#0570b0")
    ) +
    ggplot2::scale_fill_manual(values = c(
        "#ebcccc",
        "lightgrey",
        "#cde2ef"
    )) +
    ggplot2::xlab("Pearson's r") +
    ggplot2::ylab("-log10(P-value)") +
    ggplot2::theme(
        legend.position = "none",
        text = ggplot2::element_text(size = 6),
        plot.title = ggplot2::element_text(size = 8, face = "bold", hjust = 0.5)
    )

ggplot2::ggsave(here::here("doc/figures/figure_4/proteomics_correlation_MYH2.pdf"),
                units = "mm",
                height = 60,
                width = 60)

# MYH1 --------------------------------------------------------------------

correlation_matrix_MYH1 |>
    dplyr::mutate(
        color = dplyr::case_when(
            fdr < 0.05 & correlation > 0.7 ~ "positive",
            fdr < 0.05 & correlation < -0.7 ~ "negative",
            TRUE ~ "blank"
        )
    ) |>
    dplyr::mutate(
        color = factor(color, levels = c("positive", "blank", "negative"))
    ) |>
    ggplot2::ggplot(
        ggplot2::aes(
            x = correlation,
            y = -log10(p_val),
            color = color
        )
    ) +
    ggplot2::theme_classic() +
    ggplot2::geom_point(size = 1, alpha = 0.65) +
    ggrepel::geom_label_repel(
        data = correlation_matrix_MYH1 |>
            dplyr::mutate(
                label = dplyr::case_when(
                    fdr < 0.05 & correlation > 0.7 ~ gene,
                    fdr < 0.05 & correlation < -0.7 ~ gene,
                    TRUE ~ ""
                )) |>
            dplyr::mutate(
                color = dplyr::case_when(
                    fdr < 0.05 & correlation > 0.7 ~ "positive",
                    fdr < 0.05 & correlation < -0.7 ~ "negative",
                    TRUE ~ "blank"
                )
            ) |>
            dplyr::mutate(
                color = factor(color, levels = c("positive", "blank", "negative"))
            ),
        mapping = ggplot2::aes(
            x = correlation,
            y = -log10(p_val),
            fill = color,
            label = label
        ),
        max.overlaps = Inf,
        size = 2,
        color = "black",
        label.padding = 0.1,
        min.segment.length = 0.1,
        segment.size = 0.2,
        force = 20
    ) +
    ggplot2::ggtitle("Proteins correlating with MYH1") +
    ggplot2::ylim(c(0, 350)) +
    ggplot2::scale_color_manual(
        values = c("#990000",
                   "lightgray",
                   "#0570b0")
    ) +
    ggplot2::scale_fill_manual(values = c(
        "#ebcccc",
        "lightgrey",
        "#cde2ef"
    )) +
    ggplot2::xlab("Pearson's r") +
    ggplot2::ylab("-log10(P-value)") +
    ggplot2::theme(
        legend.position = "none",
        text = ggplot2::element_text(size = 6),
        plot.title = ggplot2::element_text(size = 8, face = "bold", hjust = 0.5)
    )

ggplot2::ggsave(here::here("doc/figures/figure_4/proteomics_correlation_MYH1.pdf"),
                units = "mm",
                height = 60,
                width = 60)

# Linear mixed model version ----------------------------------------------

library(scales)
library(lmerTest)
library(doParallel)
library(foreach)


# Define the mixed.anova function with covariates
mixed.anova <- function(dat, expo, outc, formel, covariates = NULL) {
    require("lmerTest")
    # Create a string of covariates if they are provided
    covariate_str <- if (!is.null(covariates)) {
        paste(covariates, collapse = " + ")
    } else {
        ""
    }
    # Run across all exposures and all outcomes of interest
    res <- lapply(outc, function(i) {
        # Loop over all exposures of interest
        tmp <- lapply(expo, function(j) {
            # Create formula
            ff <- if (covariate_str != "") {
                paste0(i, "~", j, " + ", covariate_str, formel)
            } else {
                paste0(i, "~", j, formel)
            }
            print(ff)
            # Run the model
            ml <- lmer(as.formula(ff), data=dat)
            ff <- summary(ml)$coefficients
            # Return needed information
            ff <- data.frame(outcome=i, exposure=j, ind=row.names(ff), beta=ff[,1], se=ff[,2], pval=ff[,5])
            # Reshape the data
            ff <- reshape(ff, idvar = c("outcome", "exposure"), timevar = "ind", direction = "wide")
            # Exchange names to allow combination afterwards
            names(ff) <- gsub(j, "exposure", names(ff))
            # Add number of observations
            ff$n <- nrow(na.omit(dat[, c(expo, outc)]))
            # Add ANOVA p-values
            ml <- anova(ml)
            ff[, paste0("pval.aov.", rownames(ml))] <- ml$`Pr(>F)`
            return(ff)
        })
        # Combine into one data frame
        tmp <- do.call(rbind, tmp)
        return(tmp)
    })
    # Combine into one data frame
    res <- do.call(rbind, res)
    print(head(res))
    # Get only information needed
    res <- res[, c("outcome", "n", grep("exposure|aov", names(res), value=T))]
    return(res)
}
# Define a function to perform mixed ANOVA in parallel
mixed_anova_parallel <- function(dat, expo, outc, formel, covariates = NULL) {
    res <- foreach(i = 1:length(outc), .combine = rbind, .packages = c("lmerTest", "dplyr", "reshape2"), .export = c("mixed.anova")) %dopar% {
        mixed.anova(dat, expo, outc[i], formel, covariates)  # Using mixed.anova function
    }
    return(res)
}

# Ensure MYH7 is numeric
if (!is.numeric(data_lmodel$MYH7)) {
    data_lmodel$MYH7 <- as.numeric(data_lmodel$MYH7)
}

# Check that MYH7 is numeric
if (!is.numeric(data_lmodel$MYH7)) {
    stop("MYH7 column is not numeric in fast.dat")
}

data_lmodel <- data_proteomics |>
    tibble::column_to_rownames("Genes") |>
    t() |>
    as.data.frame() |>
    tibble::rownames_to_column("sample_ID") |>
    dplyr::inner_join(
        metadata |>
            tibble::rownames_to_column("sample_ID") |>
            dplyr::select(c(sample_ID, subject, fiber_type_seurat))
    )

# Identify the other proteins
other_proteins <- setdiff(unique(data_proteomics$Genes), "MYH2")

# Set up parallel processing
num_cores <- detectCores() - 1  # Adjust the number of cores as needed
registerDoParallel(cores = num_cores)

# Call the function with MYH7 as the exposure
res.model.linear <- mixed_anova_parallel(
    data_lmodel,
    "MYH2",  # Set MYH7 as the exposure
    other_proteins,  # Other proteins as outcomes
    "+ (1|subject)",
    covariates = "fiber_type_seurat"
)
# Display the results
head(res.model.linear)

# Stop parallel processing
stopImplicitCluster()

# Join additional information from prot.label
res.model.linear <- res.model.linear |>
    dplyr::mutate(fdr.aov = p.adjust(pval.aov.MYH2, method = "BH"))

plotly::ggplotly(
res.model.linear |>
    ggplot2::ggplot(
        ggplot2::aes(
            x = beta.exposure,
            y = -log10(pval.aov.MYH2),
            color = fdr.aov < 0.05,
            names = outcome
        )
    ) +
    ggplot2::theme_bw() +
    ggplot2::geom_point()
)


