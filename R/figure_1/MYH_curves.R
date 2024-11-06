# Load your data.
# Im using a format where samples are columns and genes are rows.
# Im sending my MYH-filtered data as an example:

proteomics_MYH_filtered <- vroom::vroom(
  here::here("data-raw/raw_proteomics.csv")
) |>
  as.data.frame() |>
    dplyr::filter(!duplicated(Gene_name)) |>
    dplyr::filter(!is.na(Gene_name)) |>
  tibble::column_to_rownames("Gene_name")

# Sorting dataset with the genes that are used for fiber typing,
# plotting intensities and rank:
data_fiber_typing <- proteomics_MYH_filtered |>
  t() |>
  as.data.frame() |>
  dplyr::select("MYH7", "MYH2", "MYH1") |>
  dplyr::mutate("Fiber_ID" = colnames(proteomics_MYH_filtered)) |>
  dplyr::mutate(across(.cols = c("MYH7", "MYH2", "MYH1"), as.numeric)) |>
  dplyr::arrange(desc(MYH7)) |>
  dplyr::mutate("order_7" = seq_len(ncol(proteomics_MYH_filtered))) |>
  dplyr::arrange(desc(MYH2)) |>
  dplyr::mutate("order_2" = seq_len(ncol(proteomics_MYH_filtered)))

data_fiber_typing |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = order_7,
      y = MYH7
    )
  ) +
  ggplot2::geom_point() +
  ggplot2::theme_minimal()


# Ben and I double checked Marta Murgia's paper on SMF proteomics
# and saw that her curve figure is made using percentages
# of MYH intensities, instead of raw intensities.
# Here, we will do the same, calculate percentage
# of MYH composition and rank fibers before plotting them:

sum_of_intensities <- proteomics_MYH_filtered |>
  t() |>
  as.data.frame() |>
  dplyr::select("MYH7", "MYH2", "MYH1") |>
  dplyr::mutate(across(.cols = c("MYH7", "MYH2", "MYH1"), as.numeric)) |>
  tidyr::replace_na(list(
    "MYH7" = 0,
    "MYH2" = 0,
    "MYH1" = 0
  )) |>
  rowSums()


perc_MYHs <- proteomics_MYH_filtered |>
  t() |>
  as.data.frame() |>
  dplyr::select("MYH7", "MYH2", "MYH1") |>
  dplyr::mutate(across(.cols = c("MYH7", "MYH2", "MYH1"), as.numeric)) |>
  tidyr::replace_na(list(
    "MYH7" = 0,
    "MYH2" = 0,
    "MYH1" = 0
  )) |>
  tibble::add_column(sum_of_intensities) |>
  dplyr::mutate(across(
    .cols = c("MYH7", "MYH2", "MYH1"),
    ~ .x / sum_of_intensities * 100
  )) |>
  as.data.frame() |>
  dplyr::arrange(desc(MYH7)) |>
  dplyr::mutate("order_7" = seq_len(ncol(proteomics_MYH_filtered))) |>
  dplyr::arrange(desc(MYH2)) |>
  dplyr::mutate("order_2" = seq_len(ncol(proteomics_MYH_filtered))) |>
  dplyr::arrange(desc(MYH1)) |>
  dplyr::mutate("order_1" = seq_len(ncol(proteomics_MYH_filtered))) |>
  tibble::rownames_to_column("fiber_ID") |>
  dplyr::arrange(desc(fiber_ID)) |>
  dplyr::select(c(
    "MYH7",
    "MYH2",
    "MYH1",
    "order_7",
    "order_2",
    "order_1",
    "fiber_ID"
  )) |>
  tidyr::pivot_longer(
    cols = c("MYH7", "MYH2", "MYH1"),
    names_to = "MYHs",
    values_to = "values"
  )

write.csv(perc_MYHs,
          here::here("data/perc_MYH_proteomics.csv"))

# Individual curve for intensity of MYH7:

perc_MYHs |>
  dplyr::filter(MYHs == "MYH7") |>
  ggplot2::ggplot(
    ggplot2::aes(order_7,
      values,
      color = MYHs
    )
  ) +
  ggplot2::geom_point() +
  ggplot2::scale_color_manual(values = "#440154FF") +
  ggplot2::theme_minimal()

# MYH7, MYH2 and MYH1 curves all together ranked by MYH7:

perc_MYHs |>
  ggplot2::ggplot(ggplot2::aes(
    x = order_7,
    y = values
  )) +
  ggplot2::geom_point(ggplot2::aes(colour = MYHs)) +
  ggplot2::scale_color_manual(values = c(
    "#fdc325",
    "#5DC863FF",
    "#440154FF"
  )) +
  ggplot2::theme_minimal()


# Finding bottom knee of MYH 7 curve,
# that's the one we are using as a threshold.
# The function detects top knee so we do the inverse of the curve:

MYH_7_curve <- perc_MYHs |>
  tidyr::pivot_wider(
    names_from = MYHs,
    values_from = values
  ) |>
  tibble::column_to_rownames("fiber_ID") |>
  dplyr::select(MYH7) |>
  dplyr::mutate(across(everything(), ~ 100 - .x)) |>
  t()

MYH_7_curve <- DropletUtils::barcodeRanks(MYH_7_curve, lower = 10)

bottom_knee_MYH7 <- 100 - MYH_7_curve@metadata$knee

perc_MYHs |>
  dplyr::filter(MYHs == "MYH7") |>
  ggplot2::ggplot(
    ggplot2::aes(order_7,
      values,
      color = MYHs
    )
  ) +
  ggplot2::geom_point() +
  ggplot2::scale_color_manual(values = "#440154FF") +
  ggplot2::geom_hline(yintercept = bottom_knee_MYH7, color = "red") +
  ggplot2::theme_minimal()

# Same with MYH2:

MYH_2_curve <- perc_MYHs |>
  tidyr::pivot_wider(
    names_from = MYHs,
    values_from = values
  ) |>
  tibble::column_to_rownames("fiber_ID") |>
  dplyr::select(MYH2) |>
  dplyr::mutate(across(everything(), ~ 100 - .x)) |>
  t()

MYH_2_curve <- DropletUtils::barcodeRanks(MYH_2_curve, lower = 10)

bottom_knee_MYH2 <- 100 - MYH_2_curve@metadata$knee

perc_MYHs |>
  dplyr::filter(MYHs == "MYH2") |>
  ggplot2::ggplot(
    ggplot2::aes(order_2,
      values,
      color = MYHs
    )
  ) +
  ggplot2::geom_point() +
  ggplot2::scale_color_manual(values = "#5DC863FF") +
  ggplot2::geom_hline(yintercept = bottom_knee_MYH2, color = "red") +
  ggplot2::theme_minimal()

# Now MYH1

MYH_1_curve <- perc_MYHs |>
  tidyr::pivot_wider(
    names_from = MYHs,
    values_from = values
  ) |>
  tibble::column_to_rownames("fiber_ID") |>
  dplyr::select(MYH1) |>
  dplyr::mutate(across(everything(), ~ 100 - .x)) |>
  t()

MYH_1_curve <- DropletUtils::barcodeRanks(MYH_1_curve, lower = 10)

bottom_knee_MYH1 <- 100 - MYH_1_curve@metadata$knee

perc_MYHs |>
  dplyr::filter(MYHs == "MYH1") |>
  ggplot2::ggplot(
    ggplot2::aes(order_1,
      values,
      color = MYHs
    )
  ) +
  ggplot2::geom_point() +
  ggplot2::scale_color_manual(values = "#fdc325") +
  ggplot2::geom_hline(yintercept = bottom_knee_MYH1, color = "red") +
  ggplot2::theme_minimal()

# This is the case when conditional that I have used to assign fiber types:

data_fiber_type <- perc_MYHs |>
  tidyr::pivot_wider(
    names_from = MYHs,
    values_from = values
  ) |>
  dplyr::mutate(
    fiber_type = dplyr::case_when(
      MYH7 >= bottom_knee_MYH7 &
        MYH2 <= bottom_knee_MYH2 ~ "Type 1",
      MYH7 >= bottom_knee_MYH7 &
        MYH2 >= bottom_knee_MYH2 ~ "Hybrid 1/2A",
      MYH1 >= bottom_knee_MYH1 &
        MYH2 >= bottom_knee_MYH2 ~ "Hybrid 2A/2X",
      MYH2 >= bottom_knee_MYH2 &
        MYH7 <= bottom_knee_MYH7 ~ "Type 2A",
      TRUE ~ "Hybrid 1/2A"
    )
  ) |>
  dplyr::arrange(desc(fiber_ID))

# Again, how the % of MYH look like themselves

data_fiber_type |>
  tidyr::pivot_longer(
    cols = c("MYH7", "MYH2", "MYH1"),
    names_to = "MYHs",
    values_to = "values"
  ) |>
  ggplot2::ggplot(ggplot2::aes(
    x = order_7,
    y = values
  )) +
  ggplot2::geom_point(ggplot2::aes(colour = MYHs)) +
  ggplot2::scale_color_manual(values = c(
    "#fdc325",
    "#5DC863FF",
    "#440154FF"
  )) +
  ggplot2::theme_minimal()

# When I colour by fiber type:

data_fiber_type |>
  tidyr::pivot_longer(
    cols = c("MYH7", "MYH2", "MYH1"),
    names_to = "MYHs",
    values_to = "values"
  ) |>
  ggplot2::ggplot(ggplot2::aes(
    x = order_7,
    y = values
  )) +
  ggplot2::geom_point(ggplot2::aes(colour = fiber_type)) +
  ggplot2::scale_color_manual(values = c(
    "#3B528BFF",
    "#fdc325",
    "#440154FF",
    "#5DC863FF"
  )) +
  ggplot2::theme_minimal()

# Overall number of fibers in each fiber type:

data_fiber_type |>
  ggplot2::ggplot(ggplot2::aes(
    x = fiber_type,
    fill = fiber_type
  )) +
  ggplot2::geom_bar(na.rm = TRUE) +
  ggplot2::scale_fill_manual("Fiber types",
    values = c(
      "#3B528BFF",
      "#fdc325",
      "#440154FF",
      "#5DC863FF"
    )
  ) +
  ggplot2::theme_light() +
  ggplot2::ggtitle("Fiber types") +
  ggplot2::theme(plot.title = ggplot2::element_text(
    size = 20,
    face = "bold"
  )) +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
  ggplot2::theme(legend.title = ggplot2::element_text(
    size = 16,
    face = "bold"
  )) +
  ggplot2::theme(
    axis.title.x = ggplot2::element_text(vjust = -0.35),
    axis.title.y = ggplot2::element_text(vjust = 0.35)
  ) +
  ggplot2::labs(
    x = "Fiber Type",
    y = "Count"
  )


# Fiber_typing_curves_MYHs ------------------------------------------------

# Here I'm basing my final plot of fiber typing curves on Thibaux script:

fiber_types <- list(
  type_1 <- data_fiber_type |>
    dplyr::filter(fiber_type == "Type 1") |>
    dplyr::arrange(desc(MYH7)),
  Hybrid_1_2A <- data_fiber_type |>
    dplyr::filter(fiber_type == "Hybrid 1/2A") |>
    dplyr::arrange(desc(MYH7)),
  type_2A <- data_fiber_type |>
    dplyr::filter(fiber_type == "Type 2A") |>
    dplyr::arrange(MYH2),
  Hybrid_2A_2X <- data_fiber_type |>
    dplyr::filter(fiber_type == "Hybrid 2A/2X") |>
    dplyr::arrange(MYH1)
)
names(fiber_types) <- c(
  "type_1",
  "hybrid_1_2A",
  "type_2A",
  "hybrid_2A_2X"
)

order_matrix <- dplyr::bind_rows(
  fiber_types$type_1,
  fiber_types$hybrid_1_2A,
  fiber_types$type_2A,
  fiber_types$hybrid_2A_2X
) |>
  tibble::add_column("sample_MYH" = seq_len(nrow(data_fiber_type)))

rectangle_colors <- data.frame(
  start = c(0, 331, 414, 787),
  end = c(330, 413, 786, 974),
  fiber_type = c("Type I", "Hybrid I/IIA", "Type IIA", "Hybrid IIA/IIX")
)

order_matrix |>
  ggplot2::ggplot() +
  ggplot2::geom_rect(
    data = rectangle_colors,
    ggplot2::aes(
      xmin = start,
      xmax = end,
      ymin = -Inf,
      ymax = Inf,
      fill = fiber_type
    ),
    alpha = 0.15
  ) +
  ggplot2::scale_fill_manual(values = c(
    "#3B528BFF",
    "#fdc325",
    "#440154FF",
    "#5DC863FF"
  )) +
  ggplot2::theme(legend.position = "none") +
  ggplot2::geom_point(
    ggplot2::aes(
      x = sample_MYH,
      y = MYH7
    ),
    colour = "#440154FF",
    size = 0.25
  ) +
  ggplot2::geom_point(
    ggplot2::aes(
      x = sample_MYH,
      y = MYH2
    ),
    colour = "#5DC863FF",
    size = 0.25
  ) +
  ggplot2::geom_point(
    ggplot2::aes(
      x = sample_MYH,
      y = MYH1
    ),
    colour = "#fdc325",
    size = 0.25
  ) +
  ggplot2::ylab("% MYH isoform expressed") +
  ggplot2::xlab("Sample (ranked by MYH expression)") +
  ggplot2::ggtitle("MYH-based fiber typing (proteomics)") +
  ggplot2::geom_vline(xintercept = 330, colour = "grey30", size = 0.25) +
  ggplot2::geom_vline(xintercept = 413, colour = "grey30", size = 0.25) +
  ggplot2::geom_vline(xintercept = 786, colour = "grey30", size = 0.25) +
  ggplot2::geom_hline(
    yintercept = 31,
    linetype = "dotted",
    colour = "black",
    size = 0.25
  ) +
  ggplot2::geom_hline(
    yintercept = 8.9,
    linetype = "dotted",
    colour = "black",
    size = 0.25
  ) +
  ggplot2::geom_hline(
    yintercept = 7.9,
    linetype = "dotted",
    colour = "black",
    size = 0.25
  ) +
  ggplot2::theme_classic() +
  ggplot2::theme(
    text = ggplot2::element_text(
      face = "bold",
      size = 12,
      colour = "black"
    ),
    strip.text = ggplot2::element_text(colour = "white"),
    strip.background = ggplot2::element_rect(fill = "black"),
    legend.position = "none",
    plot.title = ggplot2::element_text(hjust = 0.5)
  ) +
  ggplot2::annotate(
    "text",
    x = 165,
    y = 115,
    label = "Type 1 \n 330 (33.9%)",
    colour = "black",
    fontface = 2,
    size = 1.7
  ) +
  ggplot2::annotate(
    "text",
    x = 599,
    y = 115,
    label = "Type 2A \n 373 (38.3%)",
    colour = "black",
    fontface = 2,
    size = 1.7
  ) +
  ggplot2::annotate(
    "text",
    x = 880,
    y = 115,
    label = "Hybrid 2A/2X \n 188 (19.3%)",
    colour = "black",
    fontface = 2,
    size = 1.7
  ) +
  ggplot2::annotate(
    "text",
    x = 371.5,
    y = 110,
    label = "Hybrid \n1/2A \n83 (8.5%)",
    colour = "black",
    fontface = 2,
    size = 1.7
  ) +
  ggplot2::annotate(
    "text",
    x = 35,
    y = 95,
    label = "MYH7",
    colour = "#440154FF",
    fontface = 2,
    size = 2
  ) +
  ggplot2::annotate(
    "text",
    x = 750,
    y = 85,
    label = "MYH2",
    colour = "#5DC863FF",
    fontface = 2,
    size = 2
  ) +
  ggplot2::annotate(
    "text",
    x = 950,
    y = 65,
    label = "MYH1",
    colour = "#fdc325",
    fontface = 2,
    size = 2
  ) +
  ggplot2::coord_cartesian(xlim = c(0, 974), clip = "off") +
  ggplot2::annotate(
    "text",
    x = 75,
    y = 35,
    label = "MYH7 threshold",
    colour = "black",
    size = 2
  ) +
  ggplot2::annotate(
    "text",
    x = 75,
    y = 13,
    label = "MYH2 threshold",
    colour = "black",
    size = 2
  ) +
  ggplot2::annotate(
    "text",
    x = 75,
    y = 5.5,
    label = "MYH1 threshold",
    colour = "black",
    size = 2
  ) +
    ggplot2::theme(text = ggplot2::element_text(size = 7)) +
    ggplot2::scale_x_continuous(limits = c(0, 974),
                                expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(0,125),
                                expand = c(0, 0),
                                breaks = c(0, 25, 50, 75, 100),
                                labels = c("0", "25", "50", "75", "100")) +
    ggplot2::theme(axis.title = ggplot2::element_text(size = 7))

ggplot2::ggsave(
  here::here("doc/figures/figure_1/fiber_type_curves_proteomics.png"),
  device = "png",
  width = 128,
  height = 60,
  units = "mm"
)

# Exporting MYH abundances for dot blot-selected fibers -------------------

dot_blot_fibers <- c(
    "P3_fibNumber33",
    "P5_fibNumber175",
    "P3_fibNumber180",
    "P4_fibNumber162",
    "P1_fibNumber192",
    "P2_fibNumber36",
    "P5_fibNumber159",
    "P2_fibNumber90",
    "P2_fibNumber164",
    "P1_fibNumber177",
    "P5_fibNumber35",
    "P5_fibNumber187",
    "P5_fibNumber128",
    "P5_fibNumber78",
    "P4_fibNumber136",
    "P4_fibNumber159",
    "P4_fibNumber188",
    "P4_fibNumber60",
    "P4_fibNumber24",
    "P3_fibNumber108",
    "P4_fibNumber108",
    "P3_fibNumber110",
    "P4_fibNumber182",
    "P2_fibNumber29",
    "P4_fibNumber152",
    "P3_fibNumber177",
    "P2_fibNumber194",
    "P4_fibNumber199",
    "P3_fibNumber88",
    "P3_fibNumber193",
    "P4_fibNumber57",
    "P3_fibNumber100",
    "P1_fibNumber100",
    "P3_fibNumber185",
    "P3_fibNumber80",
    "P1_fibNumber171",
    "P1_fibNumber198",
    "P5_fibNumber199",
    "P4_fibNumber181",
    "P5_fibNumber119"
)

data_dot_blot <- data_fiber_type  |>
    dplyr::filter(fiber_ID %in% dot_blot_fibers)  |>
    dplyr::arrange(match(fiber_ID, dot_blot_fibers)) |>
    dplyr::select(fiber_ID, MYH7, MYH2, MYH1, fiber_type) |>
    dplyr::mutate("sample_order_left_to_right" = rep(c(1:10), 4))

readr::write_csv(data_dot_blot,
                 here::here("data/dot_blot/MS_fiber_type_dot_blot_fibers.csv"))
