#helper functions

# define helper functions ----
`%ni%` <- Negate(`%in%`)

save_ggplot_formats = function(
    plt, base_plot_dir, plt_name, create_plot_subdir = TRUE,
    formats = c("png", "pdf"), units = "in", width = 20, height = 20,
    type="cairo", res = 300,
    plot_obj ="ggplot",
    ...
){
  # if theres a plot and basedir
  if(!is.null(base_plot_dir) & !is.null(plt)){
    # for each format
    for(fmt in formats){
      f_path_fmt = file.path(base_plot_dir, paste0(plt_name, ".", fmt))
      if(create_plot_subdir) dir.create(file.path(base_plot_dir,fmt),
                                        recursive = TRUE,
                                        showWarnings = FALSE)
      if(dir.exists(file.path(base_plot_dir,fmt))){
        f_path_fmt = file.path(base_plot_dir,fmt,paste0(plt_name,".",fmt))}
      if(plot_obj == "ggplot"){
        if(fmt == "png"){
          ggplot2::ggsave(
            filename = f_path_fmt,
            plot = plt,
            device = fmt,
            units = units,
            width = width,
            height = height,
            type = type, 
            limitsize = FALSE,
            ...
          )
        }else{ #pdf. changed fmt to cairo_pdf
          ggplot2::ggsave(
            filename = f_path_fmt,
            plot = plt,
            #device = cairo_pdf,
            device = pdf,
            units = units,
            width = width,
            height = height, 
            limitsize = FALSE,...
          )
        }
      }else
        if(fmt == "png"){
          png(
            filename = f_path_fmt,
            units = units,
            width = width,
            height = height,
            type = type,
            res = res
          )
          draw(plt, ...)
          dev.off()
        }else{
          pdf(
            file = f_path_fmt,
            width = width,
            height = height
          )
          draw(plt, ...)
          dev.off()
        }
    }
  }
}


# --- flexible split plot generator ---
generate_split_plot <- function(i, 
                                scrna, 
                                cols, 
                                levels, 
                                reduction, 
                                show_cell_numbers = FALSE) {
  current_level <- levels[i]
  
  # get cell number if requested
  cell_number <- if (show_cell_numbers) {
    sum(Idents(scrna) == current_level)
  } else {
    NA
  }
  
  cols.highlight <- cols[i]
  
  # title with or without cell numbers
  title <- if (show_cell_numbers) {
    glue("{current_level} ({cell_number})")
  } else {
    glue("{current_level}")
  }
  
  DimPlot(scrna, 
          label = FALSE, 
          reduction = reduction, 
          cells.highlight = WhichCells(scrna, idents = current_level),
          cols.highlight = cols.highlight,
          cols = "lightgrey",
          sizes.highlight = 0.1) +
    ggtitle(title) + mini_umap_theme & NoLegend() 
}

# --- main function ---
make_cluster_grid <- function(scrna, 
                              meta_column, 
                              reduction = "INTE_UMAP", 
                              plots_folder = "plots",
                              base_size = 3,
                              aspect_ratio = 1,
                              show_cell_numbers = FALSE,
                              prefix = "",
                              plot_cols = NULL) {
  
  Idents(scrna) <- meta_column
  cluster_levels <- levels(as.factor(scrna[[meta_column, drop = TRUE]]))
  n_clusters <- length(cluster_levels)
  
  if(is.null(plot_cols)){
    plot_cols <- Hue_Pal(n_clusters)
  }
  
  # dynamic plot size
  n_cols <- ceiling(sqrt(n_clusters))
  n_rows <- ceiling(n_clusters / n_cols)
  plt_width <- n_cols * base_size
  plt_height <- n_rows * base_size * aspect_ratio
  
  plot_list_clusters <- lapply(
    seq_along(cluster_levels),
    generate_split_plot,
    scrna = scrna,
    levels = cluster_levels,
    reduction = reduction,
    cols = plot_cols,
    show_cell_numbers = show_cell_numbers
  )
  
  plot_grid_clusters <- wrap_plots(plot_list_clusters)
  
  save_ggplot_formats(
    plt = plot_grid_clusters,
    plt_name = glue("plot_grid_{prefix}_{meta_column}"),
    base_plot_dir = plots_folder,
    width = plt_width,
    height = plt_height
  )
  
  return(plot_grid_clusters)
}

sort_scrna_columns <- function(scrna, scrna_columns, order = c("numeric", "alphabetical")) {
  order <- match.arg(order)
  
  scrna@meta.data <- scrna@meta.data %>%
    mutate(across(all_of(scrna_columns),
                  ~ case_when(
                    order == "numeric" ~ fct_reorder(., as.numeric(as.character(.)), .fun = min),
                    order == "alphabetical" ~ fct_relevel(., sort(levels(.)))
                  )))
  
  return(scrna)
}

getsheets <- function(filename,startrow){
  #get all sheets from an excel file into a list
  sheets <- openxlsx::getSheetNames(filename)
  sheetlist <- lapply(sheets,openxlsx::read.xlsx,xlsxFile = filename,startRow = startrow)
  names(sheetlist) <- sheets
  return(sheetlist)
}

generate_prop_table <- function(scrna, table_vars, name, plots_folder) {
  tbl <- table(
    scrna@meta.data[, table_vars[1]],
    scrna@meta.data[, table_vars[2]]
  )
  
  rowsums <- rowSums(tbl)
  tbl <- cbind(tbl, rowsums)
  colsums <- colSums(tbl)
  tbl <- rbind(tbl, colsums)
  
  png(file.path(plots_folder, glue("png/proportion_table_{name}.png")),
      res = 200, width = 1600, height = 450
  )
  p <- tableGrob(tbl, theme = ttheme_default())
  grid.arrange(p)
  dev.off()
  
  pdf(file.path(plots_folder, glue("pdf/proportion_table_{name}.pdf")))
  p <- tableGrob(tbl, theme = ttheme_default())
  grid.arrange(p)
  dev.off()
}

