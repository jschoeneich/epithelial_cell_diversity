#helper functions

# define helper functions ----
`%ni%` <- Negate(`%in%`)

save_ggplot_formats = function(
    plt, base_plot_dir, plt_name, create_plot_subdir = TRUE,
    formats = c("png", "pdf"), units = "in", width = 20, height = 20,
    type = "cairo", res = 300,
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
            device = cairo_pdf,
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


# Define a function to generate the plots
generate_split_plot <- function(i, scrna, cols, levels,reduction) {
  current_level <- levels[i]
  #cell_number <- table(scrna$harmony_inte_clusters)[i] # TODO: fix for automation
  #  cell_number <- table(scrna$name)[i] # TODO: fix for automation
  cols.highlight <- cols[i]
  title <- glue("{current_level}") # ({cell_number})")
  DimPlot(scrna, 
          label = FALSE, 
          reduction = reduction, 
          cells.highlight = WhichCells(scrna, idents = current_level),
          cols.highlight = cols.highlight,
          cols = "lightgrey",
          sizes.highlight = 2,
          #order = current_level,
          raster = TRUE) +
    ggtitle(title) + mini_umap_theme & NoLegend() 
} 

getsheets <- function(filename,startrow){
  #get all sheets from an excel file into a list
  sheets <- openxlsx::getSheetNames(filename)
  sheetlist <- lapply(sheets,openxlsx::read.xlsx,xlsxFile = filename,startRow = startrow)
  names(sheetlist) <- sheets
  return(sheetlist)
}