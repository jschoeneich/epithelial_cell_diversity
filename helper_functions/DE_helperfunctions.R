# define helper functions ----
calculate_results_list_presto <- function(scrna) {
  # Generate all pairwise combinations of timepoints
  combination <- combn(levels(scrna), 2)
  
  # Create list names by combining the pairs with an underscore
  list.names <- capture.output(apply(combination, 2, function(x) {
    cat(x, sep = "_", fill = TRUE)
  }))
  
  # Remove NULL entries
  list.names <- list.names[which(list.names != "NULL")]
  
  # Set identities to 'stage'
  Idents(scrna) <- "stage"
  
  # Initialize results list
  results.list <- sapply(list.names, function(x) NULL)
  
  # Run RunPresto for each combination
  for (i in seq_along(results.list)) {
    results.list[[i]] <- RunPresto(scrna, 
                                   ident.1 = combination[2, i],  # d5
                                   ident.2 = combination[1, i])  # d1
  }
  results.list <- lapply(results.list,function(x){x$symbol <- rownames(x);return(x)})
  # Return the results list
  return(results.list)
}

GeneBarPlot <- function(de.data, xlim = NULL, main = NULL) {
  #de.data = cluster.de[[id]]
  #de.data = plot_de
  if("avg_logFC" %in% names(de.data)){ ## compatible for seurat3
    de.data$avg_log2FC <- de.data$avg_logFC/log(2)
  }
  if (any(colnames(de.data) == "cluster")) {
    top5.up <- de.data %>% group_by(cluster) %>% top_n(10, avg_log2FC) %>% 
      filter(avg_log2FC > 0) %>% arrange(-avg_log2FC)
    top5.dn <- de.data %>% group_by(cluster) %>% top_n(10, -avg_log2FC) %>%
      filter(avg_log2FC < 0) %>% arrange(-avg_log2FC)
  } else {
    top5.up <- de.data  %>% top_n(10, avg_log2FC) %>%
      filter(avg_log2FC > 0) %>% arrange(-avg_log2FC)
    top5.dn <- de.data  %>% top_n(10, -avg_log2FC) %>%
      filter(avg_log2FC < 0) %>% arrange(-avg_log2FC)
  }
  top.up.dn <- rbind(top5.up, top5.dn)
  top.up.dn$gene <- make.unique(top.up.dn$gene)
  top.up.dn$type = ifelse(top.up.dn$avg_log2FC > 0, "positive", "negative")
  top.up.dn$type <- factor(top.up.dn$type, levels = c("positive", "negative"))
  
  g <- ggplot(data = top.up.dn,
              aes(x = gene, y = avg_log2FC, fill = type)) +
    geom_bar(stat="identity") +
    scale_x_discrete(limits = rev(top.up.dn$gene)) +
    theme_minimal() + 
    theme(legend.position="none", axis.text =  element_text(size=15), 
          axis.title.y = element_blank()) +
    #scale_fill_manual(values = c(positive = "#E41A1C", negative = "#377EB8")) +
    scale_fill_manual(values = c(positive = pos_color, negative = neg_color)) +
    coord_flip()
  
  if (!is.null(main)) {
    g <- g + ggtitle(main)
  } else {
    g <- g + ggtitle("Average logFC for the top 5 up and top 5 down regulated genes")
  }
  if (!is.null(xlim)) {
    # Coordinates are flipped
    g <- g + ylim(xlim)
  }
  return(g)
}


# dotplot ----------
# customized function of https://github.com/YuLab-SMU/enrichplot/blob/master/R/dotplot.R
# this function can now take in GeneRatio as a coloring option
# also changed aes_string to aes for deprecation

#' ep_str_wrap internal string wrapping function
#' @param string the string to be wrapped
#' @param width the maximum number of characters before wrapping to a new line
#' @noRd
ep_str_wrap <- function(string, width) {
  x <- gregexpr(' ', string)
  vapply(seq_along(x),
         FUN = function(i) {
           y <- x[[i]]
           n <- nchar(string[i])
           len <- (c(y,n) - c(0, y)) ## length + 1
           idx <- len > width
           j <- which(!idx)
           if (length(j) && max(j) == length(len)) {
             j <- j[-length(j)]
           }
           if (length(j)) {
             idx[j] <- len[j] + len[j+1] > width
           }
           idx <- idx[-length(idx)] ## length - 1
           start <- c(1, y[idx] + 1)
           end <- c(y[idx] - 1, n)
           words <- substring(string[i], start, end)
           paste0(words, collapse="\n")
         },
         FUN.VALUE = character(1)
  )
}

#' default_labeller
#'
#' default labeling function that uses the
#' internal string wrapping function `ep_str_wrap`
#' @noRd
default_labeller <- function(n) {
  function(str){
    str <- gsub("_", " ", str)
    ep_str_wrap(str, n)
  }
}

dotplot.enrichResult2 <- function(object, x = "geneRatio", color = "p.adjust",
                                 showCategory=10, size=NULL, split = NULL,
                                 font.size=12, title = "", orderBy="x",
                                 label_format = 30, decreasing=TRUE) {
  
  colorBy <- match.arg(color, c("pvalue", "p.adjust", "qvalue","GeneRatio"))
  if (x == "geneRatio" || x == "GeneRatio") {
    x <- "GeneRatio"
    if (is.null(size))
      size <- "Count"
  } else if (x == "count" || x == "Count") {
    x <- "Count"
    if (is.null(size))
      size <- "GeneRatio"
  } else if (is(x, "formula")) {
    x <- as.character(x)[2]
    if (is.null(size))
      size <- "Count"
  } else {
    ## message("invalid x, setting to 'GeneRatio' by default")
    ## x <- "GeneRatio"
    ## size <- "Count"
    if (is.null(size))
      size  <- "Count"
  }
  
  df <- fortify(object, showCategory = showCategory, split=split)
  ## already parsed in fortify
  ## df$GeneRatio <- parse_ratio(df$GeneRatio)
  
  if (orderBy !=  'x' && !orderBy %in% colnames(df)) {
    message('wrong orderBy parameter; set to default `orderBy = "x"`')
    orderBy <- "x"
  }
  
  if (orderBy == "x") {
    df <- dplyr::mutate(df, x = eval(parse(text=x)))
  }
  
  label_func <- default_labeller(label_format)
  if(is.function(label_format)) {
    label_func <- label_format
  }
  
  idx <- order(df[[orderBy]], decreasing = decreasing)
  df$Description <- factor(df$Description,
                           levels=rev(unique(df$Description[idx])))
  y_var <- "Description"
  ggplot(df, aes(x = .data[[x]], y = .data[[y_var]], size = .data[[size]],
                 color = .data[[colorBy]])) +
    geom_point() +
    scale_color_continuous(low ="blue", high ="red", name = color,
                           guide = guide_colorbar(reverse = TRUE)) +
    scale_y_discrete(labels = label_func) +
    ylab(NULL) + ggtitle(title) + theme_dose(font.size) +
    scale_size(range = c(3, 8)) +
    guides(size  = guide_legend(order = 1), 
           color = guide_colorbar(order = 2))
}

#function to calculate a fisher.test on a given countmatrix for tempDE and a given gene.co.exp vector (e.g. blimp,mafb)
fisher_dataframe <- function(countmatrix,cm2meanspattern.split,gene.co.exp_vec){
  require(dplyr)
  require(glue)
  require(rstatix)
  df <- data.frame(
    gene.co.exp = sapply(cm2meanspattern.split,function(x){
      sum(rownames(x) %in% gene.co.exp_vec)})
  )
  
  df$n.genes <- sapply(means.split.new,nrow)
  
  df$not.co.exp <- df$n.genes - df$gene.co.exp
  df$not.pattern <- length(gene.co.exp_vec) - df$gene.co.exp
  df$not.co.exp.not.pattern <- nrow(countmatrix) - length(gene.co.exp_vec) - df$gene.co.exp
  
  df <- df %>% relocate(n.genes,.after = last_col())
  df$p.value <- NA
  
  for(i in 1:nrow(df)){
    m <- matrix(unlist(df[i,1:4]),nrow = 2,byrow = T)
    ft <- fisher.test(m,alternative = "greater")
    df$p.value[[i]] <- ft$p.value
  }
  
  df$p_val_adj <- p.adjust(df$p.value,method = "BH")
  df <- df %>% add_significance(p.col = "p_val_adj")
  
  var.name <- substitute(gene.co.exp_vec)
  colnames(df)[c(1,2,4)] <- c(var.name,glue("not.{var.name}"),glue("not.{var.name}.not.pattern"))
  
  #df_filter <- df %>% filter(p_val_adj < 0.05 )
  #return(df_filter)
  return(df)
}

calc_ht_size = function(ht, unit = "inch") {
  pdf(NULL)
  ht = draw(ht)
  w = ComplexHeatmap:::width(ht)
  w = convertX(w, unit, valueOnly = TRUE)
  h = ComplexHeatmap:::height(ht)
  h = convertY(h, unit, valueOnly = TRUE)
  dev.off()
  
  c(w, h)
}
#GO helper functions

# Function to truncate a GO term to the last space before maximum length
truncate_go_term <- function(term,max_length = 30) {
  term <- as.character(term)
  # Check if the term length is already within the limit
  if (nchar(term) <= max_length) {
    return(term)
  } else {
    truncated_term <- substr(term, 1, max_length)
    return(paste0(truncated_term, "[...]"))
  }
}

# Function to truncate a GO term to the last space before maximum length
truncate_go_term_on_last_word <- function(term,max_length = 33) {
  term <- as.character(term)
  # Check if the term length is already within the limit
  if (nchar(term) <= max_length) {
    return(term)
  } else {
    # Find the last space before the maximum length
    last_space <- max(gregexpr("\\s", substr(term, 1, max_length))[[1]])
    truncated_term <- substr(term, 1, last_space -1)
    return(paste0(truncated_term, "[...]"))
  }
}

rename_columns <- function(df) {
  colnames_to_replace <- c("log2FoldChange" = "avg_log2FC", 
                           "logFC" = "avg_log2FC", 
                           "adj.P.Val" = "p_val_adj",  
                           "padj" = "p_val_adj", 
                           "gene" = "symbol",
                           "gene_name" = "symbol")
  
  # Check and replace only matching column names
  colnames(df) <- ifelse(colnames(df) %in% names(colnames_to_replace),
                         colnames_to_replace[colnames(df)],
                         colnames(df))
  
  return(as.data.frame(df))
}


generate_volcano_plots <- function(results_list, 
                                   pos_color,
                                   neg_color, 
                                   plots_folder, 
                                   charts_folder,
                                   filename_prefix,
                                   do_save = FALSE,
                                   subset_condition = TRUE) {
  #contribution_df
  #rename colnames, as different packages give different names
  
  results_list <- lapply(results_list, rename_columns)

  # Precompute the quantile values across all result lists
  quantile_values <- sapply(results_list, function(x) {
    #colnames(x)[c(2,5)] <- c("avg_log2FC","p_val_adj")
    quantile(x$avg_log2FC, c(0.001, 0.999))
  })

  xmin <- min(quantile_values) + 1
  xmax <- max(quantile_values) - 1

  # Precompute the p-value cutoff across all lists
  p_cutoff <- min(sapply(results_list, function(x) {
    #colnames(x)[c(2,5)] <- c("avg_log2FC","p_val_adj")
    quantile(x$p_val_adj[x$p_val_adj > 0], 0.001)
  }))
  
  results_list <- results_list[subset_condition]
  
  # Define a function to process each individual result set
  process_result <- function(condition) {
    res <- as.data.frame(results_list[[condition]])
    
    #colnames(res)[c(2,5)] <- c("avg_log2FC","p_val_adj")
    days <- strsplit(condition, "_")[[1]]
    dayneg <- days[1]
    daypos <- days[2]
    
    #top 20
    top20 <- res %>% slice_max(avg_log2FC, n = 10) %>% arrange(-avg_log2FC)
    
    #bot 20
    bot20 <- res %>% slice_min(avg_log2FC,
                               n = 10) %>% arrange(-avg_log2FC)
    
    #Only those which are significant
    top20 <- subset(top20, subset = p_val_adj < 0.05 & avg_log2FC > 1)[1:10,]
    bot20 <- subset(bot20, subset = p_val_adj < 0.05 & avg_log2FC < - 1)[1:10,]
    
    # Vectorized adjustments to avg_log2FC and p_val_adj
    res$avg_log2FC <- pmin(pmax(res$avg_log2FC, xmin), xmax)
    res$p_val_adj <- pmax(res$p_val_adj, p_cutoff)

    # Identify genes that pass the thresholds
    up <- rownames(subset(res, avg_log2FC >= xmax))
    down <- rownames(subset(res, avg_log2FC <= xmin))
    p_up <- rownames(subset(res, p_val_adj <= p_cutoff))

    # Custom shapes for points
    custom_shape <- rep(19, nrow(res)) # Default shape: solid circle
    names(custom_shape) <- rep('normal', nrow(res))

    # Apply custom shapes based on conditions
    custom_shape[which(rownames(res) %in% up)] <- -9658 # Right arrow
    custom_shape[which(rownames(res) %in% p_up)] <- 17   # Triangle up
    custom_shape[which(rownames(res) %in% down)] <- -9668 # Left arrow

    # Custom sizes for points
    custom_size <- rep(2.0, nrow(res)) # Default size
    custom_size[which(rownames(res) %in% up)] <- 4
    custom_size[which(rownames(res) %in% p_up)] <- 4
    custom_size[which(rownames(res) %in% down)] <- 4

    # Custom colors for points
    keyvals <- ifelse(
      res$avg_log2FC <= -1 & res$p_val_adj <= 0.05, neg_color,
      ifelse(res$avg_log2FC  >= 1 & res$p_val_adj <= 0.05, pos_color, 'grey')
    )
    
    keyvals[is.na(keyvals)] <- 'grey'
    names(keyvals)[keyvals == neg_color] <- glue('up in {dayneg}')
    names(keyvals)[keyvals == 'grey'] <- 'ns'
    names(keyvals)[keyvals == pos_color] <- glue("up in {daypos}")
    
    # Count the significant genes
    up_genes <- sum(res$avg_log2FC >= 1 & res$p_val_adj <= 0.05)
    down_genes <- sum(res$avg_log2FC <= -1 & res$p_val_adj <= 0.05)
    
    # contrib_genes <- as.data.frame(res) %>% 
    #   filter(symbol %in% rownames(contribution_df)) %>% filter(
    #     abs(avg_log2FC) >= 1 & p_val_adj <= 0.05) %>% pull(symbol)
    
    #contrib_intersect <- intersect(c(bot20$symbol,top20$symbol), contrib_genes)
    
    # Generate the volcano plot
    plt <- EnhancedVolcano(res,
                           lab = res$symbol,
                           x = 'avg_log2FC',
                           y = 'p_val_adj',
                           xlab = "Log2 fold change",
                           ylab = "-Log10 p.adj.",
                           selectLab = c(bot20$symbol,top20$symbol),
                           #selectLab = contrib_genes,
                           #selectLab = contrib_intersect,
                           drawConnectors = FALSE,
                           #drawConnectors = TRUE,
                           widthConnectors = 0.5,
                           arrowheads = FALSE,
                           max.overlaps = Inf,
                           colCustom = keyvals,
                           pCutoff = 0.05,
                           FCcutoff = 1.0,
                           boxedLabels = FALSE,
                           colAlpha = 0.5,
                           labSize = 2,
                           pointSize = 1,
                           axisLabSize = 5,
                           titleLabSize = 5,
                           subtitleLabSize = 5,
                           captionLabSize = 5,
                           legendLabSize = 5,
                           legendIconSize = 2,
                           #legendPosition = "right",
                           #shapeCustom = custom_shape,
                           #pointSize = custom_size,
                           #title = glue("{daypos} vs. {dayneg}"),
                           #subtitle = glue("{dayneg}:{down_genes} {daypos}:{up_genes}"),
                           #caption = paste0("total = ", nrow(res), " genes"),
                           title = NULL,
                           subtitle = NULL,
                           caption = NULL,
                           legendLabels = NULL,
                           raster = TRUE,
                           xlim = c(xmin,xmax)
                           #xlim = c(-6,6),
                           #ylim =c(0,7)
                           ) + theme(axis.line = element_line(linewidth = 0.5),
                                     axis.ticks = element_blank(),
                                     legend.position = "none",
                                     legend.background = element_blank(),
                                     panel.background = element_blank(), #transparent panel bg
                                     plot.background = element_blank(), #transparent plot bg
                                     panel.border = element_blank()
                                     )
    
    if(do_save == TRUE){
    # Save the individual plot with the custom filename prefix
    save_ggplot_formats(
      plt = plt,
      base_plot_dir = plots_folder,
      plt_name = glue("{filename_prefix}_Volcano_{daypos}_vs_{dayneg}",
                      "_quantile_cutoff_{condition}"),
      width = 2.5, height = 2.5  # Updated width and height
    )
    
    # Save the filtered DE genes to CSV and Excel
    df2 <- res %>% filter((avg_log2FC >= 1 & p_val_adj <= 0.05) | 
                            (avg_log2FC <= -1 & p_val_adj <= 0.05))
    
     write.csv(df2, file = 
                 file.path(charts_folder,
                                     glue("{filename_prefix}_DE_{daypos}",
                                          "_vs_{dayneg}_{condition}.csv")), 
               row.names = FALSE)
     write.xlsx(df2, file = file.path(
       charts_folder, 
       glue("{filename_prefix}_DE_{daypos}",
            "_vs_{dayneg}_{condition}.xlsx")), 
       rowNames = FALSE)
    }
    # Append the plot to the list
    return(plt)
  }
  # Function to determine possible rows and columns based on plot count
   adjust_plot_layout <- function(plot_count) {
    layout_options <- list()

    # Case for less than 4 plots: single row
    if (plot_count < 4) {
      layout_options[[1]] <- list(rows = 1, cols = plot_count)
    } else {
      # Case for 4 or more plots: calculate possible row x column combinations
      for (rows in 1:plot_count) {
        cols <- ceiling(plot_count / rows)
        if (rows * cols >= plot_count) {
          layout_options[[length(layout_options) + 1]] <- list(rows = rows, cols = cols)
        }
      }
    }

    # For each layout option, calculate width and height (multiples of 7)
    layout_dimensions <- lapply(layout_options, function(layout) {
      width <- layout$cols * 2.5
      height <- layout$rows * 2.5
      return(list(rows = layout$rows, cols = layout$cols, width = width, height = height))
    })

    return(layout_dimensions)
  }
  
  # Example application based on the number of plots
  plot_list <- lapply(names(results_list), function(condition) {
    process_result(condition)
  })
  
  # Get the number of plots
  plot_count <- length(plot_list)

  # Get dynamic layout options based on plot count
  layout_info <- adjust_plot_layout(plot_count)

  if(do_save == TRUE){
  # Optionally, choose a specific layout or iterate over multiple layout versions
  for (i in seq_along(layout_info)) {
    current_layout <- layout_info[[i]]

    # Create combined volcano plot with the current layout
    combined_plot <- wrap_plots(plot_list) +
      plot_layout(ncol = current_layout$cols)
        # Save the combined plot with current layout dimensions
    save_ggplot_formats(
      plt = combined_plot,
      base_plot_dir = plots_folder,
      plt_name = glue("{filename_prefix}_Combined_Volcano_Plots_Layout_{i}"),
      width = current_layout$width,
      height = current_layout$height
    )
   }
  } else{ return(plot_list)
     }
}
# Example usage:
# generate_volcano_plots(results_list = results.list, 
#                        pos_color = "#E41A1C", 
#                        neg_color = "#377EB8", 
#                        plots_folder = "path/to/plots_folder", 
#                        charts_folder = "path/to/charts_folder")



generate_volcano_plots_mark_unique <- function(results_list,
                                               unique_genes_list, 
                                               pos_color,
                                               neg_color,
                                               mark_color = "purple",
                                               plots_folder, 
                                               charts_folder, 
                                               filename_prefix,
                                               do_save = FALSE,
                                               subset_condition = TRUE) {
  
  results_list <- lapply(results_list, rename_columns)
  
  # Precompute the quantile values across all result lists
  quantile_values <- sapply(results_list, function(x) {
    #colnames(x)[c(2,5)] <- c("avg_log2FC","p_val_adj")
    quantile(x$avg_log2FC, c(0.001, 0.999))
  })
  
  xmin <- min(quantile_values) + 1
  xmax <- max(quantile_values) - 1
  
  # Precompute the p-value cutoff across all lists
  p_cutoff <- min(sapply(results_list, function(x) {
    #colnames(x)[c(2,5)] <- c("avg_log2FC","p_val_adj")
    quantile(x$p_val_adj[x$p_val_adj > 0], 0.001)
  }))
  
  results_list <- results_list[subset_condition]

  # Define a function to process each individual result set
  process_result <- function(condition, unique_genes) {
    res <- as.data.frame(results_list[[condition]])
    
    #colnames(res)[c(2,5)] <- c("avg_log2FC","p_val_adj")
    days <- strsplit(condition, "_")[[1]]
    # dayneg <- "d05"
    # daypos <- "d05i"
    
    #top 20
    top20 <- res %>% top_n(10,avg_log2FC) %>% arrange(-avg_log2FC)
    
    #bot 20
    bot20 <- res %>% slice_min(avg_log2FC,n = 10) %>% arrange(-avg_log2FC)
    
    #Only those which are significant
    top20 <- subset(top20, subset = p_val_adj < 0.05 & avg_log2FC > 1)[1:10,]
    bot20 <- subset(bot20, subset = p_val_adj < 0.05 & avg_log2FC < - 1)[1:10,]
    
    # Vectorized adjustments to avg_log2FC and p_val_adj
    res$avg_log2FC <- pmin(pmax(res$avg_log2FC, xmin), xmax)
    res$p_val_adj <- pmax(res$p_val_adj, p_cutoff)
    
    # Identify genes that pass the thresholds
    up <- rownames(subset(res, avg_log2FC >= xmax))
    down <- rownames(subset(res, avg_log2FC <= xmin))
    p_up <- rownames(subset(res, p_val_adj <= p_cutoff))
    
    # Custom shapes for points
    custom_shape <- rep(19, nrow(res)) # Default shape: solid circle
    names(custom_shape) <- rep('normal', nrow(res))
    
    # Apply custom shapes based on conditions
     # custom_shape[which(rownames(res) %in% up)] <- -9658 # Right arrow
     # custom_shape[which(rownames(res) %in% p_up)] <- 17   # Triangle up
     # custom_shape[which(rownames(res) %in% down)] <- -9668 # Left arrow
    
    # Custom sizes for points
     #custom_size <- rep(2, nrow(res)) # Default size
     #custom_size[which(rownames(res) %in% up)] <- 4
     #custom_size[which(rownames(res) %in% down)] <- 4
    
    
    # Count the significant genes
    up_genes <- sum(res$avg_log2FC >= 1 & res$p_val_adj <= 0.05)
    down_genes <- sum(res$avg_log2FC <= -1 & res$p_val_adj <= 0.05)
    unique_labels <- res %>% filter(avg_log2FC >= 1 & p_val_adj <= 0.05 &
                                      symbol %in% unique_genes) %>% pull(symbol)
    
    # Custom colors for points
    keyvals <- ifelse(
      res$avg_log2FC <= -1 & res$p_val_adj <= 0.05, neg_color,
      ifelse(res$avg_log2FC >= 1 & res$p_val_adj <= 0.05, pos_color,'grey')
    )
    
    keyvals[which(res$symbol %in% unique_labels)] <- mark_color
    
    keyvals[is.na(keyvals)] <- 'grey'
    names(keyvals)[keyvals == neg_color] <- glue('up in {dayneg}')
    names(keyvals)[keyvals == 'grey'] <- 'ns'
    names(keyvals)[keyvals == mark_color] <- glue("unique {condition} infected")
    names(keyvals)[keyvals == pos_color] <- glue("up in {daypos}")
    
    # Generate the volcano plot
    plt <- EnhancedVolcano(res,
                           lab = res$symbol,
                           x = 'avg_log2FC',
                           y = 'p_val_adj',
                           xlab = "Log2 fold change",
                           ylab = "-Log10 p.adj.",
                           selectLab = c(bot20$symbol,top20$symbol),
                           #selectLab = contrib_genes,
                           #selectLab = contrib_intersect,
                           #drawConnectors = FALSE,
                           drawConnectors = TRUE,
                           widthConnectors = 0.5,
                           arrowheads = FALSE,
                           max.overlaps = Inf,
                           colCustom = keyvals,
                           pCutoff = 0.05,
                           FCcutoff = 1.0,
                           boxedLabels = FALSE,
                           colAlpha = 0.5,
                           labSize = 2,
                           pointSize = 1,
                           axisLabSize = 5,
                           titleLabSize = 5,
                           subtitleLabSize = 5,
                           captionLabSize = 5,
                           legendLabSize = 5,
                           legendIconSize = 2,
                           #legendPosition = "right",
                           #shapeCustom = custom_shape,
                           #pointSize = custom_size,
                           #title = glue("{daypos} vs. {dayneg}"),
                           #subtitle = glue("{dayneg}:{down_genes} {daypos}:{up_genes}"),
                           #caption = paste0("total = ", nrow(res), " genes"),
                           title = NULL,
                           subtitle = NULL,
                           caption = NULL,
                           legendLabels = NULL,
                           raster = TRUE,
                           xlim = c(xmin,xmax)
                           #xlim = c(-6,6),
                           #ylim =c(0,7)
    ) + theme(axis.line = element_line(linewidth = 0.5),
              axis.ticks = element_blank(),
              legend.position = "none",
              legend.background = element_blank(),
              panel.background = element_blank(), #transparent panel bg
              plot.background = element_blank(), #transparent plot bg
              panel.border = element_blank()
    )
    
    if(do_save == TRUE){
      # Save the individual plot with the custom filename prefix
      save_ggplot_formats(
        plt = plt,
        base_plot_dir = plots_folder,
        plt_name = glue("{filename_prefix}_Volcano_{daypos}_vs_{dayneg}_quantile_cutoff_{condition}_marked"),
        width = 2.5, height = 2.5  # Updated width and height
      )
      
      # Save the filtered DE genes to CSV and Excel
      df2 <- res %>% filter((avg_log2FC >= 1 & p_val_adj <= 0.05) | 
                              (avg_log2FC <= -1 & p_val_adj <= 0.05))
      
      write.csv(df2, file = file.path(charts_folder, 
                                      glue("{filename_prefix}_DE_{daypos}",
                                           "_vs_{dayneg}_{condition}_marked.csv")),
                row.names = FALSE)
      write.xlsx(df2, file = file.path(charts_folder,
                                       glue("{filename_prefix}_DE_{daypos}",
                                            "_vs_{dayneg}_{condition}:marked.xlsx")), 
                 rowNames = FALSE)
    }
    # Append the plot to the list
    return(plt)
  }
  # Function to determine possible rows and columns based on plot count
  adjust_plot_layout <- function(plot_count) {
    layout_options <- list()
    
    # Case for less than 4 plots: single row
    if (plot_count < 4) {
      layout_options[[1]] <- list(rows = 1, cols = plot_count)
    } else {
      # Case for 4 or more plots: calculate possible row x column combinations
      for (rows in 1:plot_count) {
        cols <- ceiling(plot_count / rows)
        if (rows * cols >= plot_count) {
          layout_options[[length(layout_options) + 1]] <- list(rows = rows, cols = cols)
        }
      }
    }
    
    # For each layout option, calculate width and height (multiples of 7)
    layout_dimensions <- lapply(layout_options, function(layout) {
      width <- layout$cols * 2.5
      height <- layout$rows * 2.5
      return(list(rows = layout$rows, 
                  cols = layout$cols, 
                  width = width, 
                  height = height)
             )
    })
    
    return(layout_dimensions)
  }
  
  # Example application based on the number of plots
  plot_list <- lapply(names(results_list), function(condition) {
    unique_genes <- unlist(unique_genes_list)  # Fetch unique genes for each condition
    process_result(condition, unique_genes)
  })
  
  # Get the number of plots
  plot_count <- length(plot_list)
  
  # Get dynamic layout options based on plot count
  layout_info <- adjust_plot_layout(plot_count)
  
  if(do_save == TRUE){
    # Optionally, choose a specific layout or iterate over multiple layout versions
    for (i in seq_along(layout_info)) {
      current_layout <- layout_info[[i]]
      
      # Create combined volcano plot with the current layout
      combined_plot <- wrap_plots(plot_list) +
        plot_layout(ncol = current_layout$cols)
      # Save the combined plot with current layout dimensions
      save_ggplot_formats(
        plt = combined_plot,
        base_plot_dir = plots_folder,
        plt_name = glue("{filename_prefix}_Combined_Volcano_Layout_marked_{i}"),
        width = current_layout$width,
        height = current_layout$height
      )
    }
  } else{ return(plot_list)
  }
}



# Function to get gene descriptions using biomaRt
get_gene_description <- function(symbols) {
  require(biomaRt)
  # Connect to Ensembl biomart
  mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  # Query biomart for gene descriptions
  result <- getBM(
    attributes = c("mgi_symbol", "description"),
    filters = "mgi_symbol",
    values = symbols,
    mart = mart
  )
  
  # Create a named vector for mapping
  description_map <- setNames(result$description, result$mgi_symbol)
  
  # Return descriptions or "Unknown" for missing values
  sapply(symbols, function(sym) description_map[[sym]] %||% "Unknown")
}
