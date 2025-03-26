
plot_aggregate_expression <- function(object,gene_names,meta_name = "aggr_exp"){
  # Are all gene names in the Seurat object?
  # Use only those which are, otherwise there is an error!
  gene_query <- gene_names %in% rownames(object@assays$MAGIC_RNA)
  if(sum(gene_query)> 0){
    gene_names <- gene_names[gene_query]
    if(length(gene_names) < 5){
      meta_name <- paste(c(meta_name,gene_names), collapse = ".")
    }
    sums <- rowSums(object@assays$MAGIC_RNA@data[gene_names,])
    gene_names <- gene_names[which(sums!=0)]
    
    #Get the Aggregate expression
    object <- MetaFeature(
      object = object,
      features = gene_names,
      meta.name = meta_name
    )
    FeaturePlot(object,features = meta_name,reduction = "INTE_UMAP",
                min.cutoff = "q05",
                max.cutoff = "q95", 
                label = FALSE, order = TRUE,cols = c("lightgrey", "red"))
    
    VlnPlot(object, 
            features = meta_name,
            pt.size = 0,
            group.by = "stage",
            sort = TRUE) +
      theme(axis.title.x = element_blank()) + NoLegend() 
    
  } #error capture
  else {return(capture.output(cat("Error! Could not find any of the following genes",
                                  gene_names)))
  }
}

plot_aggregate_expression_split <- function(object,gene_names,meta_name = "aggr_exp"){
  # Are all gene names in the Seurat object?
  # Use only those which are, otherwise there is an error!
  gene_query <- gene_names %in% rownames(object@assays$RNA)
  if(sum(gene_query)> 0){
    gene_names <- gene_names[gene_query]
    if(length(gene_names) < 5){
      meta_name <- paste(c(meta_name,gene_names), collapse = ".")
    }
    sums <- rowSums(object@assays$RNA@data[gene_names,])
    gene_names <- gene_names[which(sums!=0)]
    
    #Get the Aggregate expression
    object <- MetaFeature(
      object = object,
      features = gene_names,
      meta.name = meta_name
    )
    
    #feature_list <- list(gene_names)
    
    #scrna <- AddModuleScore(object = scrna, features = feature_list, name = "prot_not_cm")
    
    FeaturePlot(object,features = meta_name,reduction = "INTE_UMAP",
                min.cutoff = "q05",
                max.cutoff = "q95", 
                label = FALSE, order = TRUE,cols = c("lightgrey", "red"),
                split.by = "stage")
  } #error capture
  else {return(capture.output(cat("Error! Could not find any of the following genes",
                                  gene_names)))
  }
}

get_aggregate_expression <- function(object,gene_names,meta_name = "aggr_exp"){
  # Are all gene names in the Seurat object?
  # Use only those which are, otherwise there is an error!
  gene_query <- gene_names %in% rownames(object@assays$RNA)
  if(sum(gene_query)> 0){
    gene_names <- gene_names[gene_query]
    if(length(gene_names) < 5){
      meta_name <- paste(c(meta_name,gene_names), collapse = ".")
    }
    sums <- rowSums(object@assays$RNA@data[gene_names,])
    gene_names <- gene_names[which(sums!=0)]
    
    #Get the Aggregate expression
    object <- MetaFeature(
      object = object,
      features = gene_names,
      meta.name = meta_name
    )
    return(object)
  } #error capture
  else {return(capture.output(cat("Error! Could not find any of the following genes",
                                  gene_names)))
  }
}
# 
# ## and the better, more sophisticated version:
# capwords <- function(s, strict = FALSE) {
#   cap <- function(s) paste(toupper(substring(s, 1, 1)),
#                            {s <- substring(s, 2); if(strict) tolower(s) else s},
#                            sep = "", collapse = " " )
#   sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
# }
