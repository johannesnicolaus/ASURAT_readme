#--------------------------------------------------------------------------------
#
#--------------------------------------------------------------------------------
process_001_seurat <- function(obj, nfeatures){
  #--------------------------------------------------
  # Normalize the data
  #--------------------------------------------------
  obj <- NormalizeData(obj, normalization.method = "LogNormalize")
  #--------------------------------------------------
  # Perform a variance stabilizing transform (VST).
  # As mentioned in Cruz and Wishart, Cancer Inform. 2, 59-77, 2006.
  # the sample-per-variable feature ratio is set as 5:1.
  # n <- round(0.2 * obj@assays[["RNA"]]@counts@Dim[1])
  #--------------------------------------------------
  n <- nfeatures
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = n)
  #--------------------------------------------------
  # Scale the data
  #--------------------------------------------------
  obj <- ScaleData(obj)
  #--------------------------------------------------
  # Principal component analysis
  #--------------------------------------------------
  obj <- RunPCA(obj, features = VariableFeatures(obj))

  return(obj)
}
#--------------------------------------------------------------------------------
#
#--------------------------------------------------------------------------------
process_002_seurat <- function(obj, pc, resolution){
  set.seed(8)
  #--------------------------------------------------
  # Sample clustering
  #--------------------------------------------------
  obj <- FindNeighbors(obj, reduction = "pca", dim = 1:pc)
  obj <- FindClusters(obj, resolution = resolution)
  #--------------------------------------------------
  # t-SNE
  #--------------------------------------------------
  obj <- RunTSNE(
    obj, dims.use = 1:2, reduction = "pca", dims = 1:pc,
    do.fast = FALSE, perplexity = 30
  )
  #--------------------------------------------------
  # UMAP
  #--------------------------------------------------
  obj <- RunUMAP(obj, dims = 1:pc)

  return(obj)
}
#--------------------------------------------------------------------------------
#
#--------------------------------------------------------------------------------
do_enrichGO_seurat <- function(obj, qvalue_cutoff, orgdb){
  #--------------------------------------------------
  # Definitions
  #--------------------------------------------------
  library(clusterProfiler)
  tmp <- obj@misc[["markers"]]
  category_names <- c("MF", "BP", "CC", "ALL")
  res <- list()
  cluster_names <- as.character(unique(tmp$cluster))
  #--------------------------------------------------
  # enrichGO
  #--------------------------------------------------
  for(cluster in cluster_names){
    for(category in category_names){
      df <- tmp[which(tmp$cluster == cluster),]
      goi <- df[which(df$p_val_adj <= qvalue_cutoff), ]$gene
      goi <- bitr(goi, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = orgdb)
      ego <- enrichGO(
        gene = goi$ENTREZID, OrgDb = orgdb, ont = category,
        pAdjustMethod = "BH", qvalueCutoff = 0.05,  readable = TRUE
      )
      res[[cluster]][[category]] <- ego
    }
  }
  obj@misc[["enrichGO"]] <- res

  return(obj)
}
#--------------------------------------------------------------------------------
#
#--------------------------------------------------------------------------------
do_simplifyGO_seurat <- function(obj, cutoff){
  #--------------------------------------------------
  # Definitions
  #--------------------------------------------------
  library(clusterProfiler)
  tmp <- obj@misc[["enrichGO"]]
  cluster_names <- names(tmp)
  #--------------------------------------------------
  # simplify
  #--------------------------------------------------
  for(cluster in cluster_names){
    df <- obj@misc[["enrichGO"]][[cluster]]
    category_names_woALL <- setdiff(names(df), "ALL")
    for(category in category_names_woALL){
      obj@misc[["enrichGO"]][[cluster]][[category]] <- simplify(
        x = df[[category]],
        cutoff = cutoff,
        by = "p.adjust",
        select_fun = min,
        measure = "Wang",
        semData = NULL
      )
    }
  }
  
  return(obj)
}
#--------------------------------------------------------------------------------
#
#--------------------------------------------------------------------------------
plot_variableFeatures_seurat <- function(obj, n, title, title_size, xlabel, ylabel){
  #--------------------------------------------------
  # Definitions
  #--------------------------------------------------
  top_n <- head(VariableFeatures(obj), n = n)
  #--------------------------------------------------
  # Plot
  #--------------------------------------------------
  p <- LabelPoints(plot=VariableFeaturePlot(obj), points=top_n, repel=TRUE, xnudge=0, ynudge=0) +
    labs(title=title, x=xlabel, y=ylabel) +
    theme_classic(base_size=16, base_family="Helvetica") +
    theme(plot.title=element_text(hjust=0.5, size=title_size),
      plot.margin=unit(c(0.5,2,0.5,0.5), "lines"), legend.position="right")

  return(p)
}
#--------------------------------------------------------------------------------
#
#--------------------------------------------------------------------------------
plot_pcaEigen_seurat <- function(obj, title, title_size, xlabel, ylabel){
  #--------------------------------------------------
  # Definitions
  #--------------------------------------------------
  pca_eigen <- apply(obj@reductions[["pca"]]@cell.embeddings, 2, sd)
  I <- length(pca_eigen)
  df <- data.frame(dim = 1:I, eigen = pca_eigen)
  #--------------------------------------------------
  # Plot
  #--------------------------------------------------
  p <- ggplot() +
    geom_point(
      aes(x=df$dim, y=df$eigen),
      shape=21, color="black", fill="white", alpha=1, size=1, stroke=1
    ) +
    labs(title=title, x=xlabel, y=ylabel) +
    theme_classic(base_size=20, base_family="Helvetica") +
    theme(plot.title=element_text(hjust=0.5, size=title_size),
      plot.margin=unit(c(0.5,2,0.5,0.5), "lines"), legend.position="none")

  return(p)
}
#--------------------------------------------------------------------------------
#
#--------------------------------------------------------------------------------
plot_bargraph_seurat <- function(obj, title, title_size, xlabel, ylabel, ymax){
  #--------------------------------------------------
  # Set data frame to be output
  #--------------------------------------------------
  tmp <- as.data.frame(obj@meta.data[["seurat_clusters"]])
  names(tmp) <- "Count"
  df <- tmp %>% group_by(Count) %>% summarise (n=n())
  #--------------------------------------------------
  # Plot
  #--------------------------------------------------
  p <- ggplot(df) +
    geom_bar(aes(x=as.factor(Count), y=n),
      color="black", fill="black", stat="identity", width=0.7) +
    theme_classic(base_size=20, base_family="Helvetica") +
    theme(plot.title=element_text(hjust=0.5, size=title_size)) +
    geom_text(aes(x=as.factor(Count), y=n, label=sprintf("%d", n)), size=6, vjust=-0.5) +
    ylim(0,ymax) + labs(title=title, x=xlabel, y=ylabel)

  return(p)
}
#--------------------------------------------------------------------------------
#
#--------------------------------------------------------------------------------
plot_tsne_seurat <- function(
  obj, title, title_size, xlabel, ylabel, default_color = TRUE
){
  #--------------------------------------------------
  # Definitions
  #--------------------------------------------------
  obj@meta.data[["seurat_clusters"]] <- as.factor(as.integer(
    as.character(obj@meta.data[["seurat_clusters"]])) + 1)
  df <- data.frame(
    x = obj@reductions[["tsne"]]@cell.embeddings[,1],
    y = obj@reductions[["tsne"]]@cell.embeddings[,2],
    label = obj[["seurat_clusters"]]
  )
  colnames(df) <- c("x", "y", "label")
  n_groups <- length(unique(obj$seurat_clusters))
  #--------------------------------------------------
  # Positions of labels
  #--------------------------------------------------
  cluster_num <- as.integer(as.character(sort(unique(df$label))))
  cog <- c()  # Center of gravity
  for(i in cluster_num){
    cog <- rbind(cog, cbind(
      mean(df[which(df$label == i),]$x), mean(df[which(df$label == i),]$y)))
  }
  cog <- data.frame(x = cog[,1], y = cog[,2], label = cluster_num)
  #--------------------------------------------------
  # Color
  #--------------------------------------------------
  if(default_color){
    #------------------------------
    # Option 1: ggplot's default
    #------------------------------
    n_groups <- length(unique(df$label))
    ggColorHue <- function(n, l=65){   # To produce default colors of ggplot2
      hues <- seq(15, 375, length=n+1)
      hcl(h=hues, l=l, c=100)[1:n]
    }
    my_colors <- ggColorHue(n_groups)
  }else{
    #------------------------------
    # Option 2: rainbow colors
    #------------------------------
    n_groups <- length(unique(df$label))
    my_colors <- rainbow(n_groups)
  }
  #--------------------------------------------------
  # Plot
  #--------------------------------------------------
  p <- ggplot() +
    geom_point(aes(x=df[,1], y=df[,2], color=df[,3]), size=0.5, alpha=1.0) +
    labs(title=title, x=xlabel, y=ylabel, colour="") +
    scale_colour_manual(values=my_colors) +
    theme_classic(base_size=20, base_family="Helvetica") +
    theme(plot.title=element_text(hjust=0.5, size=title_size),
      plot.margin=unit(c(0.5,2,0.5,0.5), "lines"), legend.position="right") +
    geom_text(aes(x=cog[,1], y=cog[,2], label=sprintf("%d", cluster_num)), size=6)

  return(p)
}
#--------------------------------------------------------------------------------
#
#--------------------------------------------------------------------------------
plot_umap_seurat <- function(
  obj, title, title_size, xlabel, ylabel, default_color = TRUE
){
  #--------------------------------------------------
  # Definitions
  #--------------------------------------------------
  obj@meta.data[["seurat_clusters"]] <- as.factor(as.integer(
    as.character(obj@meta.data[["seurat_clusters"]])) + 1)
  df <- data.frame(
    x = obj@reductions[["umap"]]@cell.embeddings[,1],
    y = obj@reductions[["umap"]]@cell.embeddings[,2],
    label = obj[["seurat_clusters"]]
  )
  colnames(df) <- c("x", "y", "label")
  n_groups <- length(unique(obj$seurat_clusters))
  #--------------------------------------------------
  # Positions of labels
  #--------------------------------------------------
  cluster_num <- as.integer(as.character(sort(unique(df$label))))
  cog <- c()  # Center of gravity
  for(i in cluster_num){
    cog <- rbind(cog, cbind(
      mean(df[which(df$label == i),]$x), mean(df[which(df$label == i),]$y)))
  }
  cog <- data.frame(x = cog[,1], y = cog[,2], label = cluster_num)
  #--------------------------------------------------
  # Color
  #--------------------------------------------------
  if(default_color){
    #------------------------------
    # Option 1: ggplot's default
    #------------------------------
    n_groups <- length(unique(df$label))
    ggColorHue <- function(n, l=65){   # To produce default colors of ggplot2
      hues <- seq(15, 375, length=n+1)
      hcl(h=hues, l=l, c=100)[1:n]
    }
    my_colors <- ggColorHue(n_groups)
  }else{
    #------------------------------
    # Option 2: rainbow colors
    #------------------------------
    n_groups <- length(unique(df$label))
    my_colors <- rainbow(n_groups)
  }
  #--------------------------------------------------
  # Plot
  #--------------------------------------------------
  p <- ggplot() +
    geom_point(aes(x=df[,1], y=df[,2], color=df[,3]), size=0.5, alpha=1.0) +
    labs(title=title, x=xlabel, y=ylabel, colour="") +
    scale_colour_manual(values=my_colors) +
    theme_classic(base_size=20, base_family="Helvetica") +
    theme(plot.title=element_text(hjust=0.5, size=title_size),
      plot.margin=unit(c(0.5,2,0.5,0.5), "lines"), legend.position="right") +
    geom_text(aes(x=cog[,1], y=cog[,2], label=sprintf("%d", cluster_num)), size=10)

  return(p)
}
#--------------------------------------------------------------------------------
#
#--------------------------------------------------------------------------------
plot_tsne_nReads_seurat <- function(obj, title, title_size, xlabel, ylabel){
  #--------------------------------------------------
  # Definitions
  #--------------------------------------------------
  df <- data.frame(
    x = obj@reductions[["tsne"]]@cell.embeddings[,1],
    y = obj@reductions[["tsne"]]@cell.embeddings[,2],
    nReads = log(obj@meta.data[["nCount_RNA"]])
  )
  stdev <- sd(df$nReads)
  df$nReads <- (df$nReads - mean(df$nReads)) / stdev
  colnames(df) <- c("x", "y", "nReads")
  #--------------------------------------------------
  # Plot
  #--------------------------------------------------
  p <- ggplot() +
    geom_point(aes(x=df[,1], y=df[,2], color=df[,3]), size=0.5, alpha=1.0) +
    labs(title=title, x=xlabel, y=ylabel, color="z-log") +
    scale_colour_gradientn(colours=c("gray90","blue")) +
    theme_classic(base_size=20, base_family="Helvetica") +
    theme(plot.title=element_text(hjust=0.5, size=title_size),
      plot.margin=unit(c(0.5,2,0.5,0.5), "lines"), legend.position="right")

  return(p)
}
#--------------------------------------------------------------------------------
#
#--------------------------------------------------------------------------------
plot_tsne_signExpression_seurat <- function(
  obj, gene_name, title, title_size, xlabel, ylabel
){
  #--------------------------------------------------
  # Definitions
  #--------------------------------------------------
  mat <- as.matrix(obj@assays[["RNA"]]@data)
  mat <- mat[which(rownames(mat) == gene_name),]
  df <- data.frame(
    x = obj@reductions[["tsne"]]@cell.embeddings[,1],
    y = obj@reductions[["tsne"]]@cell.embeddings[,2],
    expr = mat
  )
  stdev <- sd(df$expr)
  df$expr <- (df$expr - mean(df$expr)) / stdev
  colnames(df) <- c("x", "y", "expr")
  #--------------------------------------------------
  # Plot
  #--------------------------------------------------
  p <- ggplot() +
    geom_point(aes(x=df[,1], y=df[,2], color=df[,3]), size=0.5, alpha=1.0) +
    labs(title=title, x=xlabel, y=ylabel, color="Sign score") +
    scale_colour_gradientn(colours=c("gray90","blue")) +
    theme_classic(base_size=20, base_family="Helvetica") +
    theme(plot.title=element_text(hjust=0.5, size=title_size),
      plot.margin=unit(c(0.5,2,0.5,0.5), "lines"), legend.position="right")

  return(p)
}
#--------------------------------------------------------------------------------
#
#--------------------------------------------------------------------------------
plot_umap_signExpression_seurat <- function(
  obj, gene_name, title, title_size, xlabel, ylabel
){
  #--------------------------------------------------
  # Definitions
  #--------------------------------------------------
  mat <- as.matrix(obj@assays[["RNA"]]@data)
  mat <- mat[which(rownames(mat) == gene_name),]
  df <- data.frame(
    x = obj@reductions[["umap"]]@cell.embeddings[,1],
    y = obj@reductions[["umap"]]@cell.embeddings[,2],
    expr = mat
  )
  stdev <- sd(df$expr)
  df$expr <- (df$expr - mean(df$expr)) / stdev
  colnames(df) <- c("x", "y", "expr")
  #--------------------------------------------------
  # Plot
  #--------------------------------------------------
  p <- ggplot() +
    geom_point(aes(x=df[,1], y=df[,2], color=df[,3]), size=0.5, alpha=1.0) +
    labs(title=title, x=xlabel, y=ylabel, color="Sign score") +
    scale_colour_gradientn(colours=c("gray90","blue")) +
    theme_classic(base_size=20, base_family="Helvetica") +
    theme(plot.title=element_text(hjust=0.5, size=title_size),
      plot.margin=unit(c(0.5,2,0.5,0.5), "lines"), legend.position="right")

  return(p)
}
#--------------------------------------------------------------------------------
#
#--------------------------------------------------------------------------------
plot_GObarplot_seurat <- function(
  obj, cluster, category, showCategory, title, title_size
){
  #--------------------------------------------------
  # Definitions
  #--------------------------------------------------
  library(clusterProfiler)
  data <- obj@misc[["enrichGO"]][[as.character(cluster)]][[category]]
  #--------------------------------------------------
  # Plot
  #--------------------------------------------------
  p <- barplot(data, showCategory = showCategory) +
    ggtitle(title) +
    theme_classic(base_size=20, base_family="Helvetica") +
    theme(plot.title=element_text(hjust=0.5, size=title_size),
      plot.margin=unit(c(0.5,2,0.5,0.5), "lines"), legend.position="right")

  return(p)
}
#--------------------------------------------------------------------------------
#
#--------------------------------------------------------------------------------
plot_GOemapplot_seurat <- function(obj, cluster, category, title){
  #--------------------------------------------------
  # Definitions
  #--------------------------------------------------
  library(clusterProfiler)
  library(enrichplot)
  data <- obj@misc[["enrichGO"]][[as.character(cluster)]][[category]]
  #--------------------------------------------------
  # Plot
  #--------------------------------------------------
  p <- emapplot(pairwise_termsim(data)) +
    ggtitle(title) +
    theme(
      plot.title=element_text(size=20, family="Helvetica", hjust=0.5),
      plot.margin=unit(c(0.5,2,0.5,0.5), "lines"), legend.position="right",
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )

  return(p)
}
