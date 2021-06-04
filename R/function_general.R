#--------------------------------------------------------------------------------
#
#--------------------------------------------------------------------------------
make_asurat_obj <- function(mat, obj_name){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  obj <- list(
    history = c(),
    variable = data.frame(
      symbol = rownames(mat),
      entrez = rep(NA, dim(mat)[1])
    ),
    sample = data.frame(
      barcode = colnames(mat)
    ),
    data = c(),
    misc = c()
  )
  obj[["data"]][["raw"]] <- as.data.frame(mat)
  #--------------------------------------------------
  # History
  #--------------------------------------------------
  slot_name <- "make_asurat_obj"
  obj[["history"]][[slot_name]][["obj_name"]] <- obj_name

  return(obj)
}
#--------------------------------------------------------------------------------
#
#--------------------------------------------------------------------------------
do_quickQC <- function(obj, min_nsamples, mitochondria_symbol){
  #--------------------------------------------------
  # History
  #--------------------------------------------------
  slot_name <- "do_quickQC"
  obj[["history"]][[slot_name]][["min_nsamples"]] <- min_nsamples
  #--------------------------------------------------
  # Reduce variables
  #--------------------------------------------------
  tmp <- as.matrix(obj[["data"]][["raw"]])
  inds <- which(apply(tmp, 1, function(x) sum(x>0)) >= min_nsamples)
  mat <- tmp[inds,]
  obj[["data"]][["raw"]] <- as.data.frame(mat)
  #--------------------------------------------------
  # Other information
  #--------------------------------------------------
  obj[["variable"]] <- data.frame(
    symbol = obj[["variable"]][["symbol"]][inds],
    entrez = obj[["variable"]][["entrez"]][inds]
  )
  obj[["sample"]] <- data.frame(
    barcode = colnames(mat),
    nReads = as.integer(apply(mat, 2, function(x) sum(x))),
    nGenes = as.integer(apply(mat, 2, function(x) sum(x > 0))),
    percent_MT = apply(
      mat[grepl(mitochondria_symbol, rownames(mat)),], 2, function(x) 100 * sum(x)
    ) / apply(mat, 2, function(x) sum(x))
  )

  return(obj)
}
#--------------------------------------------------------------------------------
#
#--------------------------------------------------------------------------------
trim_samples <- function(obj,
  min_nReads = 0, max_nReads = 1e+10,
  min_nGenes = 0, max_nGenes = 1e+10,
  min_percent_MT = 0, max_percent_MT = 100
){
  #--------------------------------------------------
  # History
  #--------------------------------------------------
  slot_name <- "trim_samples"
  obj[["history"]][[slot_name]][["min_nReads"]] <- min_nReads
  obj[["history"]][[slot_name]][["max_nReads"]] <- max_nReads
  obj[["history"]][[slot_name]][["min_nGenes"]] <- min_nGenes
  obj[["history"]][[slot_name]][["max_nGenes"]] <- max_nGenes
  obj[["history"]][[slot_name]][["min_percent_MT"]] <- min_percent_MT
  obj[["history"]][[slot_name]][["max_percent_MT"]] <- max_percent_MT
  #--------------------------------------------------
  # Reduce count matrices
  #--------------------------------------------------
  inds_1 <- which(
    (obj[["sample"]][["nReads"]] >= min_nReads) &
    (obj[["sample"]][["nReads"]] <= max_nReads)
  )
  inds_2 <- which(
    (obj[["sample"]][["nGenes"]] >= min_nGenes) &
    (obj[["sample"]][["nGenes"]] <= max_nGenes)
  )
  inds_3 <- which(
    (obj[["sample"]][["percent_MT"]] >= min_percent_MT) &
    (obj[["sample"]][["percent_MT"]] <= max_percent_MT)
  )
  inds <- intersect(intersect(inds_1, inds_2), inds_3)
  obj[["data"]][["raw"]] <- as.data.frame(obj[["data"]][["raw"]][,inds])
  #--------------------------------------------------
  # The other information
  #--------------------------------------------------
  obj[["sample"]] <- data.frame(
    barcode = obj[["sample"]][["barcode"]][inds],
    nReads = as.integer(obj[["sample"]][["nReads"]][inds]),
    nGenes = as.integer(obj[["sample"]][["nGenes"]][inds]),
    percent_MT = obj[["sample"]][["percent_MT"]][inds]
  )

  return(obj)
}
#--------------------------------------------------------------------------------
#
#--------------------------------------------------------------------------------
trim_variables <- function(obj, min_meanReads){
  #--------------------------------------------------
  # History
  #--------------------------------------------------
  slot_name <- "trim_variables"
  obj[["history"]][[slot_name]][["min_meanReads"]] <- min_meanReads
  #--------------------------------------------------
  # Reduce count matrix
  #--------------------------------------------------
  tmp <- as.matrix(obj[["data"]][["raw"]])
  inds <- which(apply(tmp, 1, mean) >= min_meanReads)
  mat <- tmp[inds, ]
  obj[["data"]][["raw"]] <- as.data.frame(mat)
  #--------------------------------------------------
  # The others
  #--------------------------------------------------
  genes <- data.frame(
    symbol  = obj[["variable"]][["symbol"]][inds],
    entrez  = obj[["variable"]][["entrez"]][inds]
  )
  obj[["variable"]] <- genes

  return(obj)
}

