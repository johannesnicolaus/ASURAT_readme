#--------------------------------------------------------------------------------
#
#--------------------------------------------------------------------------------
collect_Reactome <- function(organism, category_names){
  #--------------------------------------------------
  # Definitions
  #--------------------------------------------------
  orgname <- paste(organism, ": ", sep = "")
  genes <- c()
  for(category in category_names){
    if(category == "PATHID"){
      ids <- toTable(reactomePATHID2EXTID)
    }
    genes <- c(genes, ids$gene_id)
  }
  genes <- unique(genes)
  #--------------------------------------------------
  # collect_Reactome
  #--------------------------------------------------
  res <- list()
  for(category in category_names){
    tmp <- c()
    if(category == "PATHID"){
      tmp <- toTable(reactomePATHNAME2ID)
      tmp <- tmp[grep(orgname, iconv(tmp$path_name)), ]  # `tmp` is a data.frame
      ids <- toTable(reactomePATHID2EXTID)
    }
    colnames(tmp) <- c("ID", "Description")
    #--------------------------------------------------
    # Prepare map
    #--------------------------------------------------
    failure <- c()
    map <- data.frame(
      matrix(ncol = 3, nrow = 0, dimnames = list(NULL, c(
        "ID",
        "Description",
        "NCBI_geneID"
      ))))
    I <- length(tmp$ID)
    for(i in 1:I){
      #--------------------------------------------------
      # Print the current process
      #--------------------------------------------------
      if((i == 1) || (i %% floor(0.25 * I) == 0 && i < 0.95 * I) || (i == I)){
        text <- paste("Now processing ", i, "/", I, " for ", category, "...\n", sep = "")
        cat(text)
      }
      #------------------------------
      # Note: tmp and ids are not always comparable
      #------------------------------
      if(length(which(ids$DB_ID == tmp$ID[i])) == 0){
        failure <- rbind(failure, tmp$ID[i])
        next
      }
      gene_ncbiids <- as.matrix(ids[which(ids$DB_ID == tmp$ID[i]),]$gene_id)
      add <- data.frame(
        ID = tmp$ID[i],
        Description = tmp$Description[i],
        NCBI_geneID = gene_ncbiids
      )
      map <- rbind(map, add)
    }
    res[[category]][["success"]] <- map
    res[[category]][["failure"]] <- data.frame(Failed_DB_ID = failure)
  }

  return(res)
}
#--------------------------------------------------------------------------------
#
#--------------------------------------------------------------------------------
format_Reactome <- function(dict){
  #--------------------------------------------------
  # Definitions
  #--------------------------------------------------
  category_names <- names(dict)

  #--------------------------------------------------
  # Reformat
  #--------------------------------------------------
  res <- list()
  for(category in category_names){
    tmp <- dict[[category]]
    map <- unique(data.frame(
      ID = tmp[["ID"]],
      Description = tmp[["Description"]],
      Count = NA,
      Gene = NA,
      GeneID = NA
    ))
    for(i in 1:nrow(map)){
      #------------------------------
      # Gene and Count
      #------------------------------
      genes <- unique(tmp[which(tmp[["ID"]] == map[["ID"]][i]),]$NCBI_geneID)
      map$GeneID[i] <- paste(genes, collapse = "/")
      map$Count[i] <- length(genes)
    }
    rownames(map) <- 1:nrow(map)
    res[[category]] <- map
  }

  return(res)
}
