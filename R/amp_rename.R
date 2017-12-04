#' Tidy up taxonomy
#'
#' Used internally in other ampvis functions.
#'
#' @usage amp_rename(data)
#'
#' @param data (\emph{required}) Data list as loaded with \code{\link{amp_load}}.
#' @param tax_empty How to show OTUs without taxonomic information. One of the following:
#' \itemize{
#'    \item \code{"remove"}: Remove OTUs without taxonomic information.
#'    \item \code{"best"}: (\emph{default}) Use the best classification possible. 
#'    \item \code{"OTU"}: Display the OTU name.
#'    }
#' @param tax_class Converts a specific phylum to class level instead, e.g. \code{"p__Proteobacteria"}.
#' @param tax_level The taxonomic level to remove OTUs with no assigned taxonomy, only used when \code{tax_empty = "remove"}. (\emph{default:} \code{"Genus"})
#' 
#' @return A list with 3 dataframes (4 if reference sequences are provided).
#' 
#' @import dplyr
#' @export
#' 
#' @author Kasper Skytte Andersen \email{kasperskytteandersen@@gmail.com}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

amp_rename <- function(data, tax_class = NULL, tax_empty = "best", tax_level = "Genus", tax_string_type = "greengenes"){
tax_string <- c()
  if (tax_string_type == "greengenes"){
    tax_string <- c(Kingdom = "k__", Phylum = "p__", Class = "c__", Order = "o__", Family = "f__", Genus = "g__", Species = "s__")
  } else if (tax_string_type == "silva"){
    tax_string <- c(Kingdom = "D_0__", Phylum = "D_1__", Class = "D_2__", Order = "D_3__", Family = "D_4__", Genus = "D_5__", Species = "D_6__")
  }
  tax = data[["tax"]]
  
  ## First make sure that all entries are strings
  for ( i in 1:ncol(tax) ){
    tax[,i] <- as.character(tax[,i])  
  }
  
  ## Change a specific phylum to class level
  if(!is.null(tax_class)){
    for (i in 1:nrow(tax)){
      if (!is.na(tax$Phylum[i]) & tax$Phylum[i] %in% tax_class){
        tax$Phylum[i] <- tax$Class[i]   
      }
    }
  }
  
  ## Remove the underscore classifier from the data  
  tax$Kingdom <- gsub(sprintf("^%s*", tax_string[["Kingdom"]]), "", tax$Kingdom)
  tax$Phylum <- gsub(sprintf("^%s*", tax_string[["Phylum"]]), "", tax$Phylum)
  tax$Class <- gsub(sprintf("^%s*", tax_string[["Class"]]), "", tax$Class)
  tax$Order <- gsub(sprintf("^%s*", tax_string[["Order"]]), "", tax$Order)
  tax$Family <- gsub(sprintf("^%s*", tax_string[["Family"]]), "", tax$Family)
  tax$Genus <- gsub(sprintf("^%s*", tax_string[["Genus"]]), "", tax$Genus)
  tax$Kingdom <- gsub("uncultured", "", tax$Kingdom)
  tax$Phylum <- gsub("uncultured", "", tax$Phylum)
  tax$Class <- gsub("uncultured", "", tax$Class)
  tax$Order <- gsub("uncultured", "", tax$Order)
  tax$Family <- gsub("uncultured", "", tax$Family)
  tax$Genus <- gsub("uncultured", "", tax$Genus)
  
  ## Check if there is a species level otherwise add it for consistency
  if (!is.null(tax$Species)){
    tax$Species <- gsub(sprintf("^%s*", tax_string[["Species"]]), "", tax$Species)
  } else {
    tax$Species <- ""
  }
  
  tax[is.na(tax)] <- ""
  
  ## How to handle empty taxonomic assignments
  if (tax_empty == "OTU"){
    for (i in 1:nrow(tax)) {
      if (tax[i,"Species"] == "") {tax[i,"Species"] <- rownames(tax)[i]}
      if (tax[i,"Genus"] == "") {tax[i,"Genus"] <- rownames(tax)[i]}
      if (tax[i,"Family"] == "") {tax[i,"Family"] <- rownames(tax)[i]}
      if (tax[i,"Order"] == "") {tax[i,"Order"] <- rownames(tax)[i]}
      if (tax[i,"Class"] == "") {tax[i,"Class"] <- rownames(tax)[i]}
      if (tax[i,"Phylum"] == "") {tax[i,"Phylum"] <- rownames(tax)[i]}
    }
  }
  
  ## Handle empty taxonomic strings
  rn <- rownames(tax) #damn rownames are silently dropped by mutate()
  if(tax_empty == "best"){
    tax <- mutate(tax, Kingdom, Kingdom = ifelse(Kingdom == "", "Unclassified", Kingdom)) %>%
      mutate(Phylum, Phylum = ifelse(Phylum == "", paste(sprintf("%s", tax_string[["Kingdom"]]), Kingdom, "_", rownames(tax), sep = ""), Phylum)) %>%
      mutate(Class, Class = ifelse(Class == "", ifelse(grepl("__", Phylum), Phylum, paste(sprintf("%s", tax_string[["Phylum"]]), Phylum, "_", rownames(tax), sep = "")), Class)) %>%
      mutate(Order, Order = ifelse(Order == "", ifelse(grepl("__", Class), Class, paste(sprintf("%s", tax_string[["Class"]]), Class, "_", rownames(tax), sep = "")), Order)) %>%
      mutate(Family, Family = ifelse(Family == "", ifelse(grepl("__", Order), Order, paste(sprintf("%s", tax_string[["Order"]]), Order, "_", rownames(tax), sep = "")), Family)) %>%
      mutate(Genus, Genus = ifelse(Genus == "", ifelse(grepl("__", Family), Family, paste(sprintf("%s", tax_string[["Family"]]), Family, "_", rownames(tax), sep = "")), Genus)) %>%
      mutate(Species, Species = ifelse(Species == "", ifelse(grepl("__", Genus), Genus, paste(sprintf("%s", tax_string[["Genus"]]), Genus, "_", rownames(tax), sep = "")), Species))
  }
  rownames(tax) <- rn
  
  if(tax_empty == "remove"){
    abund <- data[["abund"]]
    tax <- subset(tax, tax[,tax_level] != "")
    abund <- subset(abund, rownames(abund) %in% rownames(tax))
    data[["abund"]] <- abund
  }
  data[["tax"]] <- tax
  rownames(data[["tax"]]) <- rownames(tax)
  
  return(data)
}
