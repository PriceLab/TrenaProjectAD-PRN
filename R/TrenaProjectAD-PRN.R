#---------------------------------------------------------------------------------------------------
#' @import methods
#' @import TrenaProject
#' @importFrom AnnotationDbi select
#' @import org.Hs.eg.db
#'
#' @title TrenaProjectAD-PRN-class
#'
#' @name TrenaProjectAD-PRN-class
#' @rdname TrenaProjectAD-PRN-class
#' @aliases TrenaProjectAD-PRN
#' @exportClass TrenaProjectAD-PRN
#'

.TrenaProjectAD-PRN <- setClass("TrenaProjectAD-PRN",
                                  contains="TrenaProjectHG38")

#----------------------------------------------------------------------------------------------------
#' Define an object of class TrenaProjectAD-PRN
#'
#' @description
#' Expression, variant and covariate data for the genes of interest (perhaps unbounded) for pre-term birth studies
#'
#' @rdname TrenaProjectAD-PRN-class
#'
#' @export
#'
#' @return An object of the TrenaProjectAD-PRN class
#'

TrenaProjectAD-PRN <- function(quiet=TRUE)

{
   genomeName <- "hg38"

   directory <- system.file(package="TrenaProjectAD-PRN", "extdata", "geneSets")
   geneSet.files <- list.files(directory)
   geneSets <- list()

   for(file in geneSet.files){
      full.path <- file.path(directory, file)
      genes <- scan(full.path, sep="\t", what=character(0), quiet=TRUE, comment.char="#")
      geneSet.name <- sub(".txt", "", file)
      geneSets[[geneSet.name]] <- genes
      }

   footprintDatabaseNames <- c("brain_hint_16",  "brain_hint_20", "brain_wellington_16", "brain_wellington_20")
   footprintDatabaseHost <- "khaleesi.systemsbiology.net"
   footprintDatabasePort <- 5432
   dataDirectory <- system.file(package="TrenaProjectAD-PRN", "extdata")

   covariatesFile <- NA_character_;

   .TrenaProjectAD-PRN(TrenaProjectHG38(projectName="TrenaProjectAD-PRN",
                                       supportedGenes=geneSets[[1]],
                                       footprintDatabaseHost=footprintDatabaseHost,
                                       footprintDatabasePort=footprintDatabasePort,
                                       footprintDatabaseNames=footprintDatabaseNames,
                                       packageDataDirectory=dataDirectory,
                                       quiet=quiet
                                       ))

} # TrenaProjectAD-PRN, the constructor
#----------------------------------------------------------------------------------------------------
