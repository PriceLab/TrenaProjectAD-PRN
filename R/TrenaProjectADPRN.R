#---------------------------------------------------------------------------------------------------
#' @import methods
#' @import TrenaProject
#' @importFrom AnnotationDbi select
#' @import org.Hs.eg.db
#'
#' @title TrenaProjectADPRN-class
#'
#' @name TrenaProjectADPRN-class
#' @rdname TrenaProjectADPRN-class
#' @aliases TrenaProjectADPRN
#' @exportClass TrenaProjectADPRN
#'

.TrenaProjectADPRN <- setClass("TrenaProjectADPRN",
                                  contains="TrenaProjectHG38")

#----------------------------------------------------------------------------------------------------
#' Define an object of class TrenaProjectADPRN
#'
#' @description
#' Expression, variant and covariate data for the genes of interest (perhaps unbounded) for pre-term birth studies
#'
#' @rdname TrenaProjectADPRN-class
#'
#' @export
#'
#' @return An object of the TrenaProjectADPRN class
#'

TrenaProjectADPRN <- function(quiet=TRUE)

{
   genomeName <- "hg38"

   directory <- system.file(package="TrenaProjectADPRN", "extdata", "geneSets")
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
   dataDirectory <- system.file(package="TrenaProjectADPRN", "extdata")

   covariatesFile <- NA_character_;

   .TrenaProjectADPRN(TrenaProjectHG38(projectName="TrenaProjectADPRN",
                                       supportedGenes=geneSets[[1]],
                                       footprintDatabaseHost=footprintDatabaseHost,
                                       footprintDatabasePort=footprintDatabasePort,
                                       footprintDatabaseNames=footprintDatabaseNames,
                                       packageDataDirectory=dataDirectory,
                                       quiet=quiet
                                       ))

} # TrenaProjectADPRN, the constructor
#----------------------------------------------------------------------------------------------------
