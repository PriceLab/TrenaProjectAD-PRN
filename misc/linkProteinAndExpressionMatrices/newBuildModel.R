library(TrenaProjectAD)
library(org.Hs.eg.db)
library(trenaSGM)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("tp"))
   tp <- TrenaProjectAD();

if(!exists("tbl.idMap"))
   tbl.idMap <- get(load("~/github/TrenaProjectADPRN/misc/linkProteinAndExpressionMatrices/tbl.rna.prot.map.RData"))

if(!exists("mtx.prot")){
   mtx.prot <- get(load("~/github/TrenaProjectADPRN/inst/extdata/expression/c2-8817x400.RData"))
   mtx.prot <- mtx.prot[, tbl.idMap$proteinID]
   mtx.prot[is.na(mtx.prot)] <- 0
   dim(mtx.prot)
   colnames(mtx.prot) <- tbl.idMap$individualID
   }

if(!exists("mtx.rna")){
   mtx.rna <- getExpressionMatrix(tp, "rosmap.14235x632")
   mtx.rna <- mtx.rna[, tbl.idMap$mrna]
   mtx.rna[is.na(mtx.rna)] <- 0
   dim(mtx.rna)
   colnames(mtx.rna) <- tbl.idMap$individualID
   }


#------------------------------------------------------------------------------------------------------------------------
mef2c.model <- function(display=FALSE)
{
   genome <- "hg38"
   targetGene <- "MEF2C"
   setTargetGene(tp, targetGene)
   tbl.info <- getTranscriptsTable(tp)
   tss <- tbl.info$tss
   chrom <- tbl.info$chrom

      # get genehancer regulatory regions, all tissues, then demonstrate
      # how further filtering can be done if you wish to do so

   tbl.regulatoryRegions <- getGeneRegulatoryRegions(tp, tissues="all")
   checkTrue(nrow(tbl.regulatoryRegions) > 8)

   getExpressionMatrixNames(tp)
   mtx <- mtx.rna
   mtx["MEF2C",] <- mtx.prot["MEF2C",]


   recipe <- list(title="gh all tissue demo on MEF2C",
                  type="footprint.database",
                  regions=tbl.regions,
                  geneSymbol=targetGene,
                  tss=tss,
                  matrix=mtx,
                  db.host="khaleesi.systemsbiology.net",
                  db.port=5432,
                  databases=list("brain_hint_20"),
                  annotationDbFile=dbfile(org.Hs.eg.db),
                  motifDiscovery="builtinFimo",
                  tfPool=allKnownTFs(),
                  tfMapping="MotifDB",
                  tfPrefilterCorrelation=0.1,
                  orderModelByColumn="rfScore",
                  solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"))

   fpBuilder <- FootprintDatabaseModelBuilder(genome, targetGene, recipe, quiet=FALSE)

   msg <- sprintf("gh all tissues demo on TREM2, fp model %d regions", nrow(tbl.regulatoryRegions))
   x <- build(fpBuilder)

   checkTrue("SPI1" %in% head(x.trem2$model$gene))
   if(display){
      browser()
      require(igvR)
      igv <- igvR()
      setGenome(igv, "hg38")
      chrom <- tbl.info$chrom
      shoulder <- 1000
      start <- min(tbl.regulatoryRegions$start) - shoulder
      end   <- max(tbl.regulatoryRegions$end) + shoulder
      showGenomicRegion(igv, sprintf("%s:%d-%d", chrom, start, end))
      track <- DataFrameQuantitativeTrack("gh", tbl.regulatoryRegions[, c("chrom", "start", "end", "combinedscore")],
                                          color="blue", autoscale=TRUE)
      displayTrack(igv, track)
      }

} # mef2c.model
#------------------------------------------------------------------------------------------------------------------------
