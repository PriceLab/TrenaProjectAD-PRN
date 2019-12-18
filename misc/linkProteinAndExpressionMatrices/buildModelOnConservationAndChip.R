library(TrenaValidator)
library(TrenaProjectHG38.generic)
library(igvR)
library(GenomicScores)
library(phastCons7way.UCSC.hg38); phast.7 <- phastCons7way.UCSC.hg38
library(TrenaProjectAD)
library(trenaSGM)
targetGene <- "MEF2C"
if(!exists("tp")){
   tp <- TrenaProjectAD();
   setTargetGene(tp, targetGene)
   tbl.geneInfo <- getTranscriptsTable(tp)
   tbl.benchmark <- get(load(system.file(package="TrenaValidator", "extdata", "tbl.A.RData")))
   tbl.benchmark$pubmed.count <- unlist(lapply(strsplit(tbl.benchmark$pubmedID_from_curated_resources, ","), length))
   tv <- TrenaValidator(TF="TWIST1", targetGene, tbl.benchmark);
   }

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
if(!exists("igv")){
   igv <- igvR()
   setGenome(igv, "hg38")
   tbl.regions <- with(tbl.geneInfo, data.frame(chrom=chrom, start=tss-5000, end=tss+5000, stringsAsFactors=FALSE))
   with(tbl.regions, showGenomicRegion(igv, sprintf("%s:%d-%d", chrom, start, end)))
   }
#------------------------------------------------------------------------------------------------------------------------
geneHancerTrack <- function()
{
   tbl.gh <- getEnhancers(tp)
   tbl.gh$width <- with(tbl.gh, 1 + end - start)
   #tbl.gh <- subset(tbl.gh, width < 5000)
   track <- DataFrameQuantitativeTrack("gh", tbl.gh[, c(1,2,3,11)], color="blue", autoscale=FALSE, min=0, max=50)
   displayTrack(igv, track)
   return(tbl.gh)

} # geneHancerTrack
#------------------------------------------------------------------------------------------------------------------------
conservationTrack <- function(tbl.regions=NA)
{
   loc <- getGenomicRegion(igv)
   starts <- with(loc, seq(start, end, by=5))
   ends <- starts + 5
   count <- length(starts)
   tbl.blocks <- data.frame(chrom=rep(loc$chrom, count), start=starts, end=ends, stringsAsFactors=FALSE)

   if(!is.na(tbl.regions)){
      tbl.overlaps <- as.data.frame(findOverlaps(GRanges(tbl.blocks), GRanges(tbl.regions)))
      tbl.blocks <- tbl.blocks[tbl.overlaps[,1],]
      }

   tbl.cons7 <- as.data.frame(gscores(phast.7, GRanges(tbl.blocks)), stringsAsFactors=FALSE)
   tbl.cons7$chrom <- as.character(tbl.cons7$seqnames)
   tbl.cons7 <- tbl.cons7[, c("chrom", "start", "end", "default")]
   track <- DataFrameQuantitativeTrack("phast7", tbl.cons7, autoscale=TRUE, color="red")
   displayTrack(igv, track)

   return(tbl.cons7)

} # conservationTrack
#------------------------------------------------------------------------------------------------------------------------
chipTrack <- function(tf, tissueRestriction=NA_character_, peaksAlso=FALSE)
{
   chrom.loc <- getGenomicRegion(igv)
   tbl.chip <- with(chrom.loc, getChipSeq(tp.hg38, chrom, start, end, tf))

   if(nrow(tbl.chip) == 0){
      printf("--- no ChIP found for %s", tf)
      return()
      }

   if(!is.na(tissueRestriction))
      tbl.chip <- tbl.chip[grep(tissueRestriction, tbl.chip$name),]

   if(nrow(tbl.chip) == 0){
      printf("--- no ChIP found for %s", tf)
      return()
      }

   printf("tf %s,  %d peaks", tf,  nrow(tbl.chip))

   if(nrow(tbl.chip) == 0){
      #printf("no ChIP for TF %s in %d bases", tf, chrom.loc$end - chrom.loc$start)
      return(data.frame())
      }
   tbl.track <- tbl.chip[, c("chrom", "start", "endpos", "name")]
   trackName <- sprintf("Ch-%s", tf)
   track <- DataFrameAnnotationTrack(trackName, tbl.track, color="random", displayMode="squished", trackHeight=25)
   displayTrack(igv, track)

   if(peaksAlso){
      tbl.track <- tbl.chip[, c("chrom", "peakStart", "peakEnd", "name")]
      trackName <- sprintf("peaks-%s", tf)
      track <- DataFrameAnnotationTrack(trackName, tbl.track, color="random", displayMode="squished", trackHeight=25)
      displayTrack(igv, track)
      } # peaksAlso

} # chipTrack
#------------------------------------------------------------------------------------------------------------------------
run.mef2c <- function()
{
   tbl.gh <- geneHancerTrack()
   tbl.cons <- conservationTrack(tbl.gh)
   motifs <- query(MotifDb, "hsapiens", c("jaspar2018"))
   meme.file <- "human.jaspar2018.meme"
   export(motifs, con=meme.file, format="meme")
   tbl.tfbs <- getTFBS.fimo(tv, tbl.gh, fimo.threshold=1e-4, conservation.threshold=0.9, meme.file)
   candidate.tfs <- unique(tbl.tfbs$tf)

   genome <- "hg38"

   mtx <- mtx.rna
   mtx["MEF2C",] <- mtx.prot["MEF2C",]

   candidate.tfs <- intersect(candidate.tfs, rownames(mtx))
   length(candidate.tfs)
   build.spec <- list(title="mef2c.noDNA.conserved.gh.tfs",
                      type="noDNA.tfsSupplied",
                      matrix=mtx,
                      candidateTFs=candidate.tfs,
                      tfPool=allKnownTFs(),
                      tfPrefilterCorrelation=0.0,
                      annotationDbFile=dbfile(org.Hs.eg.db),
                      orderModelByColumn="rfScore",
                      solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman", "xgboost"),
                      quiet=TRUE)

   builder <- NoDnaModelBuilder(genome, targetGene,  build.spec, quiet=TRUE)
   x <- build(builder)


} # run.mef2c
#------------------------------------------------------------------------------------------------------------------------
