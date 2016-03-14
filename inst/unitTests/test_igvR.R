library(igvR)
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#------------------------------------------------------------------------------------------------------------------------
if(!exists("igv"))
    igv <- igvR()
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_connected()
   clearAllTracks(igv)
   test_sampleBedFile()
   test_sendRawCommand()
   test_trackVisibility()
   test_goto()
   test_displayBedTable()
   test_displayGWASTable()
   test_displayScoredFeatures()
   test_displayVcfRegion()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_constructorArgs <- function()
{
   printf("--- test_constructorArgs")
   disconnect(igv)
   vcfDirectory <- system.file(package="igvR", "inst", "extdata", "vcf")
   checkTrue(file.axists(vcfDirectory))
   genome <- "hg38"

   igv <<- igvR(host="localhost", port=60151, genome=genome, vcfDirectory)
   checkTrue(connected(igv))
   checkEquals(getVcfDirectory(igv), vcfDirectory);
   checkEquals(getGenome(igv), genome)



} # test_connected
#------------------------------------------------------------------------------------------------------------------------
test_connected <- function()
{
    printf("--- test_connected")
    checkTrue(connected(igv))


} # test_connected
#------------------------------------------------------------------------------------------------------------------------
# the sample file is lifted directly from
# https://genome.ucsc.edu/FAQ/FAQformat.html#format1
# chr22 1000 5000 cloneA 960 + 1000 5000 0 2 567,488, 0,3512
# chr22 2000 6000 cloneB 900 - 2000 6000 0 2 433,399, 0,3601
test_sampleBedFile = function ()
{
  printf("--- test_sampleBedFile");

  if(!exists("igv"))
      igv <- igvR()

  full.path <- system.file(package="igvR", "extdata", "sample.bed")
  result <- loadFile(igv, full.path)
  checkEquals(result, "OK\n");

} # test.tinyBedFile
#------------------------------------------------------------------------------------------------------------------------
test_goto <- function()
{
  printf("--- test_goto")
  if(!exists("igv"))
      igv <- igvR()

   goto(igv, "chr22", 800, 6200)

} # test_goto
#------------------------------------------------------------------------------------------------------------------------
test_displayBedTable <- function()
{
   printf("--- test_displayBedTable")
   chr <- rep("chr1", 3)
   start <- rep(100470155, 3)
   end <- rep(100470168, 3)
   name <- c("MA0679.1", "MA0756.1", "MA0757.1")
   name <- rep(" ", 3)
   score <- c(56.0, 58.6, 60.5)
   strand <- rep("+", 3)
   tbl <- data.frame(chr=chr, start=start, end=end, name=name, score=score, strand=strand, stringsAsFactors=FALSE)

   if(!exists("igv"))
      igv <- igvR()

   result <- displayBedTable(igv, tbl, name="X3")
   checkEquals(result, "OK\n")
   result <- goto(igv, "chr1", start[1] - 20, end[3] + 20)
   checkEquals(result, "OK\n")

} # test_displayBedTable
#------------------------------------------------------------------------------------------------------------------------
test_displayGWASTable <- function()
{
   printf("--- test_displayGWASTable")
   count <- 10

    # GWAS file must contain four columns (case-insensitive):
    #   CHR: chromosome (aliases chr, chromosome)
    #    BP: nucleotide location (aliases bp, pos, position)
    #   SNP: SNP identifier (aliases snp, rs, rsid, rsnum, id, marker, markername)
    #     P: p-value for the association (aliases p, pval, p-value, pvalue, p.value)

   chr <- rep("chr1", count)
   bp <- seq(from=1000, by=10, length.out=10)
   snp <- paste("rs", 1:10, sep="")
   exponents <- as.integer(runif(count, 3, 40))
   base <- runif(count, 1, 10)
   p <- unlist(lapply(1:count, function(i) base[i]^-exponents[i]))

   tbl <- data.frame(CHR=chr, BP=bp, SNP=snp, P=p, stringsAsFactors=FALSE)
   if(!exists("igv"))
      igv <- igvR()


   result <- displayGWASTable(igv, tbl, name="10 gwas snps")
   checkEquals(result, "OK\n")
   result <- goto(igv, "chr1", bp[1] - 20, bp[count] + 20)
   checkEquals(result, "OK\n")

} # test_displayGWASTable
#------------------------------------------------------------------------------------------------------------------------
test_displayScoredFeatures <- function()
{
   printf("--- test_displayScoredFeatures")
   count <- 8
   chr <- rep("chr5", count)
   start <- c(88364636, 88364636, 88364787, 88365143, 88365146, 88365146, 88365146, 88365168)
   end <-   c(88364651, 88364651, 88364821, 88365165, 88365165, 88365165, 88365165, 88365187)
   feature.names <- sprintf("FP.%02d", 1:count)
   score <- as.integer(runif(count, 1, 100))
   tbl <- data.frame(chr=chr, start=start, end=end, featureName=feature.names, test_fp=score)

   if(!exists("igv"))
      igv <- igvR()

   result <- displayScoredFeatures(igv, tbl)
   goto(igv, "chr5", start[1] - 200, end[count] + 200)

} # test_displayScoredFeatures
#------------------------------------------------------------------------------------------------------------------------
test_displayVcfRegion <- function()
{
   printf("--- test_displayVcfRegion")

   vcfDirectory <- system.file(package="igvR", "extdata", "vcf")
   vcfFile <- file.path(vcfDirectory, "chr22-sub.vcf.gz")
   checkTrue(file.exists(vcfFile))
   checkTrue(file.exists(sprintf("%s.tbi", vcfFile)))

   if(!exists("igv"))
       igv <- igvR()

   checkTrue(connected(igv))

      # a 20kb region
   start <- 49901760
   end   <- 49921796

   group.1 <- c("HG00096", "HG00101")
   group.2 <- c("HG00097", "HG00099", "HG00100")

      # display all samples in 1 track, then two more tracks with 2 and 3 samples respectively
   displayVcfRegion(igv, "chr22", start, end, vcfFile)
   displayVcfRegion(igv, "chr22", start, end, vcfFile, sampleIDs=group.1)
   displayVcfRegion(igv, "chr22", start, end, vcfFile, sampleIDs=group.2)

   goto(igv, "chr22", start, end)

} # test_displayVcfRegion
#------------------------------------------------------------------------------------------------------------------------
test_sendRawCommand <- function()
{
   printf("--- test_sendRawCommand")
   for(i in 1:3){
      sendRawCommand(igv, "goto chr1:169,436,527-169,802,250")
      sendRawCommand(igv, "goto chr1:169,510,853-169,876,576")
      } # for i

} # test_sendRawCommand
#------------------------------------------------------------------------------------------------------------------------
test_trackVisibility <- function()
{
   printf("--- test_trackVisibility")
   start <- c(100, 100, 100, 200, 300, 400)
   end <-   c(150, 150, 176, 220, 380, 500)
   count <- length(start)
   strand <- rep("+", count)
   chr <- rep("chr1", count)

   feature.names <- sprintf("FP.%02d", 1:count)
   scores <- as.integer(runif(count, 1, 100))

   tbl <- data.frame(chr=chr, start=start, end=end, name=feature.names, score=scores, strand=strand, stringsAsFactors=FALSE)
   displayBedTable(igv, tbl, "foo")
   goto(igv, "chr1", 50, 550)

   expandTrack(igv, "foo.bed")
   squishTrack(igv, "foo.bed")
   collapseTrack(igv, "foo.bed")

} # test_trackVisibility
#------------------------------------------------------------------------------------------------------------------------
