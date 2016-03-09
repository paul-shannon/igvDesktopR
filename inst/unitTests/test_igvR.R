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
   score <- c(56.0, 58.6, 60.5)
   strand <- rep("+", 3)
   tbl <- data.frame(chr=chr, start=start, end=end, name=name, score=score, strand=strand, stringsAsFactors=FALSE)

  if(!exists("igv"))
      igv <- igvR()

   result <- displayBedTable(igv, tbl, name="3 motifs")
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
   print("--- test_displayScoredFeatures")
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
private_test_displayVcfRegion <- function()
{
   printf("--- test_displayVcfRegion")

   vcfDirectory <- system.file(package="igvR", "extdata", "vcf")
   checkTrue(file.exists(file.path(vcfDirectory, "chr5-sub.vcf.gz")))
   checkTrue(file.exists(file.path(vcfDirectory, "chr5-sub.vcf.gz.tbi")))

   igv <- igvR()

   checkTrue(connected(igv))

      # a 37.5kb region around the several mef2c transcription start sites
   start <- 88873095
   stop <-  88910640

   apoe.dementia.samples <- c("002_S_1268", "002_S_4521", "003_S_1057", "06_S_4346", "006_S_4546")
   control.samples       <- c("002_S_0413", "002_S_0685", "002_S_1261", "002_S_1280", "002_S_4213")

   displayVcfRegion(igv, "chr5", start, stop, vcfDirectory, sampleIDs=apoe.dementia.samples)
   displayVcfRegion(igv, "chr5", start, stop, vcfDirectory, sampleIDs=control.samples)

} # test_displayVcfRegion
#------------------------------------------------------------------------------------------------------------------------
test_displayVcfRegion <- function()
{
   printf("--- test_displayVcfRegion")

   vcfDirectory <- system.file(package="igvR", "extdata", "vcf")
   checkTrue(file.exists(file.path(vcfDirectory, "chr22-sub.vcf.gz")))
   checkTrue(file.exists(file.path(vcfDirectory, "chr22-sub.vcf.gz.tbi")))

   if(!exists("igv"))
       igv <- igvR()

   checkTrue(connected(igv))

      # a 20kb region
   start <- 49901760
   end   <- 49921796

   group.1 <- c("HG00096", "HG00101")
   group.2 <- c("HG00097", "HG00099", "HG00100")

      # display all samples in 1 track, then two more tracks with 2 and 3 samples respectively
   displayVcfRegion(igv, "chr22", start, end, vcfDirectory)
   displayVcfRegion(igv, "chr22", start, end, vcfDirectory, sampleIDs=group.1)
   displayVcfRegion(igv, "chr22", start, end, vcfDirectory, sampleIDs=group.2)

   goto(igv, "chr22", start, end)

} # test_displayVcfRegion
#------------------------------------------------------------------------------------------------------------------------
# functions stolen from the Biocondcutor sradb package, used as an guide to early development.
# IGVsocket <- function(host='localhost', port=60151)
# {
#    sock = try(make.socket(host,port))
#    if(inherits(sock,'try-error')) {
#      stop(sprintf("Make sure that IGV is running on host %s and that IGV is set to accept commands on port %d",host,port))
#      }
#    print(sock)
#    return(sock)
# }
# #------------------------------------------------------------------------------------------------------------------------
# .socketWrite<- function(sock,string)
# {
#    print(string)
#    write.socket(sock,string)
#    printf(".socketWrite waiting for response");
#    response <- read.socket(sock)
#    printf("       response received: %s", response);
#
#    return(response)
# }
# #------------------------------------------------------------------------------------------------------------------------
# IGVload <-  function (sock, files)
# {
#    if(length(files)<1) stop('Files must be specified')
#    message(basename(files))
#    .socketWrite(sock,paste('load',paste(files,collapse=','),'\n'))
# }
# #------------------------------------------------------------------------------------------------------------------------
# IGVgoto <- function(sock,region)
# {
#    .socketWrite(sock, paste('goto',region,'\n'))
# }
# #------------------------------------------------------------------------------------------------------------------------
# IGVgenome <- function(sock,genome="hg18")
# {
#    .socketWrite(sock,paste('genome',genome,'\n'))
# }
# #------------------------------------------------------------------------------------------------------------------------
# IGVexit <- function(sock)
# {
#    .socketWrite(sock,paste('exit', '\n'))
# }
# #------------------------------------------------------------------------------------------------------------------------
# IGVping <- function(sock)
# {
#    .socketWrite(sock,paste('echo', '\n'))
# }
# #------------------------------------------------------------------------------------------------------------------------
# testRaw <- function()
# {
#    s <<- IGVsocket()
#    f <<- system.file(package="igvR", "extdata", "sample.bed")
#    IGVload(s, f)
#
#    for (i in 1:5){
#      IGVgoto(s, "chr5:80000-620000")
#      IGVgoto(s, "chr22:800-6200")
#      }
#
#
# } # testRaw
# #------------------------------------------------------------------------------------------------------------------------
#
