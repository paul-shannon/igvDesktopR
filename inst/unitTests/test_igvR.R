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
   test_sampleBedFile()
   test_goto()
   test_displayBedTable()

} # runTests
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

  full.path <- system.file(package="igvR", "extdata", "sample.bed")
  result <- loadFile(igv, full.path)
  checkEquals(result, "OK\n");

} # test.tinyBedFile
#------------------------------------------------------------------------------------------------------------------------
test_goto <- function()
{
   printf("--- test_goto")
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

   result <- displayBedTable(igv, tbl, name="3 motifs")
   checkEquals(result, "OK\n")
   result <- goto(igv, "chr1", start[1] - 20, end[1] + 20)
   checkEquals(result, "OK\n")

} # test_displayBedTable
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
