#----------------------------------------------------------------------------------------------------
.igvR <- setClass ("igvR",
                   representation = representation(socket="socket",
                                                   genome="character",
                                                   quiet="logical")
                    )

#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
setGeneric('connected',             signature='obj', function(obj) standardGeneric ('connected'))
setGeneric('disconnect',            signature='obj', function(obj) standardGeneric ('disconnect'))
setGeneric('getGenome',             signature='obj', function(obj) standardGeneric('getGenome'))
setGeneric('getVcfDirectory',       signature='obj', function(obj) standardGeneric('getVcfDirectory'))
setGeneric('loadFile',              signature='obj', function(obj, filename) standardGeneric ('loadFile'))
setGeneric('goto',                  signature='obj', function(obj, chrom, start, end) standardGeneric ('goto'))
setGeneric('gotoAll',               signature='obj', function(obj) standardGeneric ('gotoAll'))
setGeneric('clearAllTracks',        signature='obj', function(obj) standardGeneric ('clearAllTracks'))
setGeneric('displayBedTable',       signature='obj', function(obj, tbl, name) standardGeneric ('displayBedTable'))
setGeneric('displayGWASTable',      signature='obj', function(obj, tbl, name) standardGeneric ('displayGWASTable'))
setGeneric('displayScoredFeatures', signature='obj', function(obj, tbl, quiet=TRUE) standardGeneric ('displayScoredFeatures'))
setGeneric('displayVcfRegion',      signature='obj', function(obj, chrom, start, end, vcfFilename,
                                                              sampleIDs=character()) standardGeneric ('displayVcfRegion'))
#----------------------------------------------------------------------------------------------------
igvR <- function(host="localhost", port=60151, genome="hg38", quiet=TRUE)
{
  socket = try(make.socket(host,port))

  if(inherits(socket,'try-error')) {
    stop(sprintf("Make sure that IGV is running on host %s and that IGV is set to accept commands on port %d",host,port))
    }

  .igvR(socket=socket, genome=genome, quiet=quiet)

} # igvR, the constructor
#----------------------------------------------------------------------------------------------------
setMethod('connected', 'igvR',

  function (obj){
     msg <- paste("echo", "\n")
     result <- .send(obj@socket, msg)
     result == "echo\n";
     })

#----------------------------------------------------------------------------------------------------
setMethod('disconnect', 'igvR',

  function (obj){
     close(obj@socket)
     })

#----------------------------------------------------------------------------------------------------
setMethod('getGenome', 'igvR',

  function (obj){
     obj@genome
     })

#----------------------------------------------------------------------------------------------------
setMethod('loadFile', 'igvR',

  function (obj, filename){
     msg <- paste("load",filename, "\n")
     if(!obj@quiet) printf("igvR loadFile, filename %s, msg %s", filename, msg)
     result <- .send(obj@socket, msg)
     invisible(result)
     })

#----------------------------------------------------------------------------------------------------
setMethod('goto', 'igvR',

  function (obj, chrom, start, end){
     region <- sprintf("%s:%d-%d", chrom, start, end);
     msg <- paste("goto", region, "\n")
     result <- .send(obj@socket, msg)
     invisible(result)
     })

#----------------------------------------------------------------------------------------------------
setMethod('gotoAll', 'igvR',

  function (obj){
     msg <- paste("goto all","\n")
     result <- .send(obj@socket, msg)
     invisible(result)
     })

#----------------------------------------------------------------------------------------------------
setMethod('clearAllTracks', 'igvR',

  function (obj){
     msg <- paste("new","\n")
     result <- .send(obj@socket, msg)
     invisible(result)
     })

#----------------------------------------------------------------------------------------------------
setMethod('displayBedTable', 'igvR',

  function (obj, tbl, name){
        # some simple sanity checks.
     stopifnot("data.frame" %in% is(tbl))
     stopifnot(ncol(tbl) >= 3)
     stopifnot(length(grep("^chr", tbl[,1])) == nrow(tbl))
     stopifnot("numeric" %in% is(tbl[,2]))
     stopifnot("numeric" %in% is(tbl[,3]))

     name <- gsub(" ", ".", name);
     tempFile <- sprintf("%s/%s.bed", tempdir(), name)
     write.table(tbl, file=tempFile, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
     result <- loadFile(obj, tempFile)
     invisible(result)
     })

#----------------------------------------------------------------------------------------------------
# use the igv format: https://www.broadinstitute.org/software/igv/IGV
# An IGV file (.igv) is a tab-delimited text file that defines
# tracks. The first row contains column headings for chromosome, start
# location, end location, and feature followed by the name of each
# track defined in the .igv file. Each subsequent row contains a locus
# and the associated numeric values for each track. IGV interprets the
# first four columns as chromosome, start location, end location, and
# feature name regardless of the column headings in the file. IGV uses
# the column headings for the fifth and subsequent columns as track
# names. Feature names are not displayed in IGV.

setMethod('displayScoredFeatures', 'igvR',

  function (obj, tbl, quiet=TRUE){
        # some simple sanity checks.
     stopifnot("data.frame" %in% is(tbl))
     stopifnot(ncol(tbl) >= 5)
     stopifnot(length(grep("^chr", tbl[,1])) == nrow(tbl))
     stopifnot(tolower(colnames(tbl)[1]) %in% c("chr", "chrom", "chromosome", "seqname", "seqnames"))
     stopifnot(tolower(colnames(tbl)[2:3]) == c("start", "end"))
     stopifnot("numeric" %in% is(tbl[,2]))
     stopifnot("numeric" %in% is(tbl[,3]))
     stopifnot("numeric" %in% is(tbl[,5]))

     tbl <- tbl[order(tbl$start),]
     tempFile <- tempfile(fileext=".igv")

     if(!obj@quiet){
        message(sprintf("displayScoredFeatures writing to temp file %s", tempFile))
        }

     write.table(tbl, file=tempFile, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
     result <- loadFile(obj, tempFile)
     invisible(result)
     })

#----------------------------------------------------------------------------------------------------
# see http://www.broadinstitute.org/software/igv/GWAS for table format information

setMethod('displayGWASTable', 'igvR',

  function (obj, tbl, name){
        # some simple sanity checks.
     stopifnot("data.frame" %in% is(tbl))
     stopifnot(all(c("CHR", "BP", "SNP", "P") %in% colnames(tbl)))
     stopifnot(length(grep("^chr", tbl[, "CHR"])) == nrow(tbl))
     stopifnot("numeric" %in% is(tbl[,"BP"]))
     stopifnot("numeric" %in% is(tbl[, "P"]))
     stopifnot(all(tbl$P > 0))

     name <- gsub(" ", ".", name);
     tempFile <- sprintf("%s/%s.gwas", tempdir(), name)
     write.table(tbl, file=tempFile, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
     result <- loadFile(obj, tempFile)
     invisible(result)
     })

#----------------------------------------------------------------------------------------------------
setMethod('displayVcfRegion', 'igvR',

  function(obj, chrom, start, end, vcfFilename, sampleIDs=character()) {
      stopifnot(file.exists(vcfFilename))
      if(!obj@quiet) printf("  chrom filenames: '%s'", paste(filenames, collapse=' ,'))
      tbiFilename <- paste(vcfFilename, ".tbi", sep="")
      stopifnot(file.exists(tbiFilename))
      gr <- GRanges(chrom, IRanges(start, end))
      print(ranges(gr))
      params <- ScanVcfParam(which=gr, samples=sampleIDs)
      vcf <- readVcf(TabixFile(vcfFilename), obj@genome, params)
      #browser()
      tempFile <- tempfile(fileext=".vcf")
      if(!obj@quiet) printf("  about to write region of vcf to '%s'", tempFile)
      writeVcf(vcf, tempFile)
      if(!obj@quiet) printf("   about to gzip");
      tempFile.gz <- sprintf("%s.gz", tempFile)
      bgzip(tempFile, dest=tempFile.gz, overwrite=TRUE)
      if(!obj@quiet) printf("   about to index");
      tempFile.gz.tbi <- indexTabix(tempFile.gz, format="vcf")
      loadFile(obj, tempFile.gz)
      }) # vcfRegionDisplay

#------------------------------------------------------------------------------------------------------------------------
.send <- function(socket, text)
{
  write.socket(socket, text)
  response <- read.socket(socket)

  response

} # .send
#----------------------------------------------------------------------------------------------------
#
# IGVsocket <- function(host='localhost', port=60151) {
#   sock = try(make.socket(host,port))
#   if(inherits(sock,'try-error')) {
#     stop(sprintf("Make sure that IGV is running on host %s and that IGV is set to accept commands on port %d",host,port))
#   }
#   print(sock)
#   return(sock)
# }
#
# .socketWrite<-
#   function(sock,string) {
#   print(string)
#   write.socket(sock,string)
#   response <- read.socket(sock)
#   return(response)
# }
#
#
# IGVload <-
# function (sock, files) {
#   if(length(files)<1) stop('Files must be specified')
#   message(basename(files))
#   .socketWrite(sock,paste('load',paste(files,collapse=','),'\n'))
# }
#
# IGVgoto <-
#   function(sock,region) {
#     .socketWrite(sock, paste('goto',region,'\n'))
#   }
#
# IGVgenome <-
#   function(sock,genome="hg18") {
#     .socketWrite(sock,paste('genome',genome,'\n'))
#   }
#
# IGVsnapshot <-
#   function(sock,fname="",dirname=getwd()) {
# 	.socketWrite(sock,paste('snapshotDirectory',dirname,'\n'))
#     .socketWrite(sock, paste('snapshot',fname,'\n'))
#   }
#
# IGVclear <-
#   function(sock) {
#     .socketWrite(sock,paste('new \n'))
#   }
#
# ## option can be base, position, strand, quality, sample, and readGroup
# IGVsort <-
#   function(sock, option) {
#     .socketWrite(sock,paste('sort', option,'\n'))
#   }
#
#  IGVcollapse <-
#   function(sock) {
#     .socketWrite(sock,paste('collapse \n'))
#   }
#
