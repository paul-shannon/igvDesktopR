#----------------------------------------------------------------------------------------------------
.igvR <- setClass ("igvR",
                   representation = representation(socket="socket",
                                                   genome="character")
                    )

#----------------------------------------------------------------------------------------------------
setGeneric('loadFile',         signature='obj', function(obj, filename) standardGeneric ('loadFile'))
setGeneric('goto',             signature='obj', function(obj, chrom, start, end) standardGeneric ('goto'))
setGeneric('displayBedTable',  signature='obj', function(obj, tbl, name) standardGeneric ('displayBedTable'))
setGeneric('displayGWASTable',  signature='obj', function(obj, tbl, name) standardGeneric ('displayGWASTable'))
setGeneric('connected',        signature='obj', function(obj) standardGeneric ('connected'))
#----------------------------------------------------------------------------------------------------
igvR <- function(host="localhost", port=60151, genome="hg38")
{
  socket = try(make.socket(host,port))

  if(inherits(socket,'try-error')) {
    stop(sprintf("Make sure that IGV is running on host %s and that IGV is set to accept commands on port %d",host,port))
    }

  .igvR(socket=socket, genome=genome)

} # igvR, the constructor
#----------------------------------------------------------------------------------------------------
setMethod('connected', 'igvR',

  function (obj){
     msg <- paste("echo", "\n")
     result <- .send(obj@socket, msg)
     result == "echo\n";
     })

#----------------------------------------------------------------------------------------------------
setMethod('loadFile', 'igvR',

  function (obj, filename){
     msg <- paste("load",filename, "\n")
     printf("igvR loadFile, filename %s, msg %s", filename, msg)
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
setMethod('displayBedTable', 'igvR',

  function (obj, tbl, name){
        # some simple sanity checks.
     stopifnot("data.frame" %in% is(tbl))
     stopifnot(ncol(tbl) >= 3)
     stopifnot(length(grep("^chr", tbl[,1])) == nrow(tbl))
     browser()
     stopifnot("numeric" %in% is(tbl[,2]))
     stopifnot("numeric" %in% is(tbl[,3]))

     name <- gsub(" ", ".", name);
     tempFile <- sprintf("%s/%s.bed", tempdir(), name)
     write.table(tbl, file=tempFile, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
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

     name <- gsub(" ", ".", name);
     tempFile <- sprintf("%s/%s.gwas", tempdir(), name)
     write.table(tbl, file=tempFile, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
     result <- loadFile(obj, tempFile)
     invisible(result)
     })

#----------------------------------------------------------------------------------------------------
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
