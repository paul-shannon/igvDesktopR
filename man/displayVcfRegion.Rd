\name{displayVcfRegion}
\alias{displayVcfRegion}
\title{displayVcfRegion}

\description{
  Display a (small!) vcf subset (by region, by sampleIDs).

}

\arguments{
  \item{igv}{an \code{igvR} object. }
  \item{chrom}{a \code{data.frame}.}
  \item{start}{a \code{numeric}.}
  \item{end}{a \code{numeric}.}
  \item{vcfDirecotry}{a \code{character string}.}
  \item{sampleIDs}{a \code{list of character srings}.}
 }

\value{
  none.
}
\author{Paul Shannon}

\examples{
   \dontrun{
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
      }
}


\seealso{
}


\keyword{Graphics}
