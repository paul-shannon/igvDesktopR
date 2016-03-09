\name{displayGWASTable}
\alias{displayGWASTable}
\title{displayGWASTable}

\description{
  Create a Manhattan plot from a data.frame, with at least these
  columns, each strictly defined and named:

  \itemize{
    \item   CHR: chromosome (aliases chr, chromosome)
    \item    BP: nucleotide location (aliases bp, pos, position)
    \item   SNP: SNP identifier (aliases snp, rs, rsid, rsnum, id, marker, markername)
    \item     P: p-value for the association (aliases p, pval, p-value, pvalue, p.value)
    }
}

\arguments{
  \item{igv}{an \code{igvR} object. }
  \item{table}{a \code{data.frame}.}
  \item{name}{a \code{character sring}.}
}

\value{
  none.
}
\author{Paul Shannon}

\examples{
   \dontrun{
      chr <- rep("chr1", count)
      bp <- seq(from=1000, by=10, length.out=10)
      snp <- paste("rs", 1:10, sep="")
      exponents <- as.integer(runif(count, 3, 40))
      base <- runif(count, 1, 10)
      p <- unlist(lapply(1:count, function(i) base[i]^-exponents[i]))
      tbl <- data.frame(CHR=chr, BP=bp, SNP=snp, P=p, stringsAsFactors=FALSE)
      displayGWASTable(igv, tbl, name="10 gwas snps")
      goto(igv, "chr1", bp[1] - 20, bp[count] + 20)
      }
}


\seealso{
}


\keyword{Graphics}
