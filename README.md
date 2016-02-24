# igvR
a simple R interface to igv (http://igv.org) a java program for viewing genome tracks

sample use, with IGV already running on localhost
````
library(igvR)
igv <- igvR()   # localhost, standard port, hg38 are default values
connected(igv)  # TRUE
````
Load a tiny 2-line bed file included with the package, from the UCSC website BED format documentation  
````  
bed.file <- system.file(package="igvR", "extdata", "sample.bed")
stopifnot(file.exists(bed.file))
loadFile(igv, bed.file)
goto(igv, "chr22", 500, 6500)
````

Create and display a 3-row, 6-column data.frame (we will add Bioconductor GenomicRanges support soon)
````    
chr    <- rep("chr1", 3)
start  <- c(100470168, 100471168, 100472168)
end    <- c(100470193, 100471175, 100472306)
name   <- c("MA0679.1", "MA0756.1", "MA0757.1")
score  <- c(56.0, 58.6, 60.5)
strand <- rep("+", 3)
tbl <- data.frame(chr=chr, start=start, end=end, name=name, score=score, strand=strand, stringsAsFactors=FALSE)
displayBedTable(igv, tbl, name="3 motifs")
goto(igv, "chr1", start[1] - 200, end[3] + 200)
````    

Create and display a GWAS table, with 10 SNPs

````
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

   result <- displayGWASTable(igv, tbl, name="10 gwas snps")
   checkEquals(result, "OK\n")
   result <- goto(igv, "chr1", bp[1] - 20, bp[count] + 20)
