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
bed.file <- system.path(package="igvR", "extdata", "sample.bed")  # 2 lines only, from UCSC website
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
