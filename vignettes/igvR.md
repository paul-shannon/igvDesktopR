<!--
% \VignetteIndexEntry{The igvR Package}
% \VignetteDepends{igvR}
% \VignetteEngine{knitr::knitr}
-->

<script type="text/javascript">
document.addEventListener("DOMContentLoaded", function() {
  document.querySelector("h1").className = "title";
});
</script>
<script type="text/javascript">
document.addEventListener("DOMContentLoaded", function() {
  var links = document.links;  
  for (var i = 0, linksLength = links.length; i < linksLength; i++)
    if (links[i].hostname != window.location.hostname)
      links[i].target = '_blank';
});
</script>

# The igvR Package

**Package**: *[igvR](http://bioconductor.org/packages/release/bioc/html/igvR.html)*<br />
**Authors**: Paul Shannon<br />
**Modified**: 1 March, 2016<br />
**Compiled**: Tue Mar  1 12:26:04 2016

## Overview

This R package provides a simple connection between an R session and an (already running) instance of desktop *IGV*,
the [Integrated Genomics Viewer](http://www.broadinstitute.org/igv/home) from the Broad Institute, created by
Jim Robinson and colleagues.  We use the [socket interface](https://www.broadinstitute.org/igv/PortCommands).

## Session info

Here is the output of `sessionInfo()` on the system on which this
document was compiled:


```
## R version 3.2.3 (2015-12-10)
## Platform: x86_64-apple-darwin13.4.0 (64-bit)
## Running under: OS X 10.11.2 (El Capitan)
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] BiocStyle_1.8.0      knitr_1.12.3         BiocInstaller_1.20.1
## 
## loaded via a namespace (and not attached):
## [1] magrittr_1.5  formatR_1.2.1 tools_3.2.3   stringi_1.0-1 stringr_1.0.0
## [6] evaluate_0.8
```
