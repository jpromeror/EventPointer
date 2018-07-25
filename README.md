# EventPointer
*EventPointer* is an R package to identify alternative splicing events 
		that involve either simple (case-control experiment) or complex experimental designs 
		such as time course experiments and studies including paired-samples. The algorithm can
		be used to analyze data from either junction arrays (Affymetrix Arrays) or sequencing data (RNA-Seq). 
		The software returns a data.frame with the detected alternative splicing 
		events: gene name, type of event (cassette, alternative 3',...,etc), genomic 
		position, statistical significance and increment of the percent spliced in (Delta PSI) for all 
		the events.
		The algorithm can generate a series of files to visualize the detected alternative 
		splicing events in IGV. This eases the interpretation of results and the design 
		of primers for standard PCR validation.

![**Figure 1.** EventPointer pipeline ](https://github.com/jpromeror/EventPointer/blob/master/vignettes/Figure1.png)

# Installation
EventPointer can be installed from Bioconductor using the BiocManager package:

```{r, eval=FALSE}
library(BiocManager)
BiocManager::install("EventPointer")
```