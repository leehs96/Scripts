# DESeq2

<em>1. requirements :
- counts matrix
- condition file 
 
 
 
 
 
 
 
 
 ![image](https://user-images.githubusercontent.com/63722122/120759790-576fd580-c54e-11eb-8a67-f9736797cd54.png)


 
 
<em>2. optional(if batch correction necessary) :
- batch information must be in condition file
 

 # OUTPUT
 - batch corrected count matrix
 - DESeq2 results (GeneID, ensID, Log2Foldchange,padj)
 - distance plot
 - PCA plot
 - Volcanoplot
 - Top10 Foldchange plot
 - enriched KEGG pathway
 - enriched GO
 - enriched cell marker
 - enriched KEGG, GO, cellmarker (FC > 1)
 - enriched KEGG, GO, cellmarker (FC > -1)
 
 
 
 
 
# MUST SET OWN DIRECTORY PATH! 

![image](https://user-images.githubusercontent.com/63722122/120759375-ec260380-c54d-11eb-8543-6d58f1c5beb7.png)
 
![image](https://user-images.githubusercontent.com/63722122/120759089-98b3b580-c54d-11eb-81ad-ff8756d0a010.png)




[DESeq2 Document] (http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
 
[ComBat-seq github] (https://github.com/zhangyuqing/ComBat-seq)
 
[EnhancedVolcano document] (https://bioconductor.org/packages/release/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html)
 
[pheatmap] (https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/)
 
[ggplot2] (https://statkclee.github.io/R-ecology-lesson/kr/05-visualization-ggplot2.html)
 
 

