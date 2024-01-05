library(Ravages)
# Import data in a bed matrix
x <- as.bed.matrix(x=LCT.matrix.bed, fam=LCT.matrix.fam, bim=LCT.snps)
# Add population
x@ped[,c("pop", "superpop")] <- LCT.matrix.pop1000G[,c("population", "super.population")]
x.CADDregions <- set.CADDregions(x)
# Annotation of variants with adjusted CADD scores
x <- adjustedCADD.annotation(x)
# Keep only CADD regions with 2 variants and variants with a MAF greater than 1%
# and with an adjusted CADD score greater than the median
x.median <- filter.adjustedCADD(x.CADDregions, maf.threshold = 0.01, min.nb.snps = 2)
x.burden <- burden.subscores(x.median, x1.H0.burden, burden.function = WSS)
# SKAT analysis
x.SKAT <- SKAT(x.median, x1.H0.SKAT)
