'
an example of an association analysis and previous steps of data filtering using the dataset LCT.matrix
available with the package Ravages. This dataset containts data from the 1000Genome project in the locus
containing the Lactase gene. In this example, we look for an association between rare variants and the
european populations of 1000Genomes. The population of each individual is available in the dataframe
LCT.matrix.pop1000G. A classical analysis by gene is performed, and an analysis using the strategy “RAVAFIRST”.

'

# Load the Ravages package and convert data to a specific format
library(Ravages)
x <- as.bed.matrix(x = LCT.matrix.bed, fam = LCT.matrix.fam, bim = LCT.snps)

# Incorporate population information into the data
x@ped[, c("pop", "superpop")] <- LCT.matrix.pop1000G[, c("population", "super.population")]

# Filter data to retain only individuals from the European population
x <- select.inds(x, superpop == "EUR")
x@ped$pop <- droplevels(x@ped$pop)  # Drop unused levels to clean up the data

# Define genomic regions with a specific flank width
x <- set.genomic.region(x, flank.width = 500)

# Display the distribution of genomic regions before filtering rare variants
table(x@snps$genomic.region, useNA = "ifany")

# Filter out rare variants based on specified thresholds
x1 <- filter.rare.variants(x, filter = "whole", maf.threshold = 0.01, min.nb.snps = 10)

# Display the distribution of genomic regions after filtering rare variants
table(x1@snps$genomic.region, useNA = "ifany")

# Set up a null model for burden testing
x1.H0.burden <- NullObject.parameters(x1@ped$pop, ref.level = "CEU",
                                      RVAT = "burden", pheno.type = "categorical")

# Perform a burden test analysis on the filtered data
burden(x1, NullObject = x1.H0.burden, burden = "CAST", cores = 1)

# Set up a null model for SKAT testing
x1.H0.SKAT <- NullObject.parameters(x1@ped$pop, RVAT = "SKAT", pheno.type = "categorical")

# Perform a SKAT test analysis on the filtered data
SKAT(x1, x1.H0.SKAT, params.sampling = list(perm.target = 10, perm.max = 500))

