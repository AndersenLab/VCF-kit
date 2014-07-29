# Plot variables within a VCF across its length.
library(data.table)
library(reshape2)
library(tidyr)
library(dplyr)

setwd("~/Documents/git/vcf-toolbox/")

vcf <- "00_ALL_c.bcf"
plotvar <- "%QUAL"

get_sample_data <- function(vcf, var, fast_sample=F) {
  # Retrieve data for a sample variable.
  var_opts <- list(depth="%DP")
  
  if (fast_sample == T) {
    fast_sample <- "| awk 'NR % 200 == 0'"
  } else {
    fast_sample <- ""
  }
  
  # Set the variable for pre-defined terms, otherwise use user input.
  if (is.na(var_opts[tolower(var)]) == F) {
    var <- as.vector(var_opts[tolower(var)])
  }
  
  # Function for returning depth for every sample
  tmp <- tempfile()
  # Get sample names
  samples <- read.table(pipe(sprintf("bcftools query -l %s", vcf)), col.names=c("sample"), stringsAsFactors=F )
  system(sprintf("bcftools query -f '%%CHROM\t%%POS[\t%s]\n' %s %s > %s", var, vcf, fast_sample, tmp))
  r <- fread(tmp,header=F)
  setnames(r,c("CHROM", "POS", samples$sample))
  r
}

get_variant_data <- function( vcf, var, fast = F) {
  # Retrieve data for a variant variable.
  var_opts <- list(depth="%DP", quality="%QUAL")
  
  # Set the variable for pre-defined terms, otherwise use user input.
  if (is.na(var_opts[tolower(var)]) == F) {
    var <- as.vector(var_opts[tolower(var)])
  }
  
  # Sample 1/200 / For quick sample.
  if (fast_sample == T) {
    fast_sample <- "| awk 'NR % 200 == 0'"
  } else {
    fast_sample <- ""
  }
  
  comm <- sprintf("bcftools query -f '%%CHROM\t%%POS\t%s\n' %s %s", plotvar, vcf, fast_sample)
  r <- read.csv2(pipe(comm), as.is=T, sep='\t',header=F, check.names=FALSE)  
  r$V3 <- as.numeric(r$V3)
}

p <- get_variant_data(vcf, "quality")
p <- get_sample_data(vcf, "%HWE", )

plot_variant_data <- function( vcf, var, fast = F, scale = "equal") {
    representation <- list(equal = "fixed", actual = "free_x")
    r <- get_variant_data()
}


ggplot(r) +
  geom_point( aes(x=V2,y=V3, alpha=0.2)) +
  facet_grid(~V1, scales="free_x", space=representation[[scale]], margins=F) +
  stat_smooth( aes(x=V2,y=V3),  method="loess", fill = "grey50", size=2, color="#0080ff")


plot_vcf(vcf, "%QUAL", fast=T)



sample_depth <- get_sample_depth(vcf)
sample_depth$mt[sample_depth$CHROM == "chrM"] <- 1




# Depth by Autosomes vs. MT

m <- sample_depth %>%
  group_by(mt) %>%
  select(-matches("CHROM"), -matches("-POS")) %>%
  mutate(count = length(mt)) %>%
  summarise_each(funs(mean))

mito <- m[2]/m[1]
mito <- mito %>%
        select(-contains("mt"), -contains("POS"), -contains("count"))

# Reshape



# Depth by Chromosome
s <- group_by(sample_depth, CHROM) %>%
     summarise_each(funs(mean))



# Plot ratios
s <- ungroup(s) %>%
     group_by(mt)
