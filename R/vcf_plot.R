# Plot variables within a VCF across its length.
library(data.table)
library(reshape2)
library(tidyr)
library(dplyr)

setwd("~/Documents/git/vcf-toolbox/")

vcf <- "00_ALL_c.bcf"
plotvar <- "%QUAL"

get_sample_depth <- function(vcf) {
  # Function for returning depth for every sample
  tmp <- tempfile()
  # Get sample names
  samples <- read.table(pipe(sprintf("bcftools query -l %s", vcf)), col.names=c("sample"), stringsAsFactors=F )
  system(sprintf("bcftools query -f '%%CHROM\t%%POS[\t%%DP]\n' %s > %s", vcf, tmp))
  r <- fread(tmp,header=F)
  setnames(r,c("CHROM", "POS", samples$sample))
}


plot_vcf <- function( vcf, plotvar, fast = F, scale = "equal") {
  plotvar_opts <- list(depth="%DP", "%dp"="%DP", quality="%QUAL", "%qual"="%QUAL")
  representation <- list(equal = "fixed", actual = "free_x")
  
  plotvar <- as.vector(plotvar_opts[tolower(plotvar)])
  
  # Sample 1/200 of variants for faster plotting.
  fast_sample <- ""
  if (fast == T) {
    fast_sample <- "| awk 'NR % 200 == 0'"
  }
  
  comm <- sprintf("bcftools query -f '%%CHROM\t%%POS\t%s\n' %s %s", plotvar, vcf, fast_sample)
  r <- read.csv2(pipe(comm), as.is=T, sep='\t',header=F, check.names=FALSE)
  
  r$V3 <- as.numeric(r$V3)
  
  ggplot(r) +
    geom_point( aes(x=V2,y=V3, alpha=0.2)) +
    facet_grid(~V1, scales="free_x", space=representation[[scale]], margins=F) +
    stat_smooth( aes(x=V2,y=V3),  method="loess", fill = "grey50", size=2, color="#0080ff")
}

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
