# Plot variables within a VCF across its length.
library(data.table)
library(reshape2)
library(tidyr)
library(dplyr)
library(scales)
library(ggplot2)

setwd("~/Documents/git/vcf-toolbox/")

vcf <- "00_ALL_m.filter.lcr.bcf"


var <- c("DP", "QUAL")

#==================================================================#
# Get Sample Data - Retrieves variables for all samples in a vcf   #
#==================================================================#
get_sample_data <- function(vcf, var, fast_sample=F) {
  # Retrieve data for a sample variable.
  var_opts <- list(depth="%DP", dp="%DP")
  
  # Take sample rather than all data.
  if (fast_sample == T) {
    fast_sample <- "| awk 'NR % 200 == 0'"
  } else {
    fast_sample <- ""
  }
  
  # Set the variable for pre-defined terms, otherwise use user input.
  for (x in range(1,length(var))) {
    print(x)
    if (is.null(var_opts[[tolower(var[x])]]) == F) {
      var[x] <- var_opts[[tolower(var[x])]]
    }
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


#==================================================================#
# Get Variant Data - Retrieves variables for all variants in a vcf #
#==================================================================#
get_variant_data <- function( vcf, var, fast_sample = F, no_multi = T) {
  var_names <- c("CHROM","POS",make.names(gsub("%","",var)))
  var# Retrieve data for a variant variable.
  var_opts <- list(depth="%DP", "quality"="%QUAL", qual="%QUAL", "%qual" = "%QUAL", "%dp"="%DP", "dp"= "%DP")
  
  # Set the variable for pre-defined terms, otherwise use user input.
  for (x in seq(1,length(var))) {
    if (is.null(var_opts[[tolower(var[x])]]) == F) {
      var[x] <- var_opts[[tolower(var[x])]]
    }
  }
  
  # Sample 1/200 / For quick sample.
  if (fast_sample == T) {
    fast_sample <- "| awk 'NR % 200 == 0'"
  } else {
    fast_sample <- ""
  }
  
  # Filter multiallelic variants
  if (no_multi == T) {
    no_multi = sprintf("bcftools view -O b -M 2 %s |", vcf)
    vcf <- ""
  } else {
    no_multi <- ""
  }
  
  tmp <- tempfile()
  system(sprintf("%s bcftools query -f '%%CHROM\t%%POS%s\n' %s %s > %s",no_multi, paste0("\t", var, collapse=""), vcf, fast_sample, tmp))
  r <- fread(tmp)
  setnames(r, var_names)  
  # Rename Variables  
  r
}


list_variant_variables <- function( vcf ) {
  system(sprintf("bcftools view %s | head -n 1000 | grep '##INFO'", vcf), wait=F)
}

genetic_scale <- function(n) {
  paste0(n/100000, "Mb")
}

plot_variant_data <- function( vcf, var, fast_sample = F, scale = "equal") {
    representation <- list(equal = "fixed", actual = "free_x")
    r <- get_variant_data(vcf = vcf, var = var, fast_sample = fast_sample)
    ggplot(r) +
      geom_point( aes(x=V2,y=V3, alpha=0.2)) +
      facet_grid(~V1, scales="free_x", space=representation[[scale]], margins=F) +
      theme_bw() +
      labs(x="Chromosome", y = var) +
      theme(legend.position="none", panel.grid.major = element_line(colour = "#4c4c4c", linetype = "dotted")) +
      stat_smooth( aes(x=V2,y=V3), fill = "grey50", size=2, color="#0080ff") +
      scale_x_continuous(labels = genetic_scale)
}



r <- get_variant_data(vcf, var = c("%DP", "quality", "%AC"))
r2 <- get_variant_data("00_ALL_c.bcf", var = "%DP")

plot_variant_data(vcf, var="%DP")


#
# MT Analysis (Temp)
#

stat_sum_single <- function(fun, geom="point", ...) {
  stat_summary(fun.y=fun, colour="red", geom=geom, size = 2, ...)
}

# Depth by Autosomes vs. MT
ggplot(r, aes(x=AC, y = quality, alpha=0.01)) +
  geom_point() +
  stat_sum_single(mean) +
  scale_y_continuous(labels=comma) +
  labs(title="ALT allele count vs Quality")

ggplot(r, aes(x=AC, y = DP, alpha=0.01)) +
  geom_point() +
  stat_sum_single(mean) +
  scale_y_continuous(labels=comma) +
  labs(title="ALT allele count vs DP")

ggsave("ALT_vs_DP.png")
                 

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
