# Plot variables within a VCF across its length.
library(data.table)
library(reshape2)
library(tidyr)
library(dplyr)
library(scales)
library(ggplot2)

setwd("~/Documents/git/vcf-toolbox/")

vcf <- "IND_BGI1-RET2-AB4-5cbd4-90f54.bam.het.bcf"
var <- c("DP", "QUAL")

#==================================================================#
# Get Variant Data - Retrieves variables for all variants in a vcf #
#==================================================================#
get_variant_data <- function( vcf, info_var, sample_var=NA, fast_sample = F, no_multi = T) {
  
  # Set up Variable Names
  q_vars <- c(info_var, sample_var)
  q_vars <- q_vars[!is.na(q_vars)]
  var_names <- c("CHROM","POS",make.names(gsub("%","",q_vars)))
  
  var_opts <- list(depth="%DP", "quality"="%QUAL", qual="%QUAL", "%qual" = "%QUAL", "%dp"="%DP", "dp"= "%DP")
  
  # Process Variant (INFO) Variables
  
  # Set the variable for pre-defined terms, otherwise use user input.
  for (x in seq(1,length(info_var))) {
    if (is.null(var_opts[[tolower(info_var[x])]]) == F) {
      info_var[x] <- var_opts[[tolower(info_var[x])]]
    }
  }
  
  
  # Process Sample Variables
  var_opts <- list(depth="%DP", dp="%DP")
  for (x in range(1,length(var))) {
    if (is.null(var_opts[[tolower(sample_var[x])]]) == F) {
      sample_var[x] <- var_opts[[tolower(sample_var[x])]]
    }
  }
  
  if (!is.na(sample_var)) {
    sample_var <- sprintf("[\t%s]", sample_var)
  } else {
    sample_var <- ""
  }
  
  # Get Contigs
  # Fast sample - only take beginning of chromosomes.
  if (fast_sample == T) {
    # Get contigs
    tmp <- tempfile()
    system(sprintf("bcftools view %s | head -n 500 | grep '##contig' > %s", vcf, tmp), wait=T)
    contigs <- gsub("(##contig=<ID=|,.*)","",readLines(tmp))
    fast_sample <- paste("--regions ",sprintf("%s", paste0(contigs,":1-1000000", collapse=",")), collapse=" ")
  } else {
    fast_sample <- ""
  }
  
  # Filter multiallelic variants
  if (no_multi == T) {
    # If fast sampling - put here, otherwise, put later.
    no_multi = sprintf("bcftools view -O b -M 2 %s %s |", vcf, fast_sample)
    fast_sample <- ""
    vcf <- ""
  } else {
    no_multi <- ""
  }
  
  tmp <- tempfile()
  system(sprintf("%s bcftools query %s -f '%%CHROM\t%%POS%s%s\n' %s > %s",no_multi, fast_sample, paste0("\t", info_var, collapse=""), sample_var, vcf, tmp))
  print(sprintf("%s bcftools query %s -f '%%CHROM\t%%POS%s%s\n' %s > %s",no_multi, fast_sample, paste0("\t", info_var, collapse=""), sample_var, vcf, tmp))
  r <- fread(tmp)
  setnames(r, var_names)  
  # Rename Variables  
  r
}


r <- get_variant_data(vcf, info_var = c("%DP", "quality"), sample_var = "%TGT", fast_sample=T)



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
