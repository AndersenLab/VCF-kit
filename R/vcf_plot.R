library(ggplot2)
library(scales)
library(stringr)
library(dplyr)
library(grid)
library(colourlovers)


# Functions for working with bcftools + VCF
#=======================#
# Import BCFTools Table #
#=======================#
import_table <- function(t_name, f) { 
  # This function can import data from bcftools stats function
  t <- read.csv2(pipe(sprintf("egrep '^(# )%s[^,]|^%s' %s", t_name, t_name, f)), as.is=T, sep='\t',header=T, check.names=FALSE)
  names(t) <- lapply(names(t), function(x) make.names(gsub("\\[[1-9]+\\]","",x)))
  # Remove first column and return table.
  t[,c(-1)]
}

#==========================#
# Ancillary Functions/Defs #
#==========================#
number_ticks <- function(n) {function(limits) pretty(limits, n)}

repl_multiple <- function(string, replacements) {
  # Convenience Function for 
  # replacing multiple string with nothing.
  for (rep in replacements) {
    string <- sub(rep, "", string)
  }
  string
}

#=================#
# Fetch Data      #
#=================#

# Get individual VCF Stats
get_vcf_stats <- function(f) {
  # Retrieve VCF stats for a single vcf.
  tmp <- tempfile()
  f_name <- gsub("(.bcf|.vcf.gz|.vcf|.gz|.txt)","",f)
  system(sprintf('bcftools stats -s - %s > %s', f, tmp))
  stats <- list()
  tables <- system(sprintf("cat %s | grep -v '#' | cut -f 1 | uniq", tmp), intern=T)
  for (t in tables[-1]) {
    stats[[t]] <- import_table(t, tmp)
    stats[[t]]$id <- f_name
  }
  # Fix SN Group, ID
  stats[["SN"]] <- reshape(stats[["SN"]], direction="wide", timevar=c("key"), ids=c("id"))
  names(stats[["SN"]]) <- make.names(gsub("value.|number.of.","",gsub(":","", names(stats[["SN"]]))))
  stats
}


#=================#
# Singletons      #
#=================#

plot_singletons <- function(s, palette="919313", title="Singletons", xlab="Strain", ylab="Singletons") {
  # Plots singleton data.

  palette <- paste0("#",rbind(clpalette(palette)$colors))
  palette <- rep(palette, (round(nrow(s$PSC)/length(palette),0)+1))

  ggplot(s$PSC) + 
    geom_point(aes(x=sample, y=nSingletons, color=sample), colour="black", size = 6) +
    geom_point(aes(x=sample, y=nSingletons, color=sample), size=5) +
    labs(title=title, x=xlab, y=ylab) +
    theme(legend.position="top", axis.text.x = element_text(angle = 90, hjust = 1))  +
    theme_bw() +
    scale_color_manual(values=palette)
}

plot_singletons(s)