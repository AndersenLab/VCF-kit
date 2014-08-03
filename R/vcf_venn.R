library(stringr)
#=======#
# Input #
#=======#

args<-commandArgs(TRUE)
f1 <- args[2]
f2 <- args[3]
clean_names <- gsub("(.bcf|.vcf|.vcf.gz|.gz)","", c(f1, f2))
query <- args[4]

setwd("./")

f1 <- "../00_ALL_m.filter.lcr.bcf" 
f2 <- "../00_ALL_c.filter.lcr.filtered.bcf"

#===========================#
# Data Generation Functions #
#===========================#

import_table <- function(t_name, f) { 
  # This function can import data from bcftools stats function
  t <- read.csv2(pipe(sprintf("egrep '^(# )%s[^,]|^%s' %s", t_name, t_name, f)), as.is=T, sep='\t',header=T, check.names=FALSE)
  names(t) <- lapply(names(t), function(x) make.names(gsub("\\[[1-9]+\\]","",x)))
  # Remove first column and return table.
  t[,c(-1)]
}

get_pair_stats <- function(f1, f2) {
  # This function retrieves comparison data for
  # each VCF.
  tmp <- tempfile()
  system(sprintf('bcftools stats -s - %s %s > %s', f1, f2, tmp))
  SN  <- import_table("SN", tmp)
  
   # Set up names.
  ids <- gsub("(.bcf|.vcf|.vcf.gz|.gz)","",list("0"=f1,"1"=f2,"2"=sprintf("%s__%s", f1,f2)))
  
  # Fix SN Group, ID
  SN <- reshape(SN, direction="wide", timevar=c("key"), ids=c("id"))
  names(SN) <- make.names(gsub("(number.of.|value.|:)","",names(SN)))

  SN$SNP_Total[SN$id == 0] <- SN[SN$id == 0,"SNPs"] + SN[SN$id == 2,"SNPs"]  
  SN$SNP_Total[SN$id == 1] <- SN[SN$id == 1,"SNPs"] + SN[SN$id == 2,"SNPs"] 
  SN$SNP_Total[SN$id == 2] <- SN[SN$id == 2,"SNPs"] # Intersection!
  
  
  SN  
}

#==========#
# Get Data #
#==========#

f <- get_pair_stats(f1, f2)

#==================#
# Plot Data        #
#==================#
png(sprintf("%s-%s.png", clean_names))
grid.newpage()
  venn.plot <- draw.pairwise.venn(area1        = f[f$id==0,"SNP_Total"],
                                  area2        = f[f$id==1,"SNP_Total"],
                                  cross.area   = f[f$id==2,"SNP_Total"],
                                  scaled       = T,
                                  category     = clean_names,
                                  fill         = c("blue", "red"),
                                  alpha        = 0.3,
                                  lty          = "blank",
                                  cex          = 2,
                                  cat.cex      = 2,
                                  cat.pos      = c(0, 0),
                                  cat.dist     = 0.05,
                                  ext.pos      = 30,
                                  ext.dist     = -0.05,
                                  ext.length   = 0.85,
                                  ext.line.lwd = 2,
                                  ext.line.lty = "dashed")
grid.draw(venn.plot)
dev.off()