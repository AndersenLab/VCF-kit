library(data.table)
library(grid)
library(VennDiagram)

setwd("~/Documents/git/vcf-compare/")

load_vcf <- function(f) {
  tmp <- tempfile()
  cmd <- sprintf("bcftools query --include '%%TYPE=\"snp\" | %%TYPE=\"indel\"' -c both -f '%%CHROM\t%%POS\t%%TYPE\t%%REF\t%%ALT[\t%%SAMPLE=%%TGT]\n' %s | gawk -f vcf_comp.awk > %s", f, tmp)
  print(cmd)
  system(cmd)
  fread(tmp, sep="\t",verbose=T)
}

# Load vcfs, set keys, and merge.
vcf1 <- load_vcf("mmp.vcf.gz")
vcf2 <- load_vcf("12_new_caller.bcf")
setkeyv(vcf1, cols=c("CHR","POS","TYPE"))
setkeyv(vcf2, cols=c("CHR","POS","TYPE"))
jvcf <- merge(x=vcf1, y=vcf2, by=c("CHR","POS","TYPE"), all=TRUE)

venn <- function(t="SNP") {
grid.newpage()
venn.plot <- draw.pairwise.venn(area1        = sum(complete.cases(jvcf[TYPE == t ,c("REF.x"), with = F])),
                                area2        = sum(complete.cases(jvcf[TYPE == t,c("REF.y"), with = F])),
                                cross.area   = sum(complete.cases(jvcf[TYPE == t,c("REF.x","REF.y"), with = F])),
                                scaled       = T,
                                category     = c("First", "Second"),
                                fill         = c("blue", "red"),
                                alpha        = 0.3,
                                lty          = "blank",
                                cex          = 2,
                                cat.cex      = 2,
                                cat.pos      = c(0, 0),
                                cat.dist     = 0.20,
                                cat.just     = list(c(-1, -1), c(1, 1)),
                                ext.pos      = 30,
                                ext.dist     = -0.05,
                                ext.length   = 0.85,
                                ext.line.lwd = 2,
                                ext.line.lty = "dashed")
grid.draw(venn.plot)
}

venn()

# Pairwise Concordances
p <- t(as.data.frame(sapply(combn(intersect_samples,m=2, simplify=F), function(s) {
  u <- select(jvcf, starts_with(s[[1]]),starts_with(s[[2]])) 
  u <- u[!is.na(u[1]) & !is.na(u[4]),]
  match <-  length(filter(u, u[[1]] == u[[4]])[,c(1)])
  no_match <- length(filter(u, u[[1]] != u[[4]])[,c(1)])
  intersect_concordance <- (match)/(length(u[,1]))
  union_concordance <- (match)/union_variants_count
  list(s1=s[[1]], s2=s[[2]],match=match,no_match=no_match,intersect_concordance=intersect_concordance,union_concordance=union_concordance)
})))

# Individual Concordances
p <- t(as.data.frame(sapply(intersect_samples, function(s) {
  u <- select(jvcf, starts_with(s)) 
  u <- u[!is.na(u[1]) & !is.na(u[2]),]
  match <-  length(filter(u, u[[1]] == u[[2]])[,c(1)])
  no_match <- length(filter(u, u[[1]] != u[[2]])[,c(1)])
  intersect_concordance <- (match)/(length(u[,1]))
  union_concordance <- (match)/union_variants_count
  list(match=match,no_match=no_match,intersect_concordance=intersect_concordance,union_concordance=union_concordance)
})))


