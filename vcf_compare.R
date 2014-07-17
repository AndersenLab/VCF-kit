library(data.table)
library(grid)
library(stringr)
library(VennDiagram)

setwd("~/Documents/git/vcf-compare/")

load_vcf <- function(f) {
  tmp <- tempfile()
  cmd <- sprintf("bcftools query --include '%%TYPE=\"snp\" | %%TYPE=\"indel\"' -c both -f '%%CHROM\t%%POS\t%%TYPE\t%%REF\t%%ALT[\t%%SAMPLE=%%TGT]\n' %s | gawk -f vcf_comp.awk > %s", f, tmp)
  print(cmd)
  system(cmd)
  fread(tmp, sep="\t",verbose=T)
}

f1 <- "00_all_bams.txt.vcf.gz"
f2 <- "mmp.vcf.gz"

# Load vcfs, set keys, and merge.
vcf1 <- load_vcf(f1)
vcf2 <- load_vcf(f2)
setkeyv(vcf1, cols=c("CHR","POS","TYPE"))
setkeyv(vcf2, cols=c("CHR","POS","TYPE"))
jvcf <- merge(x=vcf1, y=vcf2, by=c("CHR","POS","TYPE"), all=TRUE)

venn <- function(t="SNP") {
grid.newpage()
venn.plot <- draw.pairwise.venn(area1        = sum(complete.cases(jvcf[TYPE == t ,c("REF.x"), with = F])),
                                area2        = sum(complete.cases(jvcf[TYPE == t,c("REF.y"), with = F])),
                                cross.area   = sum(complete.cases(jvcf[TYPE == t,c("REF.x","REF.y"), with = F])),
                                scaled       = F,
                                category     = str_replace_all(c(f1, f2), c("(.vcf|.gz|.txt|.bcf)"),""),
                                fill         = c("blue", "red"),
                                alpha        = 0.3,
                                lty          = "blank",
                                cex          = 2,
                                cat.cex      = 1.5,
                                cat.pos      = c(0, 0),
                                cat.dist     = 0.05,
                               # cat.just     = list(c(-1, -1), c(1, 1)),
                                ext.pos      = 20,
                                ext.dist     = -0.05,
                                ext.length   = 0.85,
                                ext.line.lwd = 5,
                                ext.line.lty = "dashed")
grid.draw(venn.plot)
}
venn()


samples <- names(jvcf)[!names(jvcf) %in% c("CHR", "POS","TYPE", "REF.x","REF.y","ALT.x","ALT.y")]
isec_set <- filter(jvcf, !is.na(REF.x) & !is.na(REF.y) & TYPE == "SNP")

# Pairwise Concordances


p <- t(as.data.frame(sapply(combn(samples,m=2, simplify=F)[0:2000], function(s) {
  u <- select(isec_set, starts_with(s[1]),starts_with(s[2]))
  match <-  nrow(filter(u, (u[[1]] == u[[2]]) ))
  no_match <- nrow(filter(u, u[[1]] != u[[2]]))
  intersect_concordance <- (match)/nrow(isec_set)
  union_concordance <- (match)/nrow(jvcf)
  list(s1=s[1], s2=s[2], match=match,no_match=no_match,intersect_concordance=intersect_concordance,union_concordance=union_concordance)
})))




