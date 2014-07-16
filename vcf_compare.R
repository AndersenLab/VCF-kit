library(data.table)


fload_vcf <- function(f) {
  tmp <- tempfile()
  system(sprintf("bcftools query  -f '%%CHROM\t%%POS\t%%REF\t%%ALT\t%%TYPE[\t%%SAMPLE=%%TGT]\n' %s | gawk -f vcf_comp.awk > %s", f, tmp))
  fread(tmp, sep="auto")  
}

vcf1 <- fload_vcf("04_mmp_strains.bcf")
vcf2 <- load_vcf("wgs.vcf.gz")
jvcf <- left_join(vcf1, vcf2, by=c("CHR", "POS", "TYPE"), all=TRUE)


union_samples <- union(names(vcf1)[6:length(vcf1)], names(vcf2)[6:length(vcf2)])
intersect_samples <- intersect(names(vcf1)[6:length(vcf1)], names(vcf2)[6:length(vcf2)])

union_variants_count <- length(vcf1$CHR) + (length(vcf2$CHR) - sum(is.na(jvcf$REF.y) == F))


# Union Samples, Union Variants

# Intersect Samples, Intersect Concordance
for (i in intersect_samples) {
  u <- mutate(u, i.c= )
}
u <- mutate()

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



    

## Union Variants
  


j <- filter(jvcf, REF.x == REF.y)

# paste(rev(str_split("C/T","/")[[1]]), collapse="/")


# Check that REF columns align
sum(z$REF.x == z$REF.y && z$TYPE == "snp",na.rm=TRUE)
sum(z$REF.x != z$REF.y,na.rm=TRUE)
sum(z$REF.x == z$ALT.y,na.rm=TRUE)
sum(z$REF.x != z$ALT.y,na.rm=TRUE)

sum(jvcf$ED3049.x==jvcf$ED3049.y, na.rm=TRUE)
(subset(jvcf, ED3049.x != ED3049.y, select = c("ED3049.x", "ED3049.y")))
