library(parallel)
setwd("~/Documents/git/vcfcompare/")
source("helper_fcns.R")

#===========================#
# Data Generation Functions #
#===========================#

get_pair_stats <- function(f1, f2, label="", lab_val, f1_loc=NA) {
  # This function retrieves comparison data for
  # each VCF.
  tmp <- tempfile()
  if (is.na(f1_loc)) {
    f1_loc <- f1
  }
  system(sprintf('bcftools stats -s - %s %s > %s', f1_loc, f2, tmp))
  SN  <- import_table("SN", tmp)
  GCsS <- import_table("GCsS", tmp)
  
  # Set up names.
  ids <- gsub("(.bcf|.vcf|.vcf.gz|.gz)","",list("0"=f1,"1"=f2,"2"=sprintf("%s__%s", f1,f2)))
    
  # Fix SN Group, ID
  SN <- reshape(SN, direction="wide", timevar=c("key"), ids=c("id"))
  names(SN) <- make.names(gsub("(number.of.|value.|:)","",names(SN)))
  
  # Create SNP_Total, Indel_Total Columns
  SN$SNP_Total[SN$id == 0] <- SN[SN$id == 0,"SNPs"] + SN[SN$id == 2,"SNPs"]  
  SN$SNP_Total[SN$id == 1] <- SN[SN$id == 1,"SNPs"] + SN[SN$id == 2,"SNPs"] 
  SN$SNP_Total[SN$id == 2] <- SN[SN$id == 2,"SNPs"] # Intersection!
  SN$Indel_Total[SN$id == 0] <- SN[SN$id == 0,"indels"] + SN[SN$id == 2,"indels"]  
  SN$Indel_Total[SN$id == 1] <- SN[SN$id == 1,"indels"] + SN[SN$id == 2,"indels"] 
  SN$Indel_Total[SN$id == 2] <- SN[SN$id == 2,"indels"]
  

  GCsS <- mutate(GCsS, matches =  RR.Hom.matches + RA.Het.matches + AA.Hom.matches) %.%
  mutate(mismatches =  RR.Hom.mismatches + RA.Het.mismatches + X.10.AA.Hom.mismatches) %.%
  mutate(isec_concordance = matches/(matches+mismatches)) %.%
  mutate(abs_concordance = (matches) / (sum(SN$SNPs)))
  
  # Clean ids
  SN$file <- ids[SN$id+1]
  
  
  # Add group label
  SN$label <- label
  GCsS$label <- label
  SN$lab_val <- as.numeric(lab_val)
  GCsS$lab_val <- as.numeric(lab_val)
  
  list("SN"=SN, "GCsS"=GCsS)
}

parse_query <- function(q_string) {
  # Parses the query String
  direction <- str_extract(q_string,"(>|<|=)")
  items <- unlist(str_split(q_string,"(>|<|=)"))
  filter <- items[1]
  values <- unlist(str_split(items[2],","))
  list(direction=direction,filter=filter,values=values)
}

filter_stats <- function(f1, f2, q_string) {
  q <- parse_query(q_string)
  r <- mclapply(q$values, mc.cores=detectCores(), function(x) { 
      tmp <- tempfile()
      system(sprintf('bcftools filter -O b --include "%s%s%s" %s > %s', q$filter, q$direction, x, f1, tmp))
      system(sprintf('bcftools index %s', tmp))
      get_pair_stats(f1, f2, "f1_loc"=tmp,  sprintf("%s%s%s", q$filter, q$direction, x), x)
  })
  
  # Bind data frames together.
  SN <- do.call(rbind.data.frame,sapply(r, function(x) { x["SN"] }))
  GCsS <- do.call(rbind.data.frame,sapply(r, function(x) { x["GCsS"] }))
  list("SN"=SN, "GCsS"=GCsS, "query"= q)
}

#===========================#
# Plotting Functions        #
#===========================#


f1 <- "04_mmp_strains.bcf"
f2 <- "andersen08_radseq.ws220.bcf"
query <- "DP<800,1000,2000,3000,4000,5000,10000,15000,20000"

f <- filter_stats(f1, f2, query)

# Plot total number of SNPs
ggplot(f$SN[f$SN$id == 0 | f$SN$id == 2 ,]) +
  geom_line( aes(x=lab_val, y=SNP_Total, group=id, color=file)) + 
  labs(title="Great", x=sprintf("%s %s", f$query$filter, f$query$direction ), y="SNPs")

clean_names <- gsub("(.bcf|.vcf|.vcf.gz|.gz)","", c(f1, f2))

# Plot Union Concordance
ggplot(f$GCsS) +
  geom_line( aes(x=lab_val, y=isec_concordance, group=sample, color=sample)) + 
  geom_point( aes(x=lab_val, y=isec_concordance, group=sample, color=sample)) +
  stat_summary(fun.y=mean, mapping = aes(x=lab_val, y = isec_concordance), geom="line", size = 2) +
  labs(title=sprintf("Intersect Concordance: %s - %s", clean_names[1], clean_names[2]), x=sprintf("%s %s", q$filter, q$direction ), y="SNPs")

ggplot(f$GCsS) +
  geom_line( aes(x=lab_val, y=abs_concordance, group=sample, color=sample)) + 
  geom_point( aes(x=lab_val, y=abs_concordance, group=sample, color=sample)) +
  stat_summary(fun.y=mean, mapping = aes(x=lab_val, y = abs_concordance), geom="line", size = 2) +
  labs(title="Union Concordance", x=sprintf("%s %s", q$filter, q$direction ), y="SNPs")
  
  
