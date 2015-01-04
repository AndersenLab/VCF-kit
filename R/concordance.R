library(gplots)
library(dplyr)
library(grid)
library(stringr)
library(reshape2)


df$Same_Sample <- factor(df$Same_Sample, levels=c(0,1), labels=c("Same","Diff"))


plot <- ggplot(df, aes(x = Concordance, fill = as.factor(Same_Sample) )) +
  geom_histogram(binwidth = 0.01, position="identity") +
  labs(x = "Concordance Rate", y = "Frequency", title="{opts.title}") +
  theme_bw() +
  xlim(c(0,1)) +
  theme(axis.text=element_text(face="bold", size=15), 
            axis.title=element_text(size=20, lineheight=50),
            axis.title.x=element_text(vjust=-2),
            axis.title.y=element_text(vjust=2),
            plot.margin = unit(c(2,2,2,2), "cm"),
            plot.title=element_text(size=25, vjust=2),
            legend.position = "none") 
 
if (length(unique(df$Same_Sample)) == 2) {{
  plot <- plot + 
    facet_grid(Same_Sample  ~ .) +
    scale_fill_manual(values=c("#999999", "#E69F00"), 
                      name="Comparison\nType",
                      labels=c("Identicle Sample", "Different Sample"))
}}


#==========#
# Heat Map #
#==========#


df_switched <- df
names(df_switched)[4:5] <- c("Sample_j", "Sample_i")

df <- unique(rbind(df, df_switched))


# Convert to Matrix
m <- melt(df, id.vars = c("Sample_i","Sample_j"), measure.vars = c("Concordance"))
df_matrix <- as.matrix(acast(m, m$Sample_i ~ m$Sample_j, m$Concordance))


png(paste0("{filename}",".heatmap.png"),    # create PNG for the heat map        
width = 15*300,        # 5 x 300 pixels
height = 15*300,
res = 600,            # 300 pixels per inch
pointsize = 8)

my_palette <- colorRampPalette(c("blue","green","yellow","red"))(n = 399)
col_breaks = c(seq(0,80,length=100), # blue
seq(80,90,length=100),  # Green
seq(90,98,length=100), # for Yellow
seq(98,100,length=100)) # for Red

heatmap.2(df_matrix*100,
          cellnote = round(df_matrix*100,1),  # same data set for cell labels
          main = "{opts.title}", # heat map title
          notecol="black",      # change font color of cell labels to black
          trace="none", 
          breaks=col_breaks,# turns off trace lines inside the heat map
          margins = c(15,15),     # widens margins around plot
          col=my_palette)
  
dev.off()
  