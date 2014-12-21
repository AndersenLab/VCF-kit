

histogram = """suppressMessages(library(ggplot2))
suppressMessages(library(data.table))

{opts.functions}

df <- as.data.frame(fread("{filename}.txt"))

ggplot(df) +
  geom_histogram(aes(x=df${var1}),{opts.binwidth}) +
  labs(x="{var1}", y="Count", title="Distribution of {var1}") +
  theme(title=element_text(size=20), 
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15)) {opts.add}

ggsave("{filename}.png", width=18, height=12)

"""

barchart = """suppressMessages(library(ggplot2))
suppressMessages(library(data.table))

{opts.functions}

df <- as.data.frame(fread("{filename}.txt"))

ggplot(df) +
  geom_bar(aes(x=df${var1}),{opts.binwidth}) +
  labs(x="{var1}", y="Count", title="Distribution of {var1}") +
  theme(title=element_text(size=20), 
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15)) {opts.add}

ggsave("{filename}.png", width=18, height=12)

"""



#=====================#
# Specialty Functions #
#=====================#

genetic_scale = """genetic_scale <- function(n) {
  paste0(n/1000000, "Mb")
}"""