# This script was generated using vcf-toolbox

suppressMessages(library(ggplot2))
suppressMessages(library(data.table))
suppressMessages(library(scales))

# Query {query}

{opts.functions}

df <- as.data.frame(fread("{filename}.txt"))

{opts.filters}

{plot}

ggsave("{filename}.png", width=18, height=12)