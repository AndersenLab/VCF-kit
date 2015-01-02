#=====================#
# Specialty Functions #
#=====================#

genetic_scale = """genetic_scale <- function(n) {
  paste0(n/1000000, "Mb")
}"""

precision_axis = """precision_axis <- function(x) {
  format(x,nsmall = 2,scientific = FALSE)
}"""