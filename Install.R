##To intall Bi_EB function
install.packages("devtools")
library("devtools")

library(roxygen2)

setwd("..parent_directory")
create("Bi.EB")

setwd("./Bi.EB")
document()

setwd("..")
install("Bi.EB")
