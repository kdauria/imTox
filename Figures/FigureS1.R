# A library of functions for processing multi-well data
source("./Scripts2/library.R")
source("./Scripts2/latexLibrary.R")
wells = parse.RTCAanalyze( metadata="./Annotations.csv", data.dir="./Data" )


subset = retrieveWells(wells, file = c("HCT8-4.txt"))
t.subset = normalize_toxin(subset, xlim = c(-2, Inf))

###### Panel A
conditions = groupWells(t.subset, group = "by.compounds") %in% c("gdTcdB", "")
p1 = plot(t.subset[conditions],xlim=c(-2,10)) + ylim(0, 1.2)