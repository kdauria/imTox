# Required packages and functions
source("./Figures/loadData.R")

# The actual code to make the figure
subset = select(wells,"TcdA[1000] & !gdTcdB",file="HCT8-4.txt",controls=TRUE)
t.subset = normalize_toxin( subset, xlim=c(-Inf, 100) )
p1 = plot(t.subset, xlim=c(-43,10), replicates=TRUE, color="concentration" )
p2 = plot(t.subset, xlim=c(-1,10), replicates=TRUE, color="concentration" )
grid.arrange( p1, p2, ncol=2 )

