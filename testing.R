# A library of functions for processing multi-well data
source("./Scripts2/library.R")
#source("./Scripts2/latexLibrary.R")
wells = parse.RTCAanalyze( metadata="./Annotations.csv", data.dir="./Data" )


# First, need to find the times at which I added the toxin


subset = retrieveWells(wells, file="J774-6.txt" )
tsubset = transform( subset, c("tcenter","normalize"), ID="toxinAdd")



plot( tsubset[1:16], xlim=c(-1,10), color="by.concentrations", replicates=TRUE )

plot( tsubset, xlim=c(-1,10), color="by.concentrations", replicates=TRUE )


# TcdA alone
a = retrieveWells( tsubset, compounds=c("TcdA"), controls=TRUE)
a = a[ !(groupWells(a,group="by.concentrations") %in% "gdTcdB-10.TcdA-1") ]
plot( a, xlim=c(-1,25), color="by.concentrations", replicates=TRUE )

# TcdA and gdTcdB
a = retrieveWells( tsubset, compounds=c("TcdA","gdTcdB"), controls=TRUE)
af = groupWells(a,group="by.concentrations")
a2 = a[ af %in% c("TcdA-1","","gdTcdB-10.TcdA-1","gdTcdB-10") ]
plot( a2, xlim=c(-1,25), color="by.concentrations", replicates=TRUE )

# TcdB alone
a = retrieveWells( tsubset, compounds=c("TcdB"), controls=TRUE)
a = a[ !(groupWells(a,group="by.concentrations") %in% c("gdTcdB-10.TcdB-0.1","gdTcdB-1.TcdB-0.1")) ]
plot( a, xlim=c(-1,25), color="by.concentrations", replicates=TRUE )

# TcdB and gdTcdB
a = retrieveWells( tsubset, compounds=c("TcdB","gdTcdB"), controls=TRUE)
af = groupWells(a,group="by.concentrations")
a2 = a[ af %in% c("gdTcdB-10.TcdB-0.1","gdTcdB-1.TcdB-0.1","gdTcdB-1","gdTcdB-10","TcdB-0.1","") ]
plot( a2, xlim=c(-0.2,5), color="by.concentrations", replicates=TRUE )

# TcdB alone versus gdTcdB alone
a = retrieveWells( tsubset, compounds=c("TcdB","gdTcdB"), controls=TRUE)
af = groupWells(a,group="by.concentrations")
a2 = a[ af %in% c("gdTcdB-1","gdTcdB-10","gdTcdB-100","TcdB-1","TcdB-10","TcdB-100","") ]
plot( a2, xlim=c(-0.2,20), color="by.concentrations", replicates=TRUE )








