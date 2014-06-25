# A library of functions for processing multi-well data
source("./Figures/loadData.R")


# HCT8.txt
# A cell number titration was used to determine the number
# of cells and time to reach maximum impedance.
subset = select(wells, file="HCT8.txt")
plot(subset, xlim=c(0,48), color="concentration", type="total", 
     replicates=TRUE, ID="cellSeed", sd=FALSE)

# Then a concentration-response experiment was performed on the
# wells seeded with 5,000 cells
A = normalize_toxin( select(subset, "HCT8[5000] & TcdA"))
B = normalize_toxin( select(subset, "HCT8[5000] & TcdB"))              
plotA = plot(A, color="concentration", ID="toxinAdd", xlim=c(-1,10) )
plotB = plot(B, color="concentration", ID="toxinAdd", xlim=c(-1,10) ) 
grid.arrange(plotA, plotB, nrow=1)  

# HCT8-2a.txt and HCT8-2b.txt
# Experiments with gdTcdB+TcdA and gdTcdB+TcdB were done
# to investigate any possible competition
subset = normalize_toxin(select(wells, file=c("HCT8-2a.txt","HCT8-2b.txt")), xlim=c(-1,10))
tox_plot = function(x, replicates=TRUE) {
  plot(x, color="concentration", replicates=replicates) + ylim(c(0,1.2))
}
p1 = tox_plot( select(subset, "gdTcdB & !(TcdA | TcdB)")) + ggtitle("gdTcdB")
p2 = tox_plot( select(subset, "TcdA")) + ggtitle("TcdA")
p3 = tox_plot( select(subset, "TcdB")) + ggtitle("TcdB")
grid.arrange(p1, p2, p3, nrow=2)

# HCT8-3.txt
# Similar experiments were repeated at different concentrations
subset = normalize_toxin( select(wells, "TcdA | TcdB | gdTcdB", 
                                 file="HCT8-3.txt"), xlim=c(-1,10) )
p1 = tox_plot( select(subset, "gdTcdB & !(TcdA | TcdB)")) + ggtitle("gdTcdB")
p2 = tox_plot( select(subset, "TcdA")) + ggtitle("TcdA")
p3 = tox_plot( select(subset, "TcdB")) + ggtitle("TcdB")
grid.arrange(p1, p2, p3, nrow=2)

# HCT8-4.txt
# Again, similar experiments were performed
subset = normalize_toxin( select(wells, file="HCT8-4.txt"), xlim=c(-1,12) )
p1 = tox_plot( select(subset, "gdTcdB & !(TcdA | TcdB)")) + ggtitle("gdTcdB")
p2 = tox_plot( select(subset, "TcdA[100]")) + ggtitle("TcdA 100 ng/ml")
p3 = tox_plot( select(subset, "TcdA[1000]")) + ggtitle("TcdA 1000 ng/ml")
p4 = tox_plot( select(subset, "TcdB[10]")) + ggtitle("TcdB 10 ng/ml")
p5 = tox_plot( select(subset, "TcdB[100]")) + ggtitle("TcdB 100 ng/ml")
p6 = tox_plot( select(subset, "TcdB[100]"), replicates=FALSE) + ggtitle("TcdB 100 ng/ml")
grid.arrange(p1, p2, p3, p4, p5, p6, nrow=3)


# CHO.txt
subset = normalize_toxin(select(wells, file="CHO.txt"))
p1 = plot(select(subset, "TcdA"), xlim=c(-1,20), color="concentration") + ggtitle("TcdA")
p2 = plot(select(subset, "TcdB"), xlim=c(-1,10), color="concentration") + ggtitle("TcdB")
grid.arrange(p1, p2, nrow=1)

# IMCE.txt
subset = normalize_toxin(select(wells, "IMCE"))
p1 = plot(select(subset, "TcdA"), color="concentration", xlim=c(-1,20), replicates=TRUE) + ggtitle("TcdA")
p2 = plot(select(subset, "TcdB"), color="concentration", xlim=c(-1,20), replicates=TRUE) + ggtitle("TcdB")
grid.arrange(p1, p2, nrow=1)

# HUVEC-a.txt & HUVEC-b.txt
subset = normalize_toxin(select(wells, "T84"))
p1 = plot(select(subset, "TcdA"), color="concentration", xlim=c(-1,20)) + ggtitle("TcdA")
p2 = plot(select(subset, "TcdB"), color="concentration", xlim=c(-1,20)) + ggtitle("TcdB")
grid.arrange(p1, p2, nrow=1)

# J774-a.txt & J774-b.txt
subset = normalize_toxin(select(wells, file=c("J774-a.txt","J774-b.txt")))
p1 = plot(select(subset,"TcdA"), xlim=c(-1,40), replicates=TRUE, color="concentration") + ggtitle("TcdA")
p2 = plot(select(subset,"TcdB"), xlim=c(-1,40), replicates=TRUE, color="concentration") + ggtitle("TcdB")
grid.arrange(p1, p2, nrow=1)

# J774-2.txt
subset = normalize_toxin(select(wells, file="J774-2.txt"))
plot(subset, replicates=TRUE, xlim=c(-1,24), color="concentration")

# J774-3.txt
subset = normalize_toxin(select(wells, file=c("J774-3a.txt","J774-3b.txt")))
tox_plot = function(x, replicates=TRUE) plot(x, color="concentration", replicates=TRUE, xlim=c(-1,10))
p1 = tox_plot( select(subset, "gdTcdB & !(TcdA | TcdB)")) + ggtitle("gdTcdB")
p2 = tox_plot( select(subset, "TcdA")) + ggtitle("TcdA")
p3 = tox_plot( select(subset, "TcdB")) + ggtitle("TcdB")
grid.arrange(p1, p2, p3, nrow=2)

# J774-4.txt
subset = normalize_toxin(select(wells, file="J774-4.txt"))
tox_plot = function(x, replicates=TRUE) plot(x, color="concentration", replicates=TRUE, xlim=c(-1,24))
p1 = tox_plot( select(subset, "gdTcdB & !(TcdA | TcdB)")) + ggtitle("gdTcdB")
p2 = tox_plot( select(subset, "TcdA")) + ggtitle("TcdA")
p3 = tox_plot( select(subset, "TcdB")) + ggtitle("TcdB")
grid.arrange(p1, p2, p3, nrow=2)

# J774-5.txt
subset = normalize_toxin(select(wells, file="J774-5.txt"))
tox_plot = function(x, replicates=TRUE) plot(x, color="concentration", replicates=TRUE, xlim=c(-1,24))
p1 = tox_plot( select(subset, "gdTcdB & !(TcdA | TcdB)")) + ggtitle("gdTcdB")
p2 = tox_plot( select(subset, "TcdA")) + ggtitle("TcdA")
p3 = tox_plot( select(subset, "TcdB")) + ggtitle("TcdB")
grid.arrange(p1, p2, p3, nrow=2)

# J774-6.txt
subset = normalize_toxin(select(wells, file="J774-6.txt"))
tox_plot = function(x, replicates=TRUE, xlim=c(-1,24), ...) 
  plot(x, color="concentration", replicates=replicates, xlim=xlim, ...)
p1 = tox_plot( select(subset, "gdTcdB & !(TcdA | TcdB)")) + ggtitle("gdTcdB")
p2 = tox_plot( select(subset, "TcdA[1]")) + ggtitle("TcdA 1 ng/ml")
p3 = tox_plot( select(subset, "TcdA & !TcdA[10] & !gdTcdB")) + ggtitle("TcdA")
p4 = tox_plot( select(subset, "TcdB[0.1]"), xlim=c(0,5)) + ggtitle("TcdB 0.1 ng/ml")
p5 = tox_plot( select(subset, "TcdB & !gdTcdB")) + ggtitle("TcdB")
grid.arrange(p1, p2, p3, p4, p5, nrow=3)

## Next is PMNs
# PMN-2a.txt & PMN-2b.txt
subset = select(wells, file=c("PMN-2a.txt", "PMN-2b.txt"))
subset2 = transform(subset, c("tcenter","level"), ID="toxinAdd")
p1 = tox_plot( select(subset2, "TcdA")) + ggtitle("TcdA")
p2 = tox_plot( select(subset2, "TcdB")) + ggtitle("TcdB")
grid.arrange(p1, p2, nrow=1)

# PMN-a.txt & PMN-b.txt
subset = select(wells, file=c("PMN-a.txt","PMN-b.txt"))
subset2 = transform(subset, c("tcenter","level"), ID="toxinAdd")
p1 = tox_plot( select(subset2, "!TcdB")) + ggtitle("TcdA")
p2 = tox_plot( select(subset2, "!TcdA"), xlim=c(-1,10)) + ggtitle("TcdB")
grid.arrange(p1, p2, nrow=1)

# PMN-3.txt
subset = select(wells, file="PMN-3.txt")
subset2 = transform(subset, c("tcenter","level"), ID="toxinAdd")
p1 = tox_plot( select(subset2, "(TcdA & IL8) | (IL8 & !TcdA & !TcdB)")) + ggtitle("TcdA+IL8")
p2 = tox_plot( select(subset2, "TcdA & !IL8")) + ggtitle("TcdA")
p3 = tox_plot( select(subset2, "(TcdB & IL8) | (IL8 & !TcdA & !TcdB)"), replicates=FALSE) + ggtitle("TcdB+IL8")
p4 = tox_plot( select(subset2, "TcdB & !IL8")) + ggtitle("TcdB")
grid.arrange(p1, p2, p3, p4, nrow=2)

# PMN-4.txt
subset = select(wells, file="PMN-4.txt")
subset2 = transform(subset, c("tcenter","level"), ID="toxinAdd")
p1 = tox_plot( select(subset2, "(TcdA & IL8) | (IL8 & !TcdA & !TcdB)")) + ggtitle("TcdA+IL8")
p2 = tox_plot( select(subset2, "TcdA & !IL8")) + ggtitle("TcdA")
p3 = tox_plot( select(subset2, "(TcdB & IL8) | (IL8 & !TcdA & !TcdB)"), replicates=FALSE) + ggtitle("TcdB+IL8")
p4 = tox_plot( select(subset2, "TcdB & !IL8")) + ggtitle("TcdB")
grid.arrange(p1, p2, p3, p4, nrow=2)

# LaTeX tables
fwells = split(wells, filename(wells))
fwells = fwells[sort(names(fwells))]
lapply( fwells, latex_layout )



