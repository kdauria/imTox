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
A = normalize_toxin( select(subset, "HCT8[5000] & TcdA", controls=TRUE ) )
B = normalize_toxin( select(subset, "HCT8[5000] & TcdB", controls=TRUE ) )              
plotA = plot(A, color="concentration", ID="toxinAdd", xlim=c(-1,10) )
plotB = plot(B, color="concentration", ID="toxinAdd", xlim=c(-1,10) ) 
grid.arrange(plotA, plotB, nrow=1)  

# HCT8-2a.txt and HCT8-2b.txt
# Experiments with gdTcdB+TcdA and gdTcdB+TcdB were done
# to investigate any possible competition
subset = normalize_toxin(select(wells, file=c("HCT8-2a.txt","HCT8-2b.txt")), xlim=c(-1,10))
p1 = plot( select(subset, "gdTcdB & !(TcdA | TcdB)", controls=TRUE), 
           color="concentration", replicates=TRUE ) + ylim(c(0,1.2)) + ggtitle("gdTcdB")
p2 = plot( select(subset, "TcdA", controls=TRUE), 
           color="concentration", replicates=TRUE ) + ylim(c(0,1.2)) + ggtitle("TcdA")
p3 = plot( select(subset, "TcdB", controls=TRUE), 
           color="concentration", replicates=TRUE ) + ylim(c(0,1.2)) + ggtitle("TcdB")
grid.arrange(p1, p2, p3, nrow=2)



p1 = plot(subset, color="concentration", xlim=c(-1,4), replicates=TRUE)

# The maximum rate of the curves
subset = add_smoother( subset, method="composite2", w=2, max.knots=15, global.change=1.2 )
pgon = rbind( c(0,-10), c(0,0.85), c(2.5,0.85), c(2.5,100), c(100,100), c(100,-10), c(0,-10) )
maxs = calculate_max_rate( subset, xlim=c(0.2,30), direction="negative", pgon=pgon)
maxs$concentration = group(subset,"concentration")
p2 = ggplot(maxs,aes(x=concentration,y=rate, color=concentration)) + geom_point()
grid.arrange(p1,p2,nrow=1)

# HCT8-3.txt
subset = normalize_toxin( select(wells, "TcdA | TcdB | gdTcdB", 
                                 file="HCT8-3.txt", controls=TRUE), xlim=c(-2,10) )
p1 = plot( select(subset, "gdTcdB & !(TcdA | TcdB)", controls=TRUE), 
           color="concentration", replicates=TRUE ) + ylim(c(0,1.2)) + ggtitle("gdTcdB")
p2 = plot( select(subset, "TcdA", controls=TRUE), 
           color="concentration", replicates=TRUE ) + ylim(c(0,1.2)) + ggtitle("TcdA")
p3 = plot( select(subset, "TcdB", controls=TRUE), 
           color="concentration", replicates=TRUE ) + ylim(c(0,1.2)) + ggtitle("TcdB")
grid.arrange(p1, p2, p3, nrow=2)






x = select( subset, "gdTcdB[100-1000] & !(TcdA | TcdB)", controls=TRUE)
pgon = rbind( c(0,-10), c(0,0.85), c(2.5,0.85), c(2.5,100), c(100,100), c(100,-10), c(0,-10) )
params = calculate_max_rate( x, xlim=c(0.2,30), direction="negative", pgon=pgon)
check_rates(x, params, c(0,10))

check_rates = function(x, params, xlim ) {
  p1 = plot(x,color="concentration",shape="compound",
            ID="toxinAdd", xlim=xlim, points=TRUE, smoother=TRUE) + 
    geom_point(data=params, size=5) + 
    geom_point(data=params, size=3, color="white") +
    facet_wrap(~compound) + xlim(xlim)
  p2 = plot(x,color="concentration",shape="compound", deriv=1,
            ID="toxinAdd", xlim=xlim, points=FALSE, smoother=TRUE) + 
    geom_point(data=params,aes(y=rate-0.1*diff(range(params$rate))),size=5) + 
    geom_point(data=params,aes(y=rate-0.1*diff(range(params$rate))),size=3, color="white") +
    facet_wrap(~compound) + xlim(xlim)
  grid.arrange(p1,p2,nrow=2)
}