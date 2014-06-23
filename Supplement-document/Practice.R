# Redoing the Supplement to be cleaner and faster
# Now with the wellz package
library(ggplot2)
library(stringr)
library(data.table)
library(gridExtra)
library(wellz)
wells = parse_metadata( metadata="./Annotations2.csv", 
                        data.dir="./Data",
                        parse_fun=wellz:::parse_rtca)
normalize_toxin = function(x, ID="toxinAdd", ...) {
  transform(x, c("tcenter","slice","normalize"), ID=ID, ... )
}

###### References to supplement

# Reference 1
subset = select(wells, file=c("J774-a.txt","J774-b.txt"))
t.subset = normalize_toxin(subset,xlim=c(-2,Inf))
p1 = plot( select(t.subset,"TcdA", controls=TRUE), xlim=c(-0.01,0.7), points=TRUE, 
           color="concentration", replicates=TRUE) +
     ylab("Normalize Impedance") + ggtitle("TcdA")
p2 = plot( select(t.subset,"TcdB", controls=TRUE), xlim=c(-0.01,1.2), points=TRUE, 
           color="concentration", replicates=TRUE) +
     ylab("Normalize Impedance") + ggtitle("TcdB")
grid.arrange(p1,p2,ncol=2)

##### Reproducing Figures

# Figure 1
subset = select(wells,"TcdA[1000] & !gdTcdB",file="HCT8-4.txt",controls=TRUE)
t.subset = normalize_toxin( subset, xlim=c(-Inf, 100) )
p1 = plot(t.subset, xlim=c(-43,10), replicates=TRUE, color="concentration" )
p2 = plot(t.subset, xlim=c(-1,10), replicates=TRUE, color="concentration" )
grid.arrange( p1, p2, ncol=2 )

# Figure 2

# Select the data set used for each cell type
hct8 = select(wells,"HCT8[0-6000]",file="HCT8.txt")
files = list( "CHO.txt", "IMCE.txt", c("HUVEC-a.txt","HUVEC-b.txt"), 
              c("T84-a.txt","T84-b.txt"))
subsets = lapply(files, function(f) select(wells, file = f))
subsets = c(list(hct8), subsets)

# Normalize each subset, slicing different sections of the data (xlims)
# and add smoothers to each subset
xlims = list(c(-0.1, 43), c(-1, Inf), c(-1, Inf), c(-1, 60), c(-1, 27))
n.subsets = Map(normalize_toxin, subsets, xlim = xlims)
s.subsets = Map(add_smoother, n.subsets, method="composite2",
                w=c(2,1,1,3,3),
                global.change=c(1.5,1,1,1,1),
                max.knots=c(20,10,10,12,12))

# Calculate ABC for each subset
calculate_area = function(x, lower, upper, ...) {
  area = area_under_curve(x, lower, upper, ID="toxinAdd")
  area$compound = group(x, "compound")
  area$celltype = group(x,"compound",type="total")
  area$concentration = str_extract( group(x,"concentration", ...), 
                                    "[0-9.e]+[-]?[0-9]?[0-9]?" )
  area$area = area$area - mean(area$area[ area$compound=="" ]) # subtract controls
  area
}
lower = c( 0,  0,  0,  0,  0)
upper = c(43, 40, 80, 80, 27)
areas = Map( calculate_area, s.subsets, lower, upper, ID="toxinAdd" )
areas.df = rbindlist(areas)

# Calculate the maximum area under each 
calculate_max_rate = function(x, ...) {
  maxs = max_rate(x, ...)
  maxs$compound = group(x,"compound")
  maxs$celltype = group(x,"compound",type="total")
  maxs$concentration = str_extract( group(x,"concentration", ID="toxinAdd"), 
                                    "[0-9.e]+[-]?[0-9]?[0-9]?" )
  maxs
}
xlim = list( c(0.1,30), c(0.1,30), c(0.1,30), c(0.2,30), c(0.2,30))
maxs = Map( calculate_max_rate, s.subsets, ID="toxinAdd", 
            xlim=xlim, direction="negative")
maxs.df = rbindlist(maxs)

# Why are these different?
params = max_rate(s.subsets[[4]], xlim=c(0.2,30), group=c("concentration","compound"), 
                  ID="toxinAdd", direction="negative")
maxs[[4]]



######## HCT8 are working well now
x = add_smoother(n.subsets[[1]],method="composite2",
                 w=2,max.knots=20,global.change=1.5)
params = max_rate(x, xlim=c(1,30), group=c("concentration","compound"), 
                  ID="toxinAdd", direction="negative")
params$concentration = as.numeric(str_extract(params$concentration,"[0-9.]+"))
ggplot(params, aes(x=log10(concentration),y=rate,color=compound)) + geom_point()
check_rates(x,params,c(-1,2))



######## CHO is working well now
x = n.subsets[[2]]
for( i in 1:length(x)) {
  x[[i]] = add_smoother(x[[i]], method="composite2", w=1, max.knots=10, global.change=1)
}
params = max_rate(x, xlim=c(0.1,30), group=c("concentration","compound"), 
                  ID="toxinAdd", direction="negative")
params$concentration = as.numeric(str_extract(params$concentration,"[0-9.e]+[-]?[0-9]?[0-9]?"))
ggplot(params, aes(x=log10(concentration),y=rate,color=compound)) + geom_point() +
   facet_wrap(~compound)
check_rates(x,params,c(-1,20))

######## IMCE is working well now
x = n.subsets[[3]]
for( i in 1:length(x)) {
  x[[i]] = add_smoother(x[[i]], method="composite2", w=1, max.knots=10, global.change=1)
}
params = max_rate(x, xlim=c(0.1,30), group=c("concentration","compound"), 
                  ID="toxinAdd", direction="negative")
check_rates(x,params,c(-1,20))

######## HUVEC cells working now
x = n.subsets[[4]]
#x = select( n.subsets[[4]], "TcdB[100-1000] | TcdA[300-1000]", controls=TRUE)
for( i in 1:length(x)) {
  x[[i]] = add_smoother(x[[i]], method="composite2", w=3, max.knots=12, global.change=1)
}
params = max_rate(x, xlim=c(0.2,30), group=c("concentration","compound"), 
                  ID="toxinAdd", direction="negative")
params$concentration = as.numeric(str_extract(params$concentration,"[0-9.e]+[-]?[0-9]?[0-9]?"))
params$concentration = factor(params$concentration)
ggplot(params, aes(x=log10(as.numeric(as.character(concentration))),y=rate,color=compound)) + geom_point() +
  facet_wrap(~compound)
check_rates(x,params,c(0,1))
check_rates(x,params,c(0,5))



######## T84 cells working now
x = n.subsets[[5]]
for( i in 1:length(x)) {
  x[[i]] = add_smoother(x[[i]], method="composite2", w=3, max.knots=12, global.change=1)
}
params = max_rate(x, xlim=c(0.2,30), group=c("concentration","compound"), 
                  ID="toxinAdd", direction="negative")
check_rates(x,params,c(0,5))
check_rates(x,params,c(0,30))







# Here is where the error is. There is a problem
# I think when their is missing data points yet there
# are knots? I'll need to look at this specific case
# in more detail
aa= select(x,"TcdB[500]", ID="toxinAdd")[[1]]
plot( aa, xlim=c(-5,3), 
      points=TRUE, smoother=TRUE )

# Function to look at the maximum rates of each curve
# in the context of plots of the data and derivative of the data
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

aa = wdata( select(x, "TcdB[500]", ID="toxinAdd")[[1]] )
write.table(aa,file="/Users/kd3jd/Desktop/testdatac.txt")



########## Area beneath curve for all cell types
ggplot(areas.df, aes(x=as.numeric(concentration),y=area,color=compound)) + 
  geom_point() + scale_x_log10() +
  geom_hline( data=areas.df[areas.df$compound==""], 
              aes(yintercept=area,color=compound), linetype="dashed" ) + # controls
  facet_wrap(~celltype, scale="free") +
  xlab("Concentration (ng/ml)")

########## MaxS for all cell types
ggplot(maxs.df, aes(x=as.numeric(concentration),y=rate,color=compound)) + 
  geom_point() + scale_x_log10() +
  geom_hline( data=maxs.df[maxs.df$compound==""], 
              aes(yintercept=rate,color=compound), linetype="dashed" ) + # controls
  facet_wrap(~celltype, scale="free") +
  xlab("Concentration (ng/ml)")









######### Panel A. Different cell types. Same concentrations. #########
subset = retrieveWells(allwells, compounds = "TcdA", ID = "toxinAdd", 
                       max.concentrations = 101, min.concentrations = 99)
panelA = plot(subset, xlim = c(-1, 10), se = FALSE, color = "by.total.compounds",
              linetype = "by.compounds")

######### Panel B. Same cell type. Different toxins ############
subset = retrieveWells(allwells, compounds = "IMCE")
subset2 = retrieveWells(subset, compounds = c("TcdA", "TcdB"), ID = "toxinAdd",
                        max.concentrations = c(101, 101), min.concentrations = c(99, 99))
panelB = plot(subset2, xlim = c(-0.1, 2), se = FALSE)

######## Panel C. MaxS and ABC for IMCE cells ###########
subset = retrieveWells(allwells, compounds = "IMCE")
MaxS = groupMetric(subset, ID = "toxinAdd", metric = "max.rate")
p1 = plotMetric(MaxS)
ABC = groupMetric(subset, ID = "toxinAdd", metric = "integral")
p2 = plotMetric(ABC)
panelC = arrangeGrob(p1 + theme(legend.position = "none"), 
                     p2 + theme(legend.position = "none"), ncol = 2)

####### Panel D. The MCC for each cell type ##############
# The MCC was found by plotting the ABC of each cell type over
# a range of concentrations. The concentration diverging from
# controls was considered the MCC
ABC = groupMetric(allwells, ID = "toxinAdd", metric = "integral")
ABC$cells = groupWells(allwells, group = "by.total.compounds")
plotMetric(ABC, file = FALSE) + facet_wrap(~cells, scales = "free_y")