# A library of functions for processing multi-well data
source("./Figures/loadData.R")

# Select the data set used for each cell type
hct8 = select(wells,"HCT8[0-6000]",file="HCT8.txt")
files = list( "CHO.txt", "IMCE.txt", c("HUVEC-a.txt","HUVEC-b.txt"), 
              c("T84-a.txt","T84-b.txt"))
subsets = lapply(files, function(f) select(wells, file = f))
subsets = c(list(hct8), subsets)

# Normalize each subset, slicing different sections of the data (xlims)
# and add smoothers to each subset
xlims = list(c(-1, 43), c(-1, Inf), c(-1, Inf), c(-1, 60), c(-1, 27))
n.subsets = Map(normalize_toxin, subsets, xlim = xlims)
s.subsets = Map(add_smoother, n.subsets, method="composite2",
                w=c(2,1,1,3,3),
                global.change=c(1.5,1,1,1,1),
                max.knots=c(20,10,10,12,12))

######### Panel A. Different cell types. Same concentrations. #########
allwells = do.call("c",s.subsets)
subset = select(allwells, "TcdA[100]", ID="toxinAdd", controls=TRUE)
plot(subset, xlim=c(-1,10), replicates=TRUE, sd=FALSE, 
     color="compound", type="total")

######### Panel B. Same cell type. Different toxins ############
subset = select(allwells, "IMCE & (TcdA[100] | TcdB[100])", controls=TRUE)
plot(subset, xlim = c(-0.1, 2), replicates=TRUE, sd = FALSE, color="compound")

######## Panel C. MaxS and ABC for IMCE cells ###########
# Calculate ABC for each subset
lower = c( 0,  0,  0,  0,  0)
upper = c(43, 40, 80, 80, 27)
areas = Map( calculate_area, s.subsets, lower, upper, ID="toxinAdd" )
areas.df = rbindlist(areas)

# Calculate MaxS for each subset
xlim = list( c(0.1,30), c(0.1,30), c(0.1,30), c(0.2,30), c(0.2,30))
maxs = Map( calculate_max_rate, s.subsets, ID="toxinAdd", 
            xlim=xlim, direction="negative")
maxs.df = rbindlist(maxs)

#### Make the plots for just IMCE cells
ggplot( maxs.df[celltype=="IMCE" & compound!="",],
        aes(x=as.numeric(concentration), y=rate, color=compound)) + 
  geom_hline(data=maxs.df[celltype=="IMCE" & compound=="",], 
             aes(intercept=rate), color="blue", linetype="dashed") +
  geom_point() + scale_x_log10()
  
ggplot( areas.df[celltype=="IMCE" & compound!="",],
          aes(x=as.numeric(concentration), y=area, color=compound)) + 
  geom_hline(data=areas.df[celltype=="IMCE" & compound=="",], 
             aes(intercept=area), color="blue", linetype="dashed") + 
  geom_point() + scale_x_log10()
  

####### Panel D. The MCC for each cell type ##############
# The MCC was found by plotting the ABC of each cell type over
# a range of concentrations. The concentration diverging from
# controls was considered the MCC

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




