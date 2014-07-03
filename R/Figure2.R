#' Make Figure 2
#' 
#' Makes the different panels for Figure 2 and combines them.
#' The wells object is loaded with \code{load_data} if it is
#' not provided. Before the figure is printed, two multi-panel
#' plots are also printed showing the area between curves (ABC)
#' and maximum slow (MaxS) for different concentrations and cell types.
#' 
#' @import wellz
#' @import gridExtra
#' @import ggplot2
#' @param wells the output of \code{load_data}
#' @export
figure_2 = function( wells=NULL ) {

  if(is.null(wells)) wells = load_data()
  
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
  subset = select(allwells, "TcdA[100]", ID="toxinAdd")
  p1 = plot(subset, xlim=c(-1,10), replicates=TRUE, sd=FALSE, 
            color="compound", type="total")
  
  ######### Panel B. Same cell type. Different toxins ############
  subset = select(allwells, "IMCE & (TcdA[100] | TcdB[100])")
  p2 = plot(subset, xlim = c(-0.1, 2), replicates=TRUE, sd = FALSE, color="compound")
  
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
  p3a = ggplot( maxs.df[celltype=="IMCE" & compound!="",],
                aes(x=as.numeric(concentration), y=rate, color=compound)) + 
        geom_hline(data=maxs.df[celltype=="IMCE" & compound=="",], 
               aes(intercept=rate), color="blue", linetype="dashed") +
        geom_point() + scale_x_log10()
    
  p3b = ggplot( areas.df[celltype=="IMCE" & compound!="",],
                aes(x=as.numeric(concentration), y=area, color=compound)) + 
        geom_hline(data=areas.df[celltype=="IMCE" & compound=="",], 
               aes(intercept=area), color="blue", linetype="dashed") + 
        geom_point() + scale_x_log10()
  p3 = arrangeGrob(p3a, p3b, nrow=2)
  
  ####### Panel D. The MCC for each cell type ##############
  # The MCC was found by plotting the ABC of each cell type over
  # a range of concentrations. The concentration diverging from
  # controls was considered the MCC
  
  ########## Area beneath curve for all cell types
  area.plot = ggplot(areas.df, aes(x=as.numeric(concentration),y=area,color=compound)) + 
    geom_point() + scale_x_log10() +
    geom_hline( data=areas.df[areas.df$compound==""], 
                aes(yintercept=area,color=compound), linetype="dashed" ) + # controls
    facet_wrap(~celltype, scale="free") +
    xlab("Concentration (ng/ml)") + ggtitle("Area Between Curves (ABC)")
  print(area.plot)
  
  ########## MaxS for all cell types
  maxs.plot = ggplot(maxs.df, aes(x=as.numeric(concentration),y=rate,color=compound)) + 
    geom_point() + scale_x_log10() +
    geom_hline( data=maxs.df[maxs.df$compound==""], 
                aes(yintercept=rate,color=compound), linetype="dashed" ) + # controls
    facet_wrap(~celltype, scale="free") +
    xlab("Concentration (ng/ml)") + ggtitle("Maximum Slope (MaxS)")
  print(maxs.plot)
  
  # Manually enter the MCC for each cell type
  d = data.frame( type = c( "CHO", "HCT8", "HUVEC", "IMCE", "T84" ),
                  a    = c( 1    ,  1    ,  1     ,  0.1  , 0.1   ),
                  b    = c( .001 ,  .01  ,  .01  ,  .0001 , 0.1   ) )
  p4 = ggplot( d, aes( x=log10(a), y=log10(b), label=type) ) + geom_text() +
    scale_x_reverse(limits=c(0,-4)) + scale_y_reverse(limits=c(0,-4)) + 
    geom_abline(slope=1,intercept=0) + coord_equal() + geom_point(color="red")
  
  # Show the final plot
  grid.arrange( p1, p2, p3, p4, ncol=2 )
}



