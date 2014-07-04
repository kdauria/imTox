#' Make Figure 3
#' 
#' Makes the different panels for Figure 3 and combines them.
#' The wells object is loaded with \code{load_data} if it is
#' not provided. Panel B was made manually by piecing together
#' different image files in Adobe Illustrator. The data for the
#' bar plots is in the \code{"./inst/extdata/CCK8"} directory and
#' the annotations of the wells is in \code{"./inst/extdata/CCK8/Annotations.csv"}.
#' 
#' @import wellz
#' @import gridExtra
#' @import ggplot2
#' @import reshape
#' @param wells the output of \code{load_data}
#' @export
figure_3 = function( wells=NULL ) {
  
  require(ggplot2)
  if(is.null(wells)) wells = load_data()

  ###### Panel A
  subset = normalize_toxin(select(wells, file=c("J774-a.txt","J774-b.txt")))
  subsetA = select(subset,"!TcdA[300] & !TcdA[3] & !TcdB & !Other")
  subsetB = select(subset,"!TcdB[300] & !TcdB[3] & !TcdA & !Other")
  p1 = plot(subsetA, replicates=TRUE, xlim=c(-1,25), color="concentration") + ggtitle("TcdA") + ylim(c(0,2.9))
  p2 = plot(subsetB, replicates=TRUE, xlim=c(-1,25), color="concentration") + ggtitle("TcdB") + ylim(c(0,2.9))
  
  subset = normalize_toxin(select(wells, file=c("J774-6.txt")))
  subsetgB = select(subset, "(gdTcdB[1-100] | TcdB[10-100]) & !((TcdB | TcdA) & gdTcdB)")
  p3 = plot(subsetgB, replicates=TRUE, xlim=c(-1,25), color="concentration") + ggtitle("gdTcdB") + ylim(c(0,2.9))
  
  #devSVG(file="Figures/J774-curves.svg",width=16,height=4)
  ptop = grid.arrange(p1,p2,p3,nrow=1)
  print(ptop)
  #dev.off()
  
  ### Panel C
  x = load_cck8_data()
  
  # Normalize to controls and lysis values
  controls = cast(melt_wellList( select(x,"control") ), t~., fun.aggregate=mean )[,2]
  blanks = cast(melt_wellList( select(x,"blank") ), t~., fun.aggregate=mean )[,2]
  for( i in seq_along(x)) vdata(x[[i]]) = (vdata(x[[i]])-blanks)/(controls-blanks)
  x = add_spline(x) # need to update the spline whenever the data is modified
  plot(x,replicates=TRUE, sd=TRUE, color="compound") + facet_wrap(~compound)
  
  # Select a subset of the data and reorganize it for bar plots
  y = select(x, "(A[10-100] | B[10-100] | mB[10-100] | control | lysis) & !(A & mB) & !(B & mB)")
  dat = melt_wellList_params(y,color="concentration")
  
  #devSVG(file="./Figures/bars.svg",width=9,height=2)
  bar.plots = ggplot(data=dat, aes(x=factor(t),y=value)) +
    stat_summary(fun.y="mean", geom="bar") +
    facet_wrap(~concentration,nrow=2) +
    geom_point(aes(color=concentration)) +
    theme(strip.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major.x=element_blank(),
          panel.grid.minor.y=element_blank(),
          panel.grid.major.y=element_line(colour="lightblue"),
          panel.background=element_blank(),
          legend.position="none",
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          axis.title=element_blank())
  bar.plots
  #dev.off()
}




