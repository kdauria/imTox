#' Make Figure 5
#' 
#' Makes the different panels for Figure 5 and combines them.
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
#' @import grid
#' @import data.table
#' @param wells the output of \code{load_data}
#' @export
figure_5 = function( wells=NULL ) {
  
  require(ggplot2)
  if(is.null(wells)) wells = load_data()
  
  ### Panel A
  subset = normalize_toxin(select(wells, file=c("J774-6.txt")))
  subsetgB = select(subset, "(gdTcdB[1-100] | TcdB[10-100]) & !((TcdB | TcdA) & gdTcdB)")
  top.plot = plot(subsetgB, xlim=c(-1,25), title="gdTcdB") + ylim(c(0,2.9))
  
  ### Panel B
  x = load_cck8_data()
  
  # Normalize to controls and lysis values
  controls = cast(melt_wellList( select(x,"control") ), t~., fun.aggregate=mean )[,2]
  blanks = cast(melt_wellList( select(x,"blank") ), t~., fun.aggregate=mean )[,2]
  for( i in seq_along(x)) vdata(x[[i]]) = (vdata(x[[i]])-blanks)/(controls-blanks)
  x = add_spline(x) # need to update the spline whenever the data is modified
  
  # Select a subset of the data and reorganize it for bar plots
  y = select(x, "(mB[10-100] | control | lysis) & !(A & mB) & !(B & mB)")
  dat = melt_wellList_params(y,color="concentration")
  
  bar.plots = ggplot(data=dat, aes(x=factor(t),y=value)) +
    stat_summary(fun.y="mean", geom="bar") +
    facet_wrap(~concentration,nrow=1) +
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

  grid.arrange(top.plot, bar.plots, nrow=1)
}