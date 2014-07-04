#' Make Figure 4
#' 
#' Makes the different panels for Figure 4 and combines them.
#' The wells object is loaded with \code{load_data} if it is
#' not provided. 
#' 
#' @import wellz
#' @import gridExtra
#' @import ggplot2
#' @import reshape
#' @param wells the output of \code{load_data}
#' @export
figure_4 = function( wells=NULL ) {
  
  require(ggplot2)
  if(is.null(wells)) wells = load_data()

  subset = normalize_toxin(select(wells, file="HCT8-4.txt"))
  p1 = plot( select(subset,"TcdA[100-1000] & !gdTcdB[100] | (gdTcdB[1000] & !TcdB & !TcdA)"), 
             xlim=c(-0.2,12), replicates=TRUE, color="concentration") + ylim(c(0.1,1.05))
  p2 = plot( select(subset,"(TcdB[10] & !gdTcdB[1000]) | (gdTcdB[100] & !TcdB & !TcdA)"), 
             xlim=c(-0.2,4), replicates=TRUE, color="concentration") + ylim(c(0.1,1.05))
  
  subset = normalize_toxin(select(wells, file="J774-6.txt"))
  p3 = plot( select(subset, "TcdA[1] | gdTcdB[10] & !TcdB"),
             xlim=c(-0.2,12), replicates=TRUE, color="concentration") + ylim(c(0.9,1.7))
  p4 = plot( select(subset, "(TcdB[0.1] | gdTcdB[10]) & !TcdA & !gdTcdB[1]"), 
             xlim=c(-0.2,4), replicates=TRUE, color="concentration") + ylim(c(0.9,1.7))
  
  #devSVG(file="./Figures/comp-curves.svg", width=13, height=8)
  grid.arrange(p1, p2, p3, p4, nrow=2)
  #dev.off()
}

