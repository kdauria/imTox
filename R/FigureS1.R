#' Code to make supplemental figure
#' 
#' @import ggplot2
#' @import gridExtra
#' @import wellz
#' @export
figure_s1 = function(wells=NULL) {

  require(ggplot2)
  if(is.null(wells)) wells = load_data()
  
  subset = normalize_toxin(select(wells, file=c("J774-a.txt","J774-b.txt")))
  p1 = plot( select(subset,"TcdA"), xlim=c(-0.01,0.7), points=TRUE, 
             title="TcdA") + ylab("Normalize Impedance")
  p2 = plot( select(subset,"TcdB"), xlim=c(-0.01,1.2), points=TRUE, 
             title="TcdB") + ylab("Normalize Impedance")
  grid.arrange(p1,p2,ncol=2)
}