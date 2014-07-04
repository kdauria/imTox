#' Make Figure 1
#' 
#' @import wellz
#' @import gridExtra
#' @import ggplot2
#' @param wells a \code{wellList} object. If \code{NULL}, then the data is 
#'           is loaded with \code{load_data}.
#' @export
figure_1 = function(wells=NULL) {
  require(ggplot2)
  if(is.null(wells)) wells = load_data()
  
  subset = select(wells, "TcdA[1000] & !gdTcdB", file="HCT8-4.txt")
  t.subset = normalize_toxin( subset, xlim=c(-Inf, 100) )
  p1 = plot(t.subset, xlim=c(-43,10), replicates=TRUE, color="concentration" )
  p2 = plot(t.subset, xlim=c(-1,10), replicates=TRUE, color="concentration" )
  grid.arrange( p1, p2, ncol=2 )
}



