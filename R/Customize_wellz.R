#' Custom normalization
#' 
#' Normalizes impedance curves. By default, \code{ID="toxinAdd"}, so
#' that the data will be centered and normalized to the time at which
#' toxin was added
#' 
#' @import wellz
#' @param x a \code{wellList} object
#' @param ID a \code{character}
#' @param xlim slice the data. Default \code{NULL} means the data isn't slice.
#' @param ... passed to wellz::transform
#' @export
normalize_toxin = function(x, ID="toxinAdd", xlim=NULL, ...) {
  if(!is.null(xlim)) {
    transform(x, c("tcenter","slice","normalize"), ID=ID, xlim=xlim, ... )
  } else {
    transform(x, c("tcenter","normalize"), ID=ID, ... )
  }
}

#' Change default of wellz::select
#'
#' Changes the default arguments of some of the arguments of
#' \code{wellz::select}: \code{controls=TRUE}.
#' 
#' @import wellz
#' @param ... passed to \code{wellz::select}
select = function(...,controls=TRUE) wellz::select(...,controls=controls)

#' small theme for ggplot2
#'
#' @import ggplot2
#' @export
small_theme = function() {
  theme(legend.text = element_text(size=6),
        legend.key.size = unit(6, "points"),
        legend.margin = unit(1, "points"),
        plot.margin = unit(rep(0,4), "points"),
        axis.text = element_text(size=6),
        axis.title = element_text(size=6),
        strip.text = element_text(size=6),
        plot.title = element_text(size=8))
}

#' change ggplot to use small_theme
#' 
#' @import ggplot2
#' @export
ggplot = function(...)
  ggplot2::ggplot(...) + small_theme()

#' change qplot to use small_theme
#' 
#' @import ggplot2
#' @export
qplot = function(...)
  ggplot2::qplot(...) + small_theme()


#' change defaults of wellz::plot.wellList
#' 
#' The new defaults are given in the example call. If 
#' \code{title} is not \code{NULL}, then it must be a string
#' that can be used with \code{ggtitle}
#' 
#' @import ggplot2
#' @import wellz
#' @export
plot.wellList = function(..., color="concentration", replicates=TRUE, 
                         ID="toxinAdd", title=NULL) {
  wellz:::plot.wellList(..., color=color, replicates=replicates, ID=ID) + 
    small_theme() + if(!is.null(title)) ggtitle(title)
}

