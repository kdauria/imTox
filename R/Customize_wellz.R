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
