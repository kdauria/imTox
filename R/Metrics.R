#' Area between curves (ABC)
#' 
#' Calculate the area under each curve and return a data
#' frame with the metadata for each well (i.e., the compounds,
#' cell type, and concentration of compounds in the well). This
#' is essentially a wrapper around \code{wellz::area_under_curve}
#' 
#' @import wellz
#' @import stringr
#' @param x a \code{wellList} object
#' @param lower the lower limit of the integral
#' @param xlim slice the data. Default \code{NULL} means the data isn't slice.
#' @param ... passed to wellz::group
#' @export
calculate_area = function(x, lower, upper, ...) {
  area = area_under_curve(x, lower, upper, ID="toxinAdd")
  area$compound = group(x, "compound")
  area$celltype = group(x,"compound",type="total")
  area$concentration = str_extract( group(x,"concentration", ...), 
                                    "[0-9.e]+[-]?[0-9]?[0-9]?" )
  area$area = area$area - mean(area$area[ area$compound=="" ]) # subtract controls
  area
}


#' Maximum slope of curve (MaxS)
#' 
#' Calculate the maximum slope of a curve using \code{wellz::max_rate}.
#' The results are returned in a data frame with the rate for each 
#' well and the well metadata in columns (i.e., the compound, cell type,
#' and concentration of compounds in each well).
#' 
#' @import wellz
#' @import stringr
#' @param x a \code{wellList} object
#' @param ... passed to wellz::max_rate
#' @export
calculate_max_rate = function(x, ...) {
  maxs = max_rate(x, ...)
  maxs$compound = group(x,"compound")
  maxs$celltype = group(x,"compound",type="total")
  maxs$concentration = str_extract( group(x,"concentration", ID="toxinAdd"), 
                                    "[0-9.e]+[-]?[0-9]?[0-9]?" )
  maxs
}


