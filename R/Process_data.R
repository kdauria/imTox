#' load raw data from text files
#' 
#' This uses the metadata in \code{"extdata/Annotations.csv"}
#' to make a \code{wellList} object. The data files are in the
#' \code{"extdata"} directory of the R package. Use \code{find.package}
#' to find the package directory.
#' 
#' @import wellz
#' @export
load_data = function() {
  base.dir = find.package("imTox")
  parse_metadata( metadata=file.path(base.dir,"extdata/Annotations.csv"), 
                  data.dir=file.path(base.dir,"extdata"),
                  parse_fun=wellz:::parse_rtca)
}

