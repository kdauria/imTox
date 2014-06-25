library(ggplot2)
library(stringr)
library(data.table)
library(gridExtra)
library(wellz)
wells = parse_metadata( metadata="./Annotations2.csv", 
                        data.dir="./Data",
                        parse_fun=wellz:::parse_rtca)

normalize_toxin = function(x, ID="toxinAdd", xlim=NULL, ...) {
  if(!is.null(xlim)) {
    transform(x, c("tcenter","slice","normalize"), ID=ID, xlim=xlim, ... )
  } else {
    transform(x, c("tcenter","normalize"), ID=ID, ... )
  }
}

calculate_area = function(x, lower, upper, ...) {
  area = area_under_curve(x, lower, upper, ID="toxinAdd")
  area$compound = group(x, "compound")
  area$celltype = group(x,"compound",type="total")
  area$concentration = str_extract( group(x,"concentration", ...), 
                                    "[0-9.e]+[-]?[0-9]?[0-9]?" )
  area$area = area$area - mean(area$area[ area$compound=="" ]) # subtract controls
  area
}

calculate_max_rate = function(x, ...) {
  maxs = max_rate(x, ...)
  maxs$compound = group(x,"compound")
  maxs$celltype = group(x,"compound",type="total")
  maxs$concentration = str_extract( group(x,"concentration", ID="toxinAdd"), 
                                    "[0-9.e]+[-]?[0-9]?[0-9]?" )
  maxs
}

select = function(...) wellz::select(...,controls=TRUE)
