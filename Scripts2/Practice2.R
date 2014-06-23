library(DierckxSpline)
library(ggplot2)
x = read.table(file="/Users/kd3jd/Desktop/testdata.txt")
x = read.table(file="/Users/kd3jd/Desktop/testdatab.txt")
x = read.table(file="/Users/kd3jd/Desktop/testdatac.txt")
x = read.table(file="/Users/kd3jd/Desktop/testdatad.txt")

y = x[,3]
x = x[,2]
w=2
max.knots=12
global.change=1.2

composite2 = function(x, y, w, max.knots, global.change) {
  
  # Prepare bins
  n.bins = ceiling(diff(range(x))/w)
  bins = as.numeric(cut(x,n.bins))
  bin.knot.locations = vector(length=n.bins, mode="list")
  
  for( bin in unique(bins) ) {
    
    # subset the data for the bin
    x.bin = x[bins==bin]
    y.bin = y[bins==bin]
    
    # Calculate the number of knots depending on local change
    local.change = diff(range(y.bin))
    relative.change = local.change/global.change
    n.knots = round(relative.change*max.knots)
    
    if(n.knots<=1) {
      bin.knot.locations[[bin]] = median(x.bin)
    } else {
            
      fit.local = curfit.free.knot(x.bin, y.bin, g=n.knots )
      local.knots = knots(fit.local)
      
#       fit.local = smooth.spline(x.bin, y.bin, nknots=n.knots)$fit
#       local.knots = unique( fit.local$knot*fit.local$range + fit.local$min )
      
      local.knots = c(min(x.bin),local.knots,max(x.bin))
      bin.knot.locations[[bin]] = unique(local.knots)
    }
  }
  
  #all.knots = c(min(x),unlist(bin.knot.locations),max(x))
  all.knots = unique(c(min(x),unlist(bin.knot.locations),max(x)))
  fit = curfit(x,y,knots=all.knots)
  fit
}

fit = composite2(x, y, w=2, max.knots=20, global.change=1.5)


newx = wellz:::na_interp( wellz:::insert_na_between(x,n=10) )
ggplot(data.frame(x,y),aes(x=x,y=y)) + geom_point() +
  geom_line(data=data.frame(x=newx,y=predict(fit,newx)),color="blue",size=1.4) +
  geom_point(data=data.frame(z=knots(fit)), aes(x=z,y=1)) +
  xlim(c(-0.1,2))

ggplot(data.frame(x,y),aes(x=x,y=y)) + geom_point() +
  geom_line(data=data.frame(x=newx,y=deriv(fit,newx)),color="blue",size=1.4) +
  geom_point(data=data.frame(z=knots(fit)), aes(x=z,y=1)) +
  xlim(c(-0.1,2))


table(as.character(cut( knots(fit), x)))





