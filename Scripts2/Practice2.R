library(DierckxSpline)
library(ggplot2)
x = read.table(file="/Users/kd3jd/Desktop/testdata.txt")
x = read.table(file="/Users/kd3jd/Desktop/testdatab.txt")

y = x[,3]
x = x[,2]


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
    n.knots = max(1, round(relative.change*max.knots))
    
    # Fit a smoother with that many knots
    n.curfit.knots = n.knots + 8 # There must be (k+1)*2 additional knots (I think???)
                      # where k is the degree of the spline (here k=3)
    fit = curfit(x=x.bin,y=y.bin,n=n.curfit.knots)
    bin.knot.locations[[bin]] = knots(fit)  
  }
  
  all.knots = unlist(bin.knot.locations)
  fit = curfit(x,y,knots=all.knots)
  fit
}

fit = composite2(x, y, w=2, max.knots=20, global.change=1.5)


ggplot(data.frame(x,y),aes(x=x,y=y)) + geom_point() +
  geom_line(data=data.frame(x,y=predict(fit,x)),color="blue",size=1.4) +
  geom_point(data=data.frame(z=knots(fit)), aes(x=z,y=1)) +
  xlim(c(0,5))








