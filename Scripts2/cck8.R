rm(list=ls())
detach("package:wellz", unload=TRUE)
library(wellz)
library(reshape)
library(ggplot2)
library(data.table)
library(stringr)

################### Parse annotations and data
my_parse_fun = function( fpath ) {
  a = read.csv(fpath, sep="\t")
  am = melt(a, measure.vars="A")
  out = cast(am, hour~well )
  out$i = 1:nrow(out)
  colnames(out)[1] = "t"
  out
}

x = parse_metadata(metadata="/Users/kd3jd/Desktop/imTox/Data/CCK8/Annotations.csv",
                   data.dir="/Users/kd3jd/Desktop/imTox/Data/CCK8/",
                   parse_fun = my_parse_fun, spline=TRUE)
plot(x,replicates=TRUE, sd=FALSE, color="compound")

# Normalize to controls and lysis values
controls = cast(melt_wellList( select(x,"control") ), t~., fun.aggregate=mean )[,2]
blanks = cast(melt_wellList( select(x,"blank") ), t~., fun.aggregate=mean )[,2]
for( i in seq_along(x)) vdata(x[[i]]) = (vdata(x[[i]])-blanks)/(controls-blanks)
x = add_spline(x) # need to update the spline whenever the data is modified
plot(x,replicates=TRUE, sd=TRUE, color="compound") + facet_wrap(~compound)

# Select a subset of the data
y = select(x, "A[10-100] | B[10-100] | mB[10-100] | control | lysis")
z = select(y, "!(A & mB) & !(B & mB)")

# Reorganize the data so that bar plots can be made
dat = melt_wellList_params(z,color="concentration")

pdf("./bars.pdf",width=3,height=7)
ggplot(data=dat, aes(x=factor(t),y=value)) +
  stat_summary(fun.y="mean", geom="bar") +
  facet_wrap(~concentration,ncol=2) +
  geom_point(aes(color=concentration)) +
  theme(strip.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_line(colour="lightblue"),
        panel.background=element_blank(),
        legend.position="none",
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.title=element_blank())
dev.off()








