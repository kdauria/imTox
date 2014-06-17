rm(list=ls())
detach("package:wellz", unload=TRUE)
library(wellz)
library(reshape)
library(ggplot2)
library(data.table)

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
plot(x,replicates=TRUE, sd=FALSE, color="compound")


############## Reorganize the data
dat = melt_wellList(x)
p = param_matrix(x)
p$A.conc = group(x,"concentration",compound="A")
p$B.conc = group(x,"concentration",compound="B")
p$mB.conc = group(x,"concentration",compound="mB")
p$names = group(x,"compound")
p$tox.conc = p$A.conc + p$B.conc + p$mB.conc*as.numeric(p$names=="mB")

p$comps = p$names
temp = p$names; yn = p$names %in% c("A+mB","B+mB")
p$comps[yn] = paste( p$names[yn], p$mB.conc[yn], sep="." )

dat = dat[p]


############  Plot the data
p = ggplot(dat, aes(x=factor(log10(tox.conc)),y=value,fill=comps, color=comps)) + facet_wrap(~t) + 
  stat_summary(fun.y="mean", geom="bar", position="dodge", alpha=0.4) +
  geom_point(position=position_dodge(width=0.8))
p

############# Plot subsets of the data
dat2 = dat[ dat$tox.conc>1, ]
p %+% dat2

dat3 = dat[ dat$tox.conc<=1, ]
p %+% dat3

########## The plotting is ugly because
# position_dodge drops a bar if it isn't in one of the
# factors on the x-axis. Therefore, I'm going to need to
# calculate the summary statistics beforehand and then do the
# plotting

# What are the keys that separate the data values into groups?
# tox.conc, comps, and t
sum.stats = as.data.frame( dat[,mean(value),by=list(tox.conc,comps,t)] )
sum.stats[order(sum.stats$t,sum.stats$comps),]


groups = with(sum.stats, interaction(comps,tox.conc,t))
length( levels(factor(sum.stats$comps))) # 9
length( levels(factor(sum.stats$tox.conc))) # 6
length( levels(factor(sum.stats$t))) # 3
# Meh this is possibly 9*6*3 = 162 groups. It makes
# more sense to pick and choose what I want to show
# because showing all is just too complicated
















