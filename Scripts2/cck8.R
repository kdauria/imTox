source("./Scripts2/library.R")

################### Parse annotations and data

f = parse.metadata.file(fpath="./Data/CCK8/Annotations.csv")
a = df.to.actions(f$actions)
wells = actions.to.wells( a )
class(wells) = c("wellList","list")

g = groupWells(wells,group="by.concentrations")

# Read in the data file
a = read.csv("./Data/CCK8/Absorbances.csv", sep="\t")
dat = cast(a, hour ~ well, value="A" )

# loop through the column of each data matrix, asssigning data to wells
file = getfiles(wells)[1]
for( j in 2:ncol(dat) ) {
  
  loc = expand.wells( colnames(dat)[j] )
  idx = well.index(wells,loc,file)
  data = as.matrix(as.data.frame(dat[,j,drop=FALSE]))
  df = data.frame( sweep = 1:nrow(data), 
                   time=dat$hour, 
                   values=as.numeric(data))
  if( is.na(idx) ) {
    message(str_c("Well ",location," in file ",file, " is not in the protocol. Well skipped."))
  } else {
    wells[[idx]]$data = df
  }
}

#################### Make data matrix
ggmat = function(wells) {
  groups = list(cond="by.concentrations",toxin="by.compounds")
  mdat = melt.wellList.wgroups(wells,groupings=groups,ID="final")
  mdat = mdat[order(mdat$cond,mdat$time),c("toxin","time","cond","value")]
  mdat = na.omit(mdat)
}

# Take a look at the raw absorbance data for CCK8
ggplot(ggmat(wells),aes(x=factor(time),y=value)) + geom_point() + facet_wrap( ~ cond)

subset = retrieveWells(wells,compounds=c("A","lysis","control"))
ggplot(ggmat(subset),aes(x=factor(time),y=value)) + geom_point() + facet_wrap( ~ cond)

subset = retrieveWells(wells,compounds=c("B","lysis","control"))
ggplot(ggmat(subset),aes(x=factor(time),y=value)) + geom_point() + facet_wrap( ~ cond)

subset = retrieveWells(wells,compounds=c("mB","lysis","control"))
ggplot(ggmat(subset),aes(x=factor(time),y=value)) + geom_point() + facet_wrap( ~ cond)



