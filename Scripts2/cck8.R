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

#################### Look at the raw data
ggmat = function(wells) {
  groups = list(cond="by.concentrations",toxin="by.compounds")
  mdat = melt.wellList.wgroups(wells,groupings=groups,ID="final")
  mdat = mdat[order(mdat$cond,mdat$time),c("toxin","time","cond","value")]
  mdat = na.omit(mdat)
}

# Take a look at the raw absorbance data for CCK8
p = ggplot(ggmat(wells),aes(x=factor(time),y=value)) + geom_point() + facet_wrap( ~ cond)

subset = retrieveWells(wells,compounds=c("A","lysis","control"))
p %+% ggmat(subset)

subset = retrieveWells(wells,compounds=c("B","lysis","control"))
p %+% ggmat(subset)

subset = retrieveWells(wells,compounds=c("mB","lysis","control"))
p %+% ggmat(subset)


############### Transform the data given the controls
x = ggmat(wells)

# make a replicate # column
allcond = interaction(x$time,x$cond)
x$reps = unlist( mapply( seq, 1, rle(as.character(allcond))$lengths ) )

# cast the data with one column per time point (#rows = #wells)
y = cast(x, cond + reps + toxin ~ time )

# sweep out the lysis buffer controls
vc = 4:6 # value columns
y[,vc] = sweep(y[,vc],2,colMeans( y[y$cond=="lysis-1",vc] ))

# normalize to controls
y[,vc] = sweep( y[,vc], 2, colMeans( y[y$cond=="control-1",vc] ), "/" )

############## Plot the transformed data
z = melt(y,measure.vars=c("1","4","24"))
p %+% z

a = z[ z$toxin %in% c("A","lysis","control"), ]
p %+% a

a = z[ z$toxin %in% c("B","lysis","control"), ]
p %+% a

a = z[ z$toxin %in% c("mB","lysis","control"), ]
p %+% a







