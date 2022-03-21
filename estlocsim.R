require(circular)
require(secr)

ncall.mean = 2
ncall = 0:10
pcall = dpois(ncall,ncall.mean)
plot(ncall,pcall)
segments(ncall,rep(0,length(ncall)),ncall,pcall)

traplocs = data.frame(x=c(0,250,500), y=c(0,433,0)) # triangular array
traplocs = data.frame(x=c(-500,0,500), y=c(0,0,0)) # triangular array
posts = read.traps(data=traplocs, detector="proximity", binary.usage=FALSE)
buffer = 2500
mask = make.mask(posts, buffer=buffer, type="trapbuffer")
plot(mask,border=0)
plot(posts, border=0, gridlines=FALSE,add=TRUE)

seed = 58
pop = sim.popn(2/100, posts, buffer=buffer, model2D="even", seed=seed)
plot(pop,add=TRUE)

g0 = 1
sigma = 700
ncalls = 10
ch = sim.capthist(posts, pop, detectpar=list(g0=g0,sigma=sigma),noccasions=ncalls, renumber=FALSE, seed=seed)

plotdets = function(ch,traps,pop,animals=NULL) {
  detected = as.integer(row.names(ch))
  if(is.null(animals)) animals = 1:length(detected)
  for(i in animals) for(j in 1:dim(ch)[3]) for(k in 1:dim(ch)[2]) {
    points(pop$x[detected[i]],pop$y[detected[i]], col=i)
    if(ch[i,k,j]>0) segments(pop$x[detected[i]],pop$y[detected[i]],traps$x[j],traps$y[j],col=i)
  }
}

getangle = function(start,end) as.numeric(atan((start$y-end$y)/(start$x-end$x)))
  
get.intersection = function(post, theta) {
  require(PlaneGeometry)
  line1 <- Line$new(A=c(post$x[1],post$y[1]), B=c(post$x[1]+100*cos(theta[1]),post$y[1]+100*sin(theta[1])))
  line2 <- Line$new(A=c(post$x[2],post$y[2]), B=c(post$x[2]+100*cos(theta[2]),post$y[2]+100*sin(theta[2])))
  pt = intersectionLineLine(line1, line2)
  return(data.frame(x=pt[1],y=pt[2]))
}

plot(mask,border=0)
plot(posts, border=0, gridlines=FALSE,add=TRUE)

sd.angle = 4*pi/180 # angle error std dev, in radians
nsim = 20 # number of random errors in angle from which to calculate mean estimated location

detected = as.integer(row.names(ch)) # get indices of detected animals
# plot the locations of the detected animals:
text(pop[detected,]$x,pop[detected,]$y,labels=as.character(row.names(pop[detected,])),cex=0.5)

ndets = length(detected) # number of animals detected
maxposts = nrow(posts)
maxpairs = choose(maxposts,2) # maximum number of pairs of recorded angles
estlocs = array(rep(NA,ndets*maxpairs*nsim*2), dim=c(ndets,maxpairs,nsim,2), 
                dimnames=list(animal=row.names(pop[detected,]), pair=1:maxpairs, sim=1:nsim, coord=1:2))
keepanimal = rep(FALSE,ndets) # which animals detected by more than one post
for(i in 1:ndets) { # for each animal
  end = data.frame(x=pop$x[detected[i]],y=pop$y[detected[i]]) # animal location
  detfrom = which(apply(ch[i,,],2,sum)>0) # posts from which it was detected
  npostdet = length(detfrom) # number of posts that detected it
  if(npostdet > 1) { # only consider animals detected from at least 2 posts
    keepanimal[i] = TRUE
    start = data.frame(x=posts$x[detfrom[1]],y=posts$y[detfrom[1]]) # location of first post
    theta = getangle(start[1,],end) # angle from first post to animal
    for(post in 2:npostdet) { # get locations and detection angles for all other posts
      start = rbind(start,
                    data.frame(x=posts$x[detfrom[post]],
                               y=posts$y[detfrom[post]])) # location of post
      theta = c(theta,getangle(start[post,],end)) # angle from post to animal
    }
    # go through all pairs of detections of animal
    npts = choose(npostdet,2) # number of pairs
    pno = 0 # pair counter
    errs = matrix(rnorm(nsim*npostdet,0,sd.angle),nrow=nsim) # matrix of nsim angle errors for each post 
    for(j in 1:(npostdet-1)) for(k in (j+1):npostdet) { # go through all pairs of angles (pairs of posts)
      thetas = matrix(c(theta[j]+errs[,j],theta[k]+errs[,k]),nrow=nsim) # estimated angles
      asign = sign(theta[j]-theta[k]) # sign of difference between true angles
      rsign = sign(thetas[,1]-thetas[,2]) # sign of difference between estimated angles
      keep = (asign*rsign) > 0 # keep only estimated angles that converge on same side as true angles
      thetas = thetas[keep,,drop=FALSE] # kept angles
      nkeep = sum(keep) # number of plausible simulated locations
      if(nkeep>0) { # if have any valid estimated locations
        pno = pno+1 # increment pair counter
        for(l in 1:nkeep)  estlocs[i,pno,l,] = as.numeric(get.intersection(start[c(j,k),],thetas[l,])) # estimated locations of animal
      }
    }
  }

}


xlim = range(c(posts$x+buffer,posts$x-buffer))
ylim = range(c(posts$y+buffer,posts$y-buffer))
plot(posts$x,posts$y,xlim=xlim,ylim=ylim,xlab="UTMX",ylab="UTMY",pch="*")
text(pop[detected,]$x,pop[detected,]$y,labels=as.character(row.names(pop[detected,])),cex=0.5)

quartz()
plot(posts$x,posts$y,xlim=xlim,ylim=ylim,xlab="UTMX",ylab="UTMY",pch="*",xaxt="n",yaxt="n")
gotpos = apply(is.na(estlocs),1,sum)-(nsim*maxpairs*2)<0 # which animals detected by at least 2 posts
colr = 0
for(ind in 1:ndets) {
  if(gotpos[ind]) {
    colr = colr+1
    if(colr>9) {
      pch1 = 1
      pch2 = 19
    }else {
      pch1 = 2
      pch2 = 17
    }
    waitinput=FALSE
    for(sim in 1:dim(estlocs)[3]) for(pair in 1:dim(estlocs)[2]) {
      if(waitinput) {cat("Hit ENTER to continue: ");a <- readLines(n=1)}
      points(estlocs[ind,pair,sim,1],estlocs[ind,pair,sim,2],pch=pch1,col=ind,cex=0.5)
    }
    points(pop[detected,][ind,],pch=pch2,col="black",cex=1)
    points(pop[detected,][ind,],pch=pch2,col="white",cex=0.75)
    points(pop[detected,][ind,],pch=pch2,col=ind,cex=0.5)
  }
}
points(posts$x,posts$y,pch=19,col="white",cex=1.25)
points(posts$x,posts$y,pch=1,cex=1.25)
points(posts$x,posts$y,pch=19,col="red",cex=0.5)


