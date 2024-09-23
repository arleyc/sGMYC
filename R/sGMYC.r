# phy: single tree or sample of trees
# ntrees: number of trees to analyze from phy (1 to ntrees) (default:1)
# subsamp: number of tips to subsample from each delimited species in full analysis (default: 2)
# nreps: number of subsampling replicates (default:100)
# ncores: number of cores for parallel analysis (when analyzing multiple trees only) (default: 1)
# sGMYC(phy = para, ntrees = 1, ncores = 1, nreps=200, subsamp=2)

sGMYC <- function(phy, ntrees=1, subsamp=2, nreps=100, ncores=1) {

# check dependencies

requireNamespace("ape", quietly = TRUE)
requireNamespace("phytools", quietly = TRUE)
requireNamespace("splits", quietly = TRUE)

# check dependency for parallel analysis

if (ncores > 1) {
  requireNamespace("parallel", quietly = TRUE)
  requireNamespace("doParallel", quietly = TRUE)
  mycores<-parallel::detectCores()-1
  if (ncores < mycores) {
    mycores<-ncores}
  myCluster<-parallel::makeCluster(mycores)
  doParallel::registerDoParallel(myCluster)
  results <- list("fullGMYC","min_bootGMYC","max_bootGMYC","mean_bootGMYC")
  results<-foreach (j = 1:ntrees, .combine=rbind) %dopar% {

# from phylo to multyphylo
if (class(phy)=="phylo") {
  phy<-list(phy)
  class(phy)<-"multiPhylo"
}

# process separately if subsampling either one or more tips

if (subsamp==1) {

# analyze each tree

# bootcomp<-c()

# for (j in 1:ntrees) {

 mytree<-phy[[j]]

 # force ultrametric tree
 if (ape::is.ultrametric(mytree)==FALSE) {
    mytree<-phytools::force.ultrametric(mytree, method="extend")}

 #do gmyc
 fullgmyc<-splits::gmyc(mytree)

 tryCatch(
 expr = {
 assignments<-splits::spec.list(fullgmyc)

 # sampling one tip per delimited entity
 # bootgmyc<-c()

	toretain<-c()
	for (k in 1:max(assignments[,1])) {
	 toretain<-c(toretain, sample(assignments[assignments[,1]==k,2],1))
	 k=k+1
	}

	#remove branches
	todrop<-setdiff(mytree$tip.label, toretain)
	bootree<-ape::drop.tip(mytree, todrop)

	#do gmyc

	tryCatch(
     expr = {newassign<-splits::spec.list(splits::gmyc(bootree))
             bootgmyc<-c(bootgmyc, max(newassign[,1]))
             #bootcomp<-rbind(bootcomp, c(max(assignments[,1]), min(bootgmyc), max(bootgmyc), mean(bootgmyc), sd(bootgmyc)))
             results[[j]]<-c(max(assignments[,1]), min(bootgmyc), max(bootgmyc), mean(bootgmyc))},
	 		 error = function(e) {print("There was an error message.")})

  },
  error = function(e) {print("There was an error message.")})

print(j)

# } close for loop
# write.table(bootcomp, file="bootcomp_1tip")

print(results[[j]])

} else { #sampling more than one tip per delimited entity

# analyze each tree

# bootcomp<-c()

# for (j in 1:ntrees) {

 mytree<-phy[[j]]

 # force ultrametric tree
  if (is.ultrametric(mytree)==FALSE) {
     mytree<-phytools::force.ultrametric(mytree, method="extend")}

 #do gmyc
 fullgmyc<-splits::gmyc(mytree)

 tryCatch(
 expr = {
 assignments<-splits::spec.list(fullgmyc)

 bootgmyc<-c()

 for (i in 1:nreps) {

	toretain<-c()
	for (k in 1:max(assignments[,1])) {
		if (dim(assignments[assignments[,1]==k,])[1] > subsamp) {
			toretain<-c(toretain, sample(assignments[assignments[,1]==k,2], subsamp, replace=FALSE))
			} else {
			toretain<-c(toretain, assignments[assignments[,1]==k,][,2])
		}
	}

	#remove branches
	todrop<-setdiff(mytree$tip.label, toretain)
	bootree<-ape::drop.tip(mytree, todrop)

	#do gmyc
	tryCatch(
     expr = {newassign<-splits::spec.list(splits::gmyc(bootree))
             bootgmyc<-c(bootgmyc, max(newassign[,1]))
             results[[j]]<-c(max(assignments[,1]), min(bootgmyc), max(bootgmyc), mean(bootgmyc), sd(bootgmyc))},
	 error = function(e) {print("There was an error message.")})
 }

#bootcomp<-rbind(bootcomp, c(max(assignments[,1]), min(bootgmyc), max(bootgmyc), mean(bootgmyc), sd(bootgmyc)))


},
error = function(e) {print("There was an error message.")})

print(j)
print(results[[j]])

#}

# write.table(bootcomp, file=paste("results_",subsamp,"tips",sep=""))

} # close else
} # close dopar

stopCluster(myCluster)
print(results)
write.table(as.data.frame(results), file=paste("results_",subsamp,"tips_",mycores, "core(s)", sep=""))

} # close if parallel

else { # run with single core

# from phylo to multyphylo
if (class(phy)=="phylo") {
  phy<-list(phy)
  class(phy)<-"multiPhylo"
}

# process separately if subsampling either one or more tips

if (subsamp==1) {

#analyze each tree

bootcomp<-c()

for (j in 1:ntrees) {

 mytree<-phy[[j]]

 # force ultrametric tree
  if (is.ultrametric(mytree)==FALSE) {
     mytree<-phytools::force.ultrametric(mytree, method="extend")}

 #do gmyc
 fullgmyc<-splits::gmyc(mytree)

 tryCatch(
 expr = {
 assignments<-splits::spec.list(fullgmyc)

 #sampling one tip per delimited entity
 bootgmyc<-c()

	toretain<-c()
	for (k in 1:max(assignments[,1])) {
	 toretain<-c(toretain, sample(assignments[assignments[,1]==k,2],1))
	 k=k+1
	}

	#remove branches
	todrop<-setdiff(mytree$tip.label, toretain)
	bootree<-ape::drop.tip(mytree, todrop)

	#do gmyc

	tryCatch(
     expr = {newassign<-splits::spec.list(splits::gmyc(bootree))
             bootgmyc<-c(bootgmyc, max(newassign[,1]))
             bootcomp<-rbind(bootcomp, c(max(assignments[,1]), min(bootgmyc), max(bootgmyc), mean(bootgmyc), sd(bootgmyc), nreps, subsamp))},
	 error = function(e) {print("There was an error message.")})

  },
  error = function(e) {print("There was an error message.")})

print(j)

}

colnames(bootcomp)<-c("fullGMYC","min_bootGMYC","max_bootGMYC","mean_bootGMYC","sd_bootGMYC","reps","subsamples/sp")
#write.table(bootcomp, file="bootcomp_1tip")
boottable<-table(bootgmyc)
return(list(boottable, bootcomp))

} else { #sampling more than one tip per delimited entity

#analyze each tree

bootcomp<-c()

for (j in 1:ntrees) {

 mytree<-phy[[j]]

 # force ultrametric tree
   if (is.ultrametric(mytree)==FALSE) {
     mytree<-phytools::force.ultrametric(mytree, method="extend")}

 #do gmyc
 fullgmyc<-splits::gmyc(mytree)

 tryCatch(
 expr = {
 assignments<-splits::spec.list(fullgmyc)

 bootgmyc<-c()
 boottable<-list()

 for (i in 1:nreps) {

	toretain<-c()
	for (k in 1:max(assignments[,1])) {
		if (dim(assignments[assignments[,1]==k,])[1] > subsamp) {
			toretain<-c(toretain, sample(assignments[assignments[,1]==k,2], subsamp, replace=FALSE))
			} else {
			toretain<-c(toretain, assignments[assignments[,1]==k,][,2])
		}
	}

	#remove branches
	todrop<-setdiff(mytree$tip.label, toretain)
	bootree<-ape::drop.tip(mytree, todrop)

	#do gmyc
	tryCatch(
     expr = {newassign<-splits::spec.list(splits::gmyc(bootree))
             bootgmyc<-c(bootgmyc, max(newassign[,1]))},
	 error = function(e) {print("There was an error message.")})
 }

bootcomp<-rbind(bootcomp, c(max(assignments[,1]), min(bootgmyc), max(bootgmyc), mean(bootgmyc), sd(bootgmyc), nreps, subsamp))
boottable[[j]]<-table(bootgmyc)

},
error = function(e) {print("There was an error message.")})

print(j)

}

colnames(bootcomp)<-c("fullGMYC","min_bootGMYC","max_bootGMYC","mean_bootGMYC","sd_bootGMYC","reps", "subsamples/sp")
#write.table(bootcomp, file=paste("bootcomp_",subsamp,"tips",sep=""))
names(boottable) <- c("Number of species (first row) and number of replicates (second row)")
return(list(boottable, bootcomp))

} # close else
} # close run with single core
} # end of function
