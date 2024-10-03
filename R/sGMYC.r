#' @export

# phy: single tree or sample of trees
# ntrees: number of trees to analyze from phy (1 to ntrees) (default:1)
# subsamp: number of tips to subsample from each delimited species in full analysis (default: 2)
# nreps: number of subsampling replicates (default:100)
# ncores: number of cores for parallel analysis (when analyzing multiple trees only) (default: 1)
# sGMYC(phy = para, ntrees = 1, ncores = 1, nreps=200, subsamp=2)

sGMYC <- function(phy, ntrees=1, subsamp=2, nreps=100, ncores=1, probmat=FALSE) {

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
  results <- list()
  results<-foreach (j = 1:ntrees, .combine=rbind) %dopar% {
  
# from phylo to multyphylo
if (class(phy)=="phylo") {
  phy<-list(phy)
  class(phy)<-"multiPhylo"
}

totalsamp<-0

# process separately if subsampling either one or more tips

if (subsamp==1) {

# analyze each tree

# bootcomp<-c()

# for (j in 1:ntrees) {

 mytree<-phy[[j]]

 # force ultrametric tree
 mytree<-phytools::force.ultrametric(mytree, method="extend")

 #do gmyc
 fullgmyc<-splits::gmyc(mytree, quiet=TRUE)
 
 tryCatch(
 expr = {
 assignments<-splits::spec.list(fullgmyc)

 # sampling one tip per delimited entity
 # bootgmyc<-c()

	toretain<-c()
	for (k in 1:max(assignments[,1])) {
	 toretain<-c(toretain, as.numeric(sample(assignments[assignments[,1]==k,2],1)))
	}
	
	#remove branches
	todrop<-setdiff(c(1:length(mytree$tip.label)), toretain)
	bootree<-ape::drop.tip(mytree, as.character(todrop))

	#do gmyc
	
	tryCatch(
     expr = {newassign<-splits::spec.list(splits::gmyc(bootree, quiet=TRUE))
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
 mytree<-phytools::force.ultrametric(mytree, method="extend")

 #do gmyc
 fullgmyc<-splits::gmyc(mytree, quiet=TRUE)
 
 tryCatch(
 expr = {
 assignments<-splits::spec.list(fullgmyc)

 bootgmyc<-c()

 for (i in 1:nreps) {

	toretain<-c()
	for (k in 1:max(assignments[,1])) {
		if (dim(assignments[assignments[,1]==k,])[1] > subsamp) {
			toretain<-c(toretain, as.numeric(sample(assignments[assignments[,1]==k,2], subsamp, replace=FALSE)))
			} else {
			toretain<-c(toretain, as.numeric(assignments[assignments[,1]==k,][,2]))
		}	
	}
	
	#remove branches
	todrop<-setdiff(c(1:length(mytree$tip.label)), toretain)
	bootree<-ape::drop.tip(mytree, as.character(todrop))

	#do gmyc
	tryCatch(
     expr = {newassign<-splits::spec.list(splits::gmyc(bootree, quiet=TRUE))
#             bootgmyc<-c(bootgmyc, max(newassign[,1]))
             results[[j]]$par<-rbind(results[[j]]$par, c(i,subsamp,max(newassign[,1])))},
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
 mytree<-phytools::force.ultrametric(mytree, method="extend")

 #do gmyc
 fullgmyc<-splits::gmyc(mytree, quiet=TRUE)
 
 tryCatch(
 expr = {
 assignments<-splits::spec.list(fullgmyc)

 #sampling one tip per delimited entity
 bootgmyc<-c()

	toretain<-c()
	for (k in 1:max(assignments[,1])) {
	 toretain<-c(toretain, as.numeric(sample(assignments[assignments[,1]==k,2],1)))
	 k=k+1
	}
	
	#remove branches
	todrop<-setdiff(c(1:length(mytree$tip.label)), toretain)
	bootree<-ape::drop.tip(mytree, as.character(todrop))

	#do gmyc
	
	tryCatch(
     expr = {newassign<-splits::spec.list(splits::gmyc(bootree, quiet=TRUE))
             bootgmyc<-c(bootgmyc, max(newassign[,1]))
             bootcomp<-rbind(bootcomp, c(max(assignments[,1]), min(bootgmyc), max(bootgmyc), mean(bootgmyc), sd(bootgmyc)))},
	 error = function(e) {print("There was an error message.")})

  },
  error = function(e) {print("There was an error message.")})

print(j)

}

write.table(bootcomp, file="bootcomp_1tip")

} else { #sampling more than one tip per delimited entity

#analyze each tree

results<-vector("list", length = ntrees)
par<-c()
tree<-c()
bootcomp<-c()
boottable<-list()

for (j in 1:ntrees) {
 
 mytree<-phy[[j]]

 # force ultrametric tree
  if (ape::is.ultrametric(mytree)==FALSE) {
    mytree<-phytools::force.ultrametric(mytree, method="extend")}

 #do gmyc
 fullgmyc<-splits::gmyc(mytree, quiet=TRUE)
 
 tryCatch(
 expr = {
 assignments<-splits::spec.list(fullgmyc)

 bootree<-list()
 mrca<-list()
 newassign<-list()
 nsp<-c()
 reassign<-list()
 #boottable<-list()

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
	bootree[[i]]<-ape::drop.tip(mytree, todrop)

	#do gmyc
	tryCatch(
     expr = {newgmyc<-splits::gmyc(bootree[[i]], quiet=TRUE)
     		 mrca[[i]]<-newgmyc$MRCA
     		 newassign[[i]]<-splits::spec.list(newgmyc)
             nsp<-c(nsp, max(newassign[[i]][,1]))
             },
	 error = function(e) {print("There was an error message.")})

#	if (max(assignments[,1])>=max(newassign[[i]][,1])) {

		mysamples<-c()
		reassigntable<-c()

		for (l in 1:length(newassign[[i]][,1])) {
 			newsample<-newassign[[i]][l,2]
 			newsp<-newassign[[i]][l,1]
 			if (length(intersect(mysamples,newsample))==0) {
 				oldsp<-assignments[assignments[,2]==newsample,1]
 				mysamples<-assignments[assignments[,1]==oldsp,2]
 				reassigntable<-rbind(reassigntable,cbind(rep(newsp,length(mysamples)),mysamples))
 				} else {
 				reassigntable[reassigntable[,2]==newsample,1]<-newsp
 				}
 		}
		reassign[[i]]<-as.data.frame(reassigntable)
#	}
 }

bootcomp<-rbind(bootcomp, c(max(assignments[,1]), min(nsp), max(nsp), mean(nsp), sd(nsp), nreps, subsamp))

boottable[[j]]<-as.data.frame(table(nsp))
colnames(boottable[[j]])<-c("#species","#reps")

results[[j]]<-list("par"=cbind(c(1:nreps),rep(max(assignments[,1]),nreps),nsp),
	"tree"=mytree,"mrca"=mrca,"assignment"=newassign,"reassignment"=reassign)

cat(paste("Subsampling of tree ", j, ", Done!", "\n", sep=""))
},
error = function(e) {print("There was an error message.")})

}

cat(paste("\n","Number of species in full analysis vs. after subsampling","\n"))
colnames(bootcomp)<-c("full_GMYC","min_sGMYC","max_sGMYC","mean_sGMYC","sd_sGMYC","#reps", "#subsamp per sp")
rownames(bootcomp)<-c("Tree 1")
print(bootcomp)
cat(paste("\n","Number of species vs. number of subsampling replicates","\n"))
names(boottable)<-c("Tree 1")
print(boottable)

class(results)<-"singlebgmyc"
return(results)

} # close else
} # close run with single core
} # end of function
