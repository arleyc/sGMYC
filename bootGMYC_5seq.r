setwd("/Users/arley/Desktop/Arley/Projects/bootGMYC/simulation_results/5seq_sp/16spp_symmetric")

#simulated tree
simtrees<-ape::read.nexus(file="simtrees_theta5_tau2_5seq.nex")

#analyze each simulated tree

bootcomp<-c()

for (j in 1:1000) {

 simtree<-simtrees[[j]]

 # force ultrametric tree
 simtree<-phytools::force.ultrametric(simtree, method="extend")

 #do gmyc
 fullgmyc<-splits::gmyc(simtree)
 
 tryCatch(
 expr = {
 assignments<-splits::spec.list(fullgmyc)

 #sampling one tip per delimited entity
 bootgmyc<-c()

# for (i in 1:1) {

	toretain<-c()
	for (k in 1:max(assignments[,1])) {
	 toretain<-c(toretain, as.numeric(sample(assignments[assignments[,1]==k,2],1)))
	 k=k+1
	}
	
	#remove branches
	todrop<-setdiff(c(1:length(simtree$tip.label)), toretain)
	bootree<-ape::drop.tip(simtree, as.character(todrop))

	#do gmyc
	
	tryCatch(
     expr = {newassign<-splits::spec.list(splits::gmyc(bootree))
             bootgmyc<-c(bootgmyc, max(newassign[,1]))
             bootcomp<-rbind(bootcomp, c(max(assignments[,1]), min(bootgmyc), max(bootgmyc), mean(bootgmyc), sd(bootgmyc)))},
	 error = function(e) {print("There was an error message.")})

  },
  error = function(e) {print("There was an error message.")})

#	print(i)
#	i<-i+1
	
#
#compare number of delimited entities

#plot(hist(bootgmyc), xlim=c(8,22), main="8 spp (green), 5 seq/sp, 1 tip/entity (gray), 
#	full gymc with 40 tips, 22 entities (red)")
#arrows(x0=8, y0=5, y1=0, col="green")
#arrows(x0=max(assignments[,1]), y0=5, y1=0, col="red")

#save max, min and mode of bootstrapped tree

print(j)
j=j+1

}

write.table(bootcomp, file="bootcomp_1tip_theta5_tau2_5seq")

#sampling 2 tips per delimited entity

#simulated tree
#simtrees<-ape::read.nexus(file="simtrees.nex")

#analyze each simulated tree

bootcomp<-c()

for (j in 1:100) {

 simtree<-simtrees[[j]]

 # force ultrametric tree
 simtree<-phytools::force.ultrametric(simtree, method="extend")

 #do gmyc
 fullgmyc<-splits::gmyc(simtree)
 
 tryCatch(
 expr = {

 assignments<-splits::spec.list(fullgmyc)

 bootgmyc<-c()

 for (i in 1:100) {

	toretain<-c()
	for (k in 1:max(assignments[,1])) {
		if (dim(assignments[assignments[,1]==k,])[1] > 2) {
			toretain<-c(toretain, as.numeric(sample(assignments[assignments[,1]==k,2],2, replace=FALSE)))
			} else {
			toretain<-c(toretain, as.numeric(assignments[assignments[,1]==k,][,2]))
		}	
		k=k+1
	}
	
	#remove branches
	todrop<-setdiff(c(1:length(simtree$tip.label)), toretain)
	bootree<-ape::drop.tip(simtree, as.character(todrop))

	#do gmyc
	tryCatch(
     expr = {newassign<-splits::spec.list(splits::gmyc(bootree))
             bootgmyc<-c(bootgmyc, max(newassign[,1]))},
	 error = function(e) {print("There was an error message.")})

 i=i+1
 }

#compare number of delimited entities

#plot(hist(bootgmyc), xlim=c(0,25), main="8 spp (green), 5 seq/sp, 2 tips/entity (gray), 
#	full gymc with 40 tips, 22 entities (red)")
#arrows(x0=8, y0=20,y1=0, col="green")
#arrows(x0=max(assignments[,1]), y0=20,y1=0, col="red")

bootcomp<-rbind(bootcomp, c(max(assignments[,1]), min(bootgmyc), max(bootgmyc), mean(bootgmyc), sd(bootgmyc)))

},
error = function(e) {print("There was an error message.")})

print(j)
j=j+1

}

write.table(bootcomp, file="bootcomp_2tip_theta5_tau2_5seq")

#sampling 3 tips per delimited entity
#analyze each simulated tree

bootcomp<-c()

for (j in 1:100) {

 simtree<-simtrees[[j]]

 # force ultrametric tree
 simtree<-phytools::force.ultrametric(simtree, method="extend")

 #do gmyc
 fullgmyc<-splits::gmyc(simtree)

 tryCatch(
 expr = {

 assignments<-splits::spec.list(fullgmyc)

 bootgmyc<-c()

 for (i in 1:100) {

	toretain<-c()
	for (k in 1:max(assignments[,1])) {
		if (dim(assignments[assignments[,1]==k,])[1] > 3) {
			toretain<-c(toretain, as.numeric(sample(assignments[assignments[,1]==k,2],3, replace=FALSE)))
			} else {
			toretain<-c(toretain, as.numeric(assignments[assignments[,1]==k,][,2]))
		}	
		k=k+1
	}
	
	#remove branches
	todrop<-setdiff(c(1:length(simtree$tip.label)), toretain)
	bootree<-ape::drop.tip(simtree, todrop)

	#do gmyc
	tryCatch(
     expr = {newassign<-splits::spec.list(splits::gmyc(bootree))
             bootgmyc<-c(bootgmyc, max(newassign[,1]))},
	 error = function(e) {print("There was an error message.")})

 i=i+1
 }

#compare number of delimited entities

#plot(hist(bootgmyc), xlim=c(0,25), main="8 spp (green), 5 seq/sp, 2 tips/entity (gray), 
#	full gymc with 40 tips, 22 entities (red)")
#arrows(x0=8, y0=20,y1=0, col="green")
#arrows(x0=max(assignments[,1]), y0=20,y1=0, col="red")

bootcomp<-rbind(bootcomp, c(max(assignments[,1]), min(bootgmyc), max(bootgmyc), mean(bootgmyc), sd(bootgmyc)))

},
error = function(e) {print("There was an error message.")})

print(j)
j=j+1

}

write.table(bootcomp, file="bootcomp_3tip_theta5_tau2_5seq")

#sampling 4 tips per delimited entity
#analyze each simulated tree

bootcomp<-c()

for (j in 1:100) {

 simtree<-simtrees[[j]]

 # force ultrametric tree
 simtree<-phytools::force.ultrametric(simtree, method="extend")

 #do gmyc
 fullgmyc<-splits::gmyc(simtree)

 tryCatch(
 expr = {

 assignments<-splits::spec.list(fullgmyc)

 bootgmyc<-c()

 for (i in 1:100) {

	toretain<-c()
	for (k in 1:max(assignments[,1])) {
		if (dim(assignments[assignments[,1]==k,])[1] > 4) {
			toretain<-c(toretain, as.numeric(sample(assignments[assignments[,1]==k,2],4, replace=FALSE)))
			} else {
			toretain<-c(toretain, as.numeric(assignments[assignments[,1]==k,][,2]))
		}	
		k=k+1
	}
	
	#remove branches
	todrop<-setdiff(c(1:length(simtree$tip.label)), toretain)
	bootree<-ape::drop.tip(simtree, todrop)

	#do gmyc
	tryCatch(
     expr = {newassign<-splits::spec.list(splits::gmyc(bootree))
             bootgmyc<-c(bootgmyc, max(newassign[,1]))},
	 error = function(e) {print("There was an error message.")})

 i=i+1
 }


#compare number of delimited entities

#plot(hist(bootgmyc), xlim=c(0,25), main="8 spp (green), 5 seq/sp, 2 tips/entity (gray), 
#	full gymc with 40 tips, 22 entities (red)")
#arrows(x0=8, y0=20,y1=0, col="green")
#arrows(x0=max(assignments[,1]), y0=20,y1=0, col="red")

bootcomp<-rbind(bootcomp, c(max(assignments[,1]), min(bootgmyc), max(bootgmyc), mean(bootgmyc), sd(bootgmyc)))

},
error = function(e) {print("There was an error message.")})

print(j)
j=j+1

}

write.table(bootcomp, file="bootcomp_4tip_theta5_tau2_5seq")


############################################################

#simulated tree
simtrees<-ape::read.nexus(file="simtrees_theta5_tau02_5seq.nex")

#analyze each simulated tree

bootcomp<-c()

for (j in 1:1000) {

 simtree<-simtrees[[j]]

 # force ultrametric tree
 simtree<-phytools::force.ultrametric(simtree, method="extend")

 #do gmyc
 fullgmyc<-splits::gmyc(simtree)
 
 tryCatch(
 expr = {
 assignments<-splits::spec.list(fullgmyc)

 #sampling one tip per delimited entity
 bootgmyc<-c()

# for (i in 1:1) {

	toretain<-c()
	for (k in 1:max(assignments[,1])) {
	 toretain<-c(toretain, as.numeric(sample(assignments[assignments[,1]==k,2],1)))
	 k=k+1
	}
	
	#remove branches
	todrop<-setdiff(c(1:length(simtree$tip.label)), toretain)
	bootree<-ape::drop.tip(simtree, as.character(todrop))

	#do gmyc
	
	tryCatch(
     expr = {newassign<-splits::spec.list(splits::gmyc(bootree))
             bootgmyc<-c(bootgmyc, max(newassign[,1]))
             bootcomp<-rbind(bootcomp, c(max(assignments[,1]), min(bootgmyc), max(bootgmyc), mean(bootgmyc), sd(bootgmyc)))},
	 error = function(e) {print("There was an error message.")})

  },
  error = function(e) {print("There was an error message.")})

#	print(i)
#	i<-i+1
	
#
#compare number of delimited entities

#plot(hist(bootgmyc), xlim=c(8,22), main="8 spp (green), 5 seq/sp, 1 tip/entity (gray), 
#	full gymc with 40 tips, 22 entities (red)")
#arrows(x0=8, y0=5, y1=0, col="green")
#arrows(x0=max(assignments[,1]), y0=5, y1=0, col="red")

#save max, min and mode of bootstrapped tree

print(j)
j=j+1

}

write.table(bootcomp, file="bootcomp_1tip_theta5_tau02_5seq")

#sampling 2 tips per delimited entity

#simulated tree
#simtrees<-ape::read.nexus(file="simtrees.nex")

#analyze each simulated tree

bootcomp<-c()

for (j in 1:100) {

 simtree<-simtrees[[j]]

 # force ultrametric tree
 simtree<-phytools::force.ultrametric(simtree, method="extend")

 #do gmyc
 fullgmyc<-splits::gmyc(simtree)
 
 tryCatch(
 expr = {

 assignments<-splits::spec.list(fullgmyc)

 bootgmyc<-c()

 for (i in 1:100) {

	toretain<-c()
	for (k in 1:max(assignments[,1])) {
		if (dim(assignments[assignments[,1]==k,])[1] > 2) {
			toretain<-c(toretain, as.numeric(sample(assignments[assignments[,1]==k,2],2, replace=FALSE)))
			} else {
			toretain<-c(toretain, as.numeric(assignments[assignments[,1]==k,][,2]))
		}	
		k=k+1
	}
	
	#remove branches
	todrop<-setdiff(c(1:length(simtree$tip.label)), toretain)
	bootree<-ape::drop.tip(simtree, as.character(todrop))

	#do gmyc
	tryCatch(
     expr = {newassign<-splits::spec.list(splits::gmyc(bootree))
             bootgmyc<-c(bootgmyc, max(newassign[,1]))},
	 error = function(e) {print("There was an error message.")})

 i=i+1
 }

#compare number of delimited entities

#plot(hist(bootgmyc), xlim=c(0,25), main="8 spp (green), 5 seq/sp, 2 tips/entity (gray), 
#	full gymc with 40 tips, 22 entities (red)")
#arrows(x0=8, y0=20,y1=0, col="green")
#arrows(x0=max(assignments[,1]), y0=20,y1=0, col="red")

bootcomp<-rbind(bootcomp, c(max(assignments[,1]), min(bootgmyc), max(bootgmyc), mean(bootgmyc), sd(bootgmyc)))

},
error = function(e) {print("There was an error message.")})

print(j)
j=j+1

}

write.table(bootcomp, file="bootcomp_2tip_theta5_tau02_5seq")

#sampling 3 tips per delimited entity
#analyze each simulated tree

bootcomp<-c()

for (j in 1:100) {

 simtree<-simtrees[[j]]

 # force ultrametric tree
 simtree<-phytools::force.ultrametric(simtree, method="extend")

 #do gmyc
 fullgmyc<-splits::gmyc(simtree)

 tryCatch(
 expr = {

 assignments<-splits::spec.list(fullgmyc)

 bootgmyc<-c()

 for (i in 1:100) {

	toretain<-c()
	for (k in 1:max(assignments[,1])) {
		if (dim(assignments[assignments[,1]==k,])[1] > 3) {
			toretain<-c(toretain, as.numeric(sample(assignments[assignments[,1]==k,2],3, replace=FALSE)))
			} else {
			toretain<-c(toretain, as.numeric(assignments[assignments[,1]==k,][,2]))
		}	
		k=k+1
	}
	
	#remove branches
	todrop<-setdiff(c(1:length(simtree$tip.label)), toretain)
	bootree<-ape::drop.tip(simtree, todrop)

	#do gmyc
	tryCatch(
     expr = {newassign<-splits::spec.list(splits::gmyc(bootree))
             bootgmyc<-c(bootgmyc, max(newassign[,1]))},
	 error = function(e) {print("There was an error message.")})

 i=i+1
 }

#compare number of delimited entities

#plot(hist(bootgmyc), xlim=c(0,25), main="8 spp (green), 5 seq/sp, 2 tips/entity (gray), 
#	full gymc with 40 tips, 22 entities (red)")
#arrows(x0=8, y0=20,y1=0, col="green")
#arrows(x0=max(assignments[,1]), y0=20,y1=0, col="red")

bootcomp<-rbind(bootcomp, c(max(assignments[,1]), min(bootgmyc), max(bootgmyc), mean(bootgmyc), sd(bootgmyc)))

},
error = function(e) {print("There was an error message.")})

print(j)
j=j+1

}

write.table(bootcomp, file="bootcomp_3tip_theta5_tau02_5seq")

#sampling 4 tips per delimited entity
#analyze each simulated tree

bootcomp<-c()

for (j in 1:100) {

 simtree<-simtrees[[j]]

 # force ultrametric tree
 simtree<-phytools::force.ultrametric(simtree, method="extend")

 #do gmyc
 fullgmyc<-splits::gmyc(simtree)

 tryCatch(
 expr = {

 assignments<-splits::spec.list(fullgmyc)

 bootgmyc<-c()

 for (i in 1:100) {

	toretain<-c()
	for (k in 1:max(assignments[,1])) {
		if (dim(assignments[assignments[,1]==k,])[1] > 4) {
			toretain<-c(toretain, as.numeric(sample(assignments[assignments[,1]==k,2],4, replace=FALSE)))
			} else {
			toretain<-c(toretain, as.numeric(assignments[assignments[,1]==k,][,2]))
		}	
		k=k+1
	}
	
	#remove branches
	todrop<-setdiff(c(1:length(simtree$tip.label)), toretain)
	bootree<-ape::drop.tip(simtree, todrop)

	#do gmyc
	tryCatch(
     expr = {newassign<-splits::spec.list(splits::gmyc(bootree))
             bootgmyc<-c(bootgmyc, max(newassign[,1]))},
	 error = function(e) {print("There was an error message.")})

 i=i+1
 }


#compare number of delimited entities

#plot(hist(bootgmyc), xlim=c(0,25), main="8 spp (green), 5 seq/sp, 2 tips/entity (gray), 
#	full gymc with 40 tips, 22 entities (red)")
#arrows(x0=8, y0=20,y1=0, col="green")
#arrows(x0=max(assignments[,1]), y0=20,y1=0, col="red")

bootcomp<-rbind(bootcomp, c(max(assignments[,1]), min(bootgmyc), max(bootgmyc), mean(bootgmyc), sd(bootgmyc)))

},
error = function(e) {print("There was an error message.")})

print(j)
j=j+1

}

write.table(bootcomp, file="bootcomp_4tip_theta5_tau02_5seq")

############################################################

#simulated tree
simtrees<-ape::read.nexus(file="simtrees_theta50_tau2_5seq.nex")

#analyze each simulated tree

bootcomp<-c()

for (j in 1:1000) {

 simtree<-simtrees[[j]]

 # force ultrametric tree
 simtree<-phytools::force.ultrametric(simtree, method="extend")

 #do gmyc
 fullgmyc<-splits::gmyc(simtree)
 
 tryCatch(
 expr = {
 assignments<-splits::spec.list(fullgmyc)

 #sampling one tip per delimited entity
 bootgmyc<-c()

# for (i in 1:1) {

	toretain<-c()
	for (k in 1:max(assignments[,1])) {
	 toretain<-c(toretain, as.numeric(sample(assignments[assignments[,1]==k,2],1)))
	 k=k+1
	}
	
	#remove branches
	todrop<-setdiff(c(1:length(simtree$tip.label)), toretain)
	bootree<-ape::drop.tip(simtree, as.character(todrop))

	#do gmyc
	
	tryCatch(
     expr = {newassign<-splits::spec.list(splits::gmyc(bootree))
             bootgmyc<-c(bootgmyc, max(newassign[,1]))
             bootcomp<-rbind(bootcomp, c(max(assignments[,1]), min(bootgmyc), max(bootgmyc), mean(bootgmyc), sd(bootgmyc)))},
	 error = function(e) {print("There was an error message.")})

  },
  error = function(e) {print("There was an error message.")})

#	print(i)
#	i<-i+1
	
#
#compare number of delimited entities

#plot(hist(bootgmyc), xlim=c(8,22), main="8 spp (green), 5 seq/sp, 1 tip/entity (gray), 
#	full gymc with 40 tips, 22 entities (red)")
#arrows(x0=8, y0=5, y1=0, col="green")
#arrows(x0=max(assignments[,1]), y0=5, y1=0, col="red")

#save max, min and mode of bootstrapped tree

print(j)
j=j+1

}

write.table(bootcomp, file="bootcomp_1tip_theta50_tau2_5seq")

#sampling 2 tips per delimited entity

#simulated tree
#simtrees<-ape::read.nexus(file="simtrees.nex")

#analyze each simulated tree

bootcomp<-c()

for (j in 1:100) {

 simtree<-simtrees[[j]]

 # force ultrametric tree
 simtree<-phytools::force.ultrametric(simtree, method="extend")

 #do gmyc
 fullgmyc<-splits::gmyc(simtree)
 
 tryCatch(
 expr = {

 assignments<-splits::spec.list(fullgmyc)

 bootgmyc<-c()

 for (i in 1:100) {

	toretain<-c()
	for (k in 1:max(assignments[,1])) {
		if (dim(assignments[assignments[,1]==k,])[1] > 2) {
			toretain<-c(toretain, as.numeric(sample(assignments[assignments[,1]==k,2],2, replace=FALSE)))
			} else {
			toretain<-c(toretain, as.numeric(assignments[assignments[,1]==k,][,2]))
		}	
		k=k+1
	}
	
	#remove branches
	todrop<-setdiff(c(1:length(simtree$tip.label)), toretain)
	bootree<-ape::drop.tip(simtree, as.character(todrop))

	#do gmyc
	tryCatch(
     expr = {newassign<-splits::spec.list(splits::gmyc(bootree))
             bootgmyc<-c(bootgmyc, max(newassign[,1]))},
	 error = function(e) {print("There was an error message.")})

 i=i+1
 }

#compare number of delimited entities

#plot(hist(bootgmyc), xlim=c(0,25), main="8 spp (green), 5 seq/sp, 2 tips/entity (gray), 
#	full gymc with 40 tips, 22 entities (red)")
#arrows(x0=8, y0=20,y1=0, col="green")
#arrows(x0=max(assignments[,1]), y0=20,y1=0, col="red")

bootcomp<-rbind(bootcomp, c(max(assignments[,1]), min(bootgmyc), max(bootgmyc), mean(bootgmyc), sd(bootgmyc)))

},
error = function(e) {print("There was an error message.")})

print(j)
j=j+1

}

write.table(bootcomp, file="bootcomp_2tip_theta50_tau2_5seq")

#sampling 3 tips per delimited entity
#analyze each simulated tree

bootcomp<-c()

for (j in 1:100) {

 simtree<-simtrees[[j]]

 # force ultrametric tree
 simtree<-phytools::force.ultrametric(simtree, method="extend")

 #do gmyc
 fullgmyc<-splits::gmyc(simtree)

 tryCatch(
 expr = {

 assignments<-splits::spec.list(fullgmyc)

 bootgmyc<-c()

 for (i in 1:100) {

	toretain<-c()
	for (k in 1:max(assignments[,1])) {
		if (dim(assignments[assignments[,1]==k,])[1] > 3) {
			toretain<-c(toretain, as.numeric(sample(assignments[assignments[,1]==k,2],3, replace=FALSE)))
			} else {
			toretain<-c(toretain, as.numeric(assignments[assignments[,1]==k,][,2]))
		}	
		k=k+1
	}
	
	#remove branches
	todrop<-setdiff(c(1:length(simtree$tip.label)), toretain)
	bootree<-ape::drop.tip(simtree, todrop)

	#do gmyc
	tryCatch(
     expr = {newassign<-splits::spec.list(splits::gmyc(bootree))
             bootgmyc<-c(bootgmyc, max(newassign[,1]))},
	 error = function(e) {print("There was an error message.")})

 i=i+1
 }

#compare number of delimited entities

#plot(hist(bootgmyc), xlim=c(0,25), main="8 spp (green), 5 seq/sp, 2 tips/entity (gray), 
#	full gymc with 40 tips, 22 entities (red)")
#arrows(x0=8, y0=20,y1=0, col="green")
#arrows(x0=max(assignments[,1]), y0=20,y1=0, col="red")

bootcomp<-rbind(bootcomp, c(max(assignments[,1]), min(bootgmyc), max(bootgmyc), mean(bootgmyc), sd(bootgmyc)))

},
error = function(e) {print("There was an error message.")})

print(j)
j=j+1

}

write.table(bootcomp, file="bootcomp_3tip_theta50_tau2_5seq")

#sampling 4 tips per delimited entity
#analyze each simulated tree

bootcomp<-c()

for (j in 1:100) {

 simtree<-simtrees[[j]]

 # force ultrametric tree
 simtree<-phytools::force.ultrametric(simtree, method="extend")

 #do gmyc
 fullgmyc<-splits::gmyc(simtree)

 tryCatch(
 expr = {

 assignments<-splits::spec.list(fullgmyc)

 bootgmyc<-c()

 for (i in 1:100) {

	toretain<-c()
	for (k in 1:max(assignments[,1])) {
		if (dim(assignments[assignments[,1]==k,])[1] > 4) {
			toretain<-c(toretain, as.numeric(sample(assignments[assignments[,1]==k,2],4, replace=FALSE)))
			} else {
			toretain<-c(toretain, as.numeric(assignments[assignments[,1]==k,][,2]))
		}	
		k=k+1
	}
	
	#remove branches
	todrop<-setdiff(c(1:length(simtree$tip.label)), toretain)
	bootree<-ape::drop.tip(simtree, todrop)

	#do gmyc
	tryCatch(
     expr = {newassign<-splits::spec.list(splits::gmyc(bootree))
             bootgmyc<-c(bootgmyc, max(newassign[,1]))},
	 error = function(e) {print("There was an error message.")})

 i=i+1
 }


#compare number of delimited entities

#plot(hist(bootgmyc), xlim=c(0,25), main="8 spp (green), 5 seq/sp, 2 tips/entity (gray), 
#	full gymc with 40 tips, 22 entities (red)")
#arrows(x0=8, y0=20,y1=0, col="green")
#arrows(x0=max(assignments[,1]), y0=20,y1=0, col="red")

bootcomp<-rbind(bootcomp, c(max(assignments[,1]), min(bootgmyc), max(bootgmyc), mean(bootgmyc), sd(bootgmyc)))

},
error = function(e) {print("There was an error message.")})

print(j)
j=j+1

}

write.table(bootcomp, file="bootcomp_4tip_theta50_tau2_5seq")

############################################################

#simulated tree
simtrees<-ape::read.nexus(file="simtrees_theta50_tau02_5seq.nex")

#analyze each simulated tree

bootcomp<-c()

for (j in 1:1000) {

 simtree<-simtrees[[j]]

 # force ultrametric tree
 simtree<-phytools::force.ultrametric(simtree, method="extend")

 #do gmyc
 fullgmyc<-splits::gmyc(simtree)
 
 tryCatch(
 expr = {
 assignments<-splits::spec.list(fullgmyc)

 #sampling one tip per delimited entity
 bootgmyc<-c()

# for (i in 1:1) {

	toretain<-c()
	for (k in 1:max(assignments[,1])) {
	 toretain<-c(toretain, as.numeric(sample(assignments[assignments[,1]==k,2],1)))
	 k=k+1
	}
	
	#remove branches
	todrop<-setdiff(c(1:length(simtree$tip.label)), toretain)
	bootree<-ape::drop.tip(simtree, as.character(todrop))

	#do gmyc
	
	tryCatch(
     expr = {newassign<-splits::spec.list(splits::gmyc(bootree))
             bootgmyc<-c(bootgmyc, max(newassign[,1]))
             bootcomp<-rbind(bootcomp, c(max(assignments[,1]), min(bootgmyc), max(bootgmyc), mean(bootgmyc), sd(bootgmyc)))},
	 error = function(e) {print("There was an error message.")})

  },
  error = function(e) {print("There was an error message.")})

#	print(i)
#	i<-i+1
	
#
#compare number of delimited entities

#plot(hist(bootgmyc), xlim=c(8,22), main="8 spp (green), 5 seq/sp, 1 tip/entity (gray), 
#	full gymc with 40 tips, 22 entities (red)")
#arrows(x0=8, y0=5, y1=0, col="green")
#arrows(x0=max(assignments[,1]), y0=5, y1=0, col="red")

#save max, min and mode of bootstrapped tree

print(j)
j=j+1

}

write.table(bootcomp, file="bootcomp_1tip_theta50_tau02_5seq")

#sampling 2 tips per delimited entity

#simulated tree
#simtrees<-ape::read.nexus(file="simtrees.nex")

#analyze each simulated tree

bootcomp<-c()

for (j in 1:100) {

 simtree<-simtrees[[j]]

 # force ultrametric tree
 simtree<-phytools::force.ultrametric(simtree, method="extend")

 #do gmyc
 fullgmyc<-splits::gmyc(simtree)
 
 tryCatch(
 expr = {

 assignments<-splits::spec.list(fullgmyc)

 bootgmyc<-c()

 for (i in 1:100) {

	toretain<-c()
	for (k in 1:max(assignments[,1])) {
		if (dim(assignments[assignments[,1]==k,])[1] > 2) {
			toretain<-c(toretain, as.numeric(sample(assignments[assignments[,1]==k,2],2, replace=FALSE)))
			} else {
			toretain<-c(toretain, as.numeric(assignments[assignments[,1]==k,][,2]))
		}	
		k=k+1
	}
	
	#remove branches
	todrop<-setdiff(c(1:length(simtree$tip.label)), toretain)
	bootree<-ape::drop.tip(simtree, as.character(todrop))

	#do gmyc
	tryCatch(
     expr = {newassign<-splits::spec.list(splits::gmyc(bootree))
             bootgmyc<-c(bootgmyc, max(newassign[,1]))},
	 error = function(e) {print("There was an error message.")})

 i=i+1
 }

#compare number of delimited entities

#plot(hist(bootgmyc), xlim=c(0,25), main="8 spp (green), 5 seq/sp, 2 tips/entity (gray), 
#	full gymc with 40 tips, 22 entities (red)")
#arrows(x0=8, y0=20,y1=0, col="green")
#arrows(x0=max(assignments[,1]), y0=20,y1=0, col="red")

bootcomp<-rbind(bootcomp, c(max(assignments[,1]), min(bootgmyc), max(bootgmyc), mean(bootgmyc), sd(bootgmyc)))

},
error = function(e) {print("There was an error message.")})

print(j)
j=j+1

}

write.table(bootcomp, file="bootcomp_2tip_theta50_tau02_5seq")

#sampling 3 tips per delimited entity
#analyze each simulated tree

bootcomp<-c()

for (j in 1:100) {

 simtree<-simtrees[[j]]

 # force ultrametric tree
 simtree<-phytools::force.ultrametric(simtree, method="extend")

 #do gmyc
 fullgmyc<-splits::gmyc(simtree)

 tryCatch(
 expr = {

 assignments<-splits::spec.list(fullgmyc)

 bootgmyc<-c()

 for (i in 1:100) {

	toretain<-c()
	for (k in 1:max(assignments[,1])) {
		if (dim(assignments[assignments[,1]==k,])[1] > 3) {
			toretain<-c(toretain, as.numeric(sample(assignments[assignments[,1]==k,2],3, replace=FALSE)))
			} else {
			toretain<-c(toretain, as.numeric(assignments[assignments[,1]==k,][,2]))
		}	
		k=k+1
	}
	
	#remove branches
	todrop<-setdiff(c(1:length(simtree$tip.label)), toretain)
	bootree<-ape::drop.tip(simtree, todrop)

	#do gmyc
	tryCatch(
     expr = {newassign<-splits::spec.list(splits::gmyc(bootree))
             bootgmyc<-c(bootgmyc, max(newassign[,1]))},
	 error = function(e) {print("There was an error message.")})

 i=i+1
 }

#compare number of delimited entities

#plot(hist(bootgmyc), xlim=c(0,25), main="8 spp (green), 5 seq/sp, 2 tips/entity (gray), 
#	full gymc with 40 tips, 22 entities (red)")
#arrows(x0=8, y0=20,y1=0, col="green")
#arrows(x0=max(assignments[,1]), y0=20,y1=0, col="red")

bootcomp<-rbind(bootcomp, c(max(assignments[,1]), min(bootgmyc), max(bootgmyc), mean(bootgmyc), sd(bootgmyc)))

},
error = function(e) {print("There was an error message.")})

print(j)
j=j+1

}

write.table(bootcomp, file="bootcomp_3tip_theta50_tau02_5seq")

#sampling 4 tips per delimited entity
#analyze each simulated tree

bootcomp<-c()

for (j in 1:100) {

 simtree<-simtrees[[j]]

 # force ultrametric tree
 simtree<-phytools::force.ultrametric(simtree, method="extend")

 #do gmyc
 fullgmyc<-splits::gmyc(simtree)

 tryCatch(
 expr = {

 assignments<-splits::spec.list(fullgmyc)

 bootgmyc<-c()

 for (i in 1:100) {

	toretain<-c()
	for (k in 1:max(assignments[,1])) {
		if (dim(assignments[assignments[,1]==k,])[1] > 4) {
			toretain<-c(toretain, as.numeric(sample(assignments[assignments[,1]==k,2],4, replace=FALSE)))
			} else {
			toretain<-c(toretain, as.numeric(assignments[assignments[,1]==k,][,2]))
		}	
		k=k+1
	}
	
	#remove branches
	todrop<-setdiff(c(1:length(simtree$tip.label)), toretain)
	bootree<-ape::drop.tip(simtree, todrop)

	#do gmyc
	tryCatch(
     expr = {newassign<-splits::spec.list(splits::gmyc(bootree))
             bootgmyc<-c(bootgmyc, max(newassign[,1]))},
	 error = function(e) {print("There was an error message.")})

 i=i+1
 }


#compare number of delimited entities

#plot(hist(bootgmyc), xlim=c(0,25), main="8 spp (green), 5 seq/sp, 2 tips/entity (gray), 
#	full gymc with 40 tips, 22 entities (red)")
#arrows(x0=8, y0=20,y1=0, col="green")
#arrows(x0=max(assignments[,1]), y0=20,y1=0, col="red")

bootcomp<-rbind(bootcomp, c(max(assignments[,1]), min(bootgmyc), max(bootgmyc), mean(bootgmyc), sd(bootgmyc)))

},
error = function(e) {print("There was an error message.")})

print(j)
j=j+1

}

write.table(bootcomp, file="bootcomp_4tip_theta50_tau02_5seq")


###############################################################
###############################################################

setwd("/Users/arley/Desktop/Arley/Projects/bootGMYC/simulation_results/10seq_sp/16spp_symmetric")

#simulated tree
simtrees<-ape::read.nexus(file="simtrees_theta5_tau2_10seq.nex")

#analyze each simulated tree

bootcomp<-c()

for (j in 1:1000) {

 simtree<-simtrees[[j]]

 # force ultrametric tree
 simtree<-phytools::force.ultrametric(simtree, method="extend")

 #do gmyc
 fullgmyc<-splits::gmyc(simtree)
 
 tryCatch(
 expr = {
 assignments<-splits::spec.list(fullgmyc)

 #sampling one tip per delimited entity
 bootgmyc<-c()

# for (i in 1:1) {

	toretain<-c()
	for (k in 1:max(assignments[,1])) {
	 toretain<-c(toretain, as.numeric(sample(assignments[assignments[,1]==k,2],1)))
	 k=k+1
	}
	
	#remove branches
	todrop<-setdiff(c(1:length(simtree$tip.label)), toretain)
	bootree<-ape::drop.tip(simtree, as.character(todrop))

	#do gmyc
	
	tryCatch(
     expr = {newassign<-splits::spec.list(splits::gmyc(bootree))
             bootgmyc<-c(bootgmyc, max(newassign[,1]))
             bootcomp<-rbind(bootcomp, c(max(assignments[,1]), min(bootgmyc), max(bootgmyc), mean(bootgmyc), sd(bootgmyc)))},
	 error = function(e) {print("There was an error message.")})

  },
  error = function(e) {print("There was an error message.")})

#	print(i)
#	i<-i+1
	
#
#compare number of delimited entities

#plot(hist(bootgmyc), xlim=c(8,22), main="8 spp (green), 5 seq/sp, 1 tip/entity (gray), 
#	full gymc with 40 tips, 22 entities (red)")
#arrows(x0=8, y0=5, y1=0, col="green")
#arrows(x0=max(assignments[,1]), y0=5, y1=0, col="red")

#save max, min and mode of bootstrapped tree

print(j)
j=j+1

}

write.table(bootcomp, file="bootcomp_1tip_theta5_tau2_10seq")

#sampling 2 tips per delimited entity

#simulated tree
#simtrees<-ape::read.nexus(file="simtrees.nex")

#analyze each simulated tree

bootcomp<-c()

for (j in 1:100) {

 simtree<-simtrees[[j]]

 # force ultrametric tree
 simtree<-phytools::force.ultrametric(simtree, method="extend")

 #do gmyc
 fullgmyc<-splits::gmyc(simtree)
 
 tryCatch(
 expr = {

 assignments<-splits::spec.list(fullgmyc)

 bootgmyc<-c()

 for (i in 1:100) {

	toretain<-c()
	for (k in 1:max(assignments[,1])) {
		if (dim(assignments[assignments[,1]==k,])[1] > 2) {
			toretain<-c(toretain, as.numeric(sample(assignments[assignments[,1]==k,2],2, replace=FALSE)))
			} else {
			toretain<-c(toretain, as.numeric(assignments[assignments[,1]==k,][,2]))
		}	
		k=k+1
	}
	
	#remove branches
	todrop<-setdiff(c(1:length(simtree$tip.label)), toretain)
	bootree<-ape::drop.tip(simtree, as.character(todrop))

	#do gmyc
	tryCatch(
     expr = {newassign<-splits::spec.list(splits::gmyc(bootree))
             bootgmyc<-c(bootgmyc, max(newassign[,1]))},
	 error = function(e) {print("There was an error message.")})

 i=i+1
 }

#compare number of delimited entities

#plot(hist(bootgmyc), xlim=c(0,25), main="8 spp (green), 5 seq/sp, 2 tips/entity (gray), 
#	full gymc with 40 tips, 22 entities (red)")
#arrows(x0=8, y0=20,y1=0, col="green")
#arrows(x0=max(assignments[,1]), y0=20,y1=0, col="red")

bootcomp<-rbind(bootcomp, c(max(assignments[,1]), min(bootgmyc), max(bootgmyc), mean(bootgmyc), sd(bootgmyc)))

},
error = function(e) {print("There was an error message.")})

print(j)
j=j+1

}

write.table(bootcomp, file="bootcomp_2tip_theta5_tau2_10seq")

#sampling 3 tips per delimited entity
#analyze each simulated tree

bootcomp<-c()

for (j in 1:100) {

 simtree<-simtrees[[j]]

 # force ultrametric tree
 simtree<-phytools::force.ultrametric(simtree, method="extend")

 #do gmyc
 fullgmyc<-splits::gmyc(simtree)

 tryCatch(
 expr = {

 assignments<-splits::spec.list(fullgmyc)

 bootgmyc<-c()

 for (i in 1:100) {

	toretain<-c()
	for (k in 1:max(assignments[,1])) {
		if (dim(assignments[assignments[,1]==k,])[1] > 3) {
			toretain<-c(toretain, as.numeric(sample(assignments[assignments[,1]==k,2],3, replace=FALSE)))
			} else {
			toretain<-c(toretain, as.numeric(assignments[assignments[,1]==k,][,2]))
		}	
		k=k+1
	}
	
	#remove branches
	todrop<-setdiff(c(1:length(simtree$tip.label)), toretain)
	bootree<-ape::drop.tip(simtree, todrop)

	#do gmyc
	tryCatch(
     expr = {newassign<-splits::spec.list(splits::gmyc(bootree))
             bootgmyc<-c(bootgmyc, max(newassign[,1]))},
	 error = function(e) {print("There was an error message.")})

 i=i+1
 }

#compare number of delimited entities

#plot(hist(bootgmyc), xlim=c(0,25), main="8 spp (green), 5 seq/sp, 2 tips/entity (gray), 
#	full gymc with 40 tips, 22 entities (red)")
#arrows(x0=8, y0=20,y1=0, col="green")
#arrows(x0=max(assignments[,1]), y0=20,y1=0, col="red")

bootcomp<-rbind(bootcomp, c(max(assignments[,1]), min(bootgmyc), max(bootgmyc), mean(bootgmyc), sd(bootgmyc)))

},
error = function(e) {print("There was an error message.")})

print(j)
j=j+1

}

write.table(bootcomp, file="bootcomp_3tip_theta5_tau2_10seq")

#sampling 4 tips per delimited entity
#analyze each simulated tree

bootcomp<-c()

for (j in 1:100) {

 simtree<-simtrees[[j]]

 # force ultrametric tree
 simtree<-phytools::force.ultrametric(simtree, method="extend")

 #do gmyc
 fullgmyc<-splits::gmyc(simtree)

 tryCatch(
 expr = {

 assignments<-splits::spec.list(fullgmyc)

 bootgmyc<-c()

 for (i in 1:100) {

	toretain<-c()
	for (k in 1:max(assignments[,1])) {
		if (dim(assignments[assignments[,1]==k,])[1] > 4) {
			toretain<-c(toretain, as.numeric(sample(assignments[assignments[,1]==k,2],4, replace=FALSE)))
			} else {
			toretain<-c(toretain, as.numeric(assignments[assignments[,1]==k,][,2]))
		}	
		k=k+1
	}
	
	#remove branches
	todrop<-setdiff(c(1:length(simtree$tip.label)), toretain)
	bootree<-ape::drop.tip(simtree, todrop)

	#do gmyc
	tryCatch(
     expr = {newassign<-splits::spec.list(splits::gmyc(bootree))
             bootgmyc<-c(bootgmyc, max(newassign[,1]))},
	 error = function(e) {print("There was an error message.")})

 i=i+1
 }


#compare number of delimited entities

#plot(hist(bootgmyc), xlim=c(0,25), main="8 spp (green), 5 seq/sp, 2 tips/entity (gray), 
#	full gymc with 40 tips, 22 entities (red)")
#arrows(x0=8, y0=20,y1=0, col="green")
#arrows(x0=max(assignments[,1]), y0=20,y1=0, col="red")

bootcomp<-rbind(bootcomp, c(max(assignments[,1]), min(bootgmyc), max(bootgmyc), mean(bootgmyc), sd(bootgmyc)))

},
error = function(e) {print("There was an error message.")})

print(j)
j=j+1

}

write.table(bootcomp, file="bootcomp_4tip_theta5_tau2_10seq")

#sampling 6 tips per delimited entity
#analyze each simulated tree

bootcomp<-c()

for (j in 1:100) {

 simtree<-simtrees[[j]]

 # force ultrametric tree
 simtree<-phytools::force.ultrametric(simtree, method="extend")

 #do gmyc
 fullgmyc<-splits::gmyc(simtree)

 tryCatch(
 expr = {

 assignments<-splits::spec.list(fullgmyc)

 bootgmyc<-c()

 for (i in 1:100) {

	toretain<-c()
	for (k in 1:max(assignments[,1])) {
		if (dim(assignments[assignments[,1]==k,])[1] > 4) {
			toretain<-c(toretain, as.numeric(sample(assignments[assignments[,1]==k,2],6, replace=FALSE)))
			} else {
			toretain<-c(toretain, as.numeric(assignments[assignments[,1]==k,][,2]))
		}	
		k=k+1
	}
	
	#remove branches
	todrop<-setdiff(c(1:length(simtree$tip.label)), toretain)
	bootree<-ape::drop.tip(simtree, todrop)

	#do gmyc
	tryCatch(
     expr = {newassign<-splits::spec.list(splits::gmyc(bootree))
             bootgmyc<-c(bootgmyc, max(newassign[,1]))},
	 error = function(e) {print("There was an error message.")})

 i=i+1
 }


#compare number of delimited entities

#plot(hist(bootgmyc), xlim=c(0,25), main="8 spp (green), 5 seq/sp, 2 tips/entity (gray), 
#	full gymc with 40 tips, 22 entities (red)")
#arrows(x0=8, y0=20,y1=0, col="green")
#arrows(x0=max(assignments[,1]), y0=20,y1=0, col="red")

bootcomp<-rbind(bootcomp, c(max(assignments[,1]), min(bootgmyc), max(bootgmyc), mean(bootgmyc), sd(bootgmyc)))

},
error = function(e) {print("There was an error message.")})

print(j)
j=j+1

}

write.table(bootcomp, file="bootcomp_6tip_theta5_tau2_10seq")

#sampling 8 tips per delimited entity
#analyze each simulated tree

bootcomp<-c()

for (j in 1:100) {

 simtree<-simtrees[[j]]

 # force ultrametric tree
 simtree<-phytools::force.ultrametric(simtree, method="extend")

 #do gmyc
 fullgmyc<-splits::gmyc(simtree)

 tryCatch(
 expr = {

 assignments<-splits::spec.list(fullgmyc)

 bootgmyc<-c()

 for (i in 1:100) {

	toretain<-c()
	for (k in 1:max(assignments[,1])) {
		if (dim(assignments[assignments[,1]==k,])[1] > 4) {
			toretain<-c(toretain, as.numeric(sample(assignments[assignments[,1]==k,2],8, replace=FALSE)))
			} else {
			toretain<-c(toretain, as.numeric(assignments[assignments[,1]==k,][,2]))
		}	
		k=k+1
	}
	
	#remove branches
	todrop<-setdiff(c(1:length(simtree$tip.label)), toretain)
	bootree<-ape::drop.tip(simtree, todrop)

	#do gmyc
	tryCatch(
     expr = {newassign<-splits::spec.list(splits::gmyc(bootree))
             bootgmyc<-c(bootgmyc, max(newassign[,1]))},
	 error = function(e) {print("There was an error message.")})

 i=i+1
 }


#compare number of delimited entities

#plot(hist(bootgmyc), xlim=c(0,25), main="8 spp (green), 5 seq/sp, 2 tips/entity (gray), 
#	full gymc with 40 tips, 22 entities (red)")
#arrows(x0=8, y0=20,y1=0, col="green")
#arrows(x0=max(assignments[,1]), y0=20,y1=0, col="red")

bootcomp<-rbind(bootcomp, c(max(assignments[,1]), min(bootgmyc), max(bootgmyc), mean(bootgmyc), sd(bootgmyc)))

},
error = function(e) {print("There was an error message.")})

print(j)
j=j+1

}

write.table(bootcomp, file="bootcomp_8tip_theta5_tau2_10seq")

############################################################

#simulated tree
simtrees<-ape::read.nexus(file="simtrees_theta5_tau02_10seq.nex")

#analyze each simulated tree

bootcomp<-c()

for (j in 1:1000) {

 simtree<-simtrees[[j]]

 # force ultrametric tree
 simtree<-phytools::force.ultrametric(simtree, method="extend")

 #do gmyc
 fullgmyc<-splits::gmyc(simtree)
 
 tryCatch(
 expr = {
 assignments<-splits::spec.list(fullgmyc)

 #sampling one tip per delimited entity
 bootgmyc<-c()

# for (i in 1:1) {

	toretain<-c()
	for (k in 1:max(assignments[,1])) {
	 toretain<-c(toretain, as.numeric(sample(assignments[assignments[,1]==k,2],1)))
	 k=k+1
	}
	
	#remove branches
	todrop<-setdiff(c(1:length(simtree$tip.label)), toretain)
	bootree<-ape::drop.tip(simtree, as.character(todrop))

	#do gmyc
	
	tryCatch(
     expr = {newassign<-splits::spec.list(splits::gmyc(bootree))
             bootgmyc<-c(bootgmyc, max(newassign[,1]))
             bootcomp<-rbind(bootcomp, c(max(assignments[,1]), min(bootgmyc), max(bootgmyc), mean(bootgmyc), sd(bootgmyc)))},
	 error = function(e) {print("There was an error message.")})

  },
  error = function(e) {print("There was an error message.")})

#	print(i)
#	i<-i+1
	
#
#compare number of delimited entities

#plot(hist(bootgmyc), xlim=c(8,22), main="8 spp (green), 5 seq/sp, 1 tip/entity (gray), 
#	full gymc with 40 tips, 22 entities (red)")
#arrows(x0=8, y0=5, y1=0, col="green")
#arrows(x0=max(assignments[,1]), y0=5, y1=0, col="red")

#save max, min and mode of bootstrapped tree

print(j)
j=j+1

}

write.table(bootcomp, file="bootcomp_1tip_theta5_tau02_10seq")

#sampling 2 tips per delimited entity

#simulated tree
#simtrees<-ape::read.nexus(file="simtrees.nex")

#analyze each simulated tree

bootcomp<-c()

for (j in 1:100) {

 simtree<-simtrees[[j]]

 # force ultrametric tree
 simtree<-phytools::force.ultrametric(simtree, method="extend")

 #do gmyc
 fullgmyc<-splits::gmyc(simtree)
 
 tryCatch(
 expr = {

 assignments<-splits::spec.list(fullgmyc)

 bootgmyc<-c()

 for (i in 1:100) {

	toretain<-c()
	for (k in 1:max(assignments[,1])) {
		if (dim(assignments[assignments[,1]==k,])[1] > 2) {
			toretain<-c(toretain, as.numeric(sample(assignments[assignments[,1]==k,2],2, replace=FALSE)))
			} else {
			toretain<-c(toretain, as.numeric(assignments[assignments[,1]==k,][,2]))
		}	
		k=k+1
	}
	
	#remove branches
	todrop<-setdiff(c(1:length(simtree$tip.label)), toretain)
	bootree<-ape::drop.tip(simtree, as.character(todrop))

	#do gmyc
	tryCatch(
     expr = {newassign<-splits::spec.list(splits::gmyc(bootree))
             bootgmyc<-c(bootgmyc, max(newassign[,1]))},
	 error = function(e) {print("There was an error message.")})

 i=i+1
 }

#compare number of delimited entities

#plot(hist(bootgmyc), xlim=c(0,25), main="8 spp (green), 5 seq/sp, 2 tips/entity (gray), 
#	full gymc with 40 tips, 22 entities (red)")
#arrows(x0=8, y0=20,y1=0, col="green")
#arrows(x0=max(assignments[,1]), y0=20,y1=0, col="red")

bootcomp<-rbind(bootcomp, c(max(assignments[,1]), min(bootgmyc), max(bootgmyc), mean(bootgmyc), sd(bootgmyc)))

},
error = function(e) {print("There was an error message.")})

print(j)
j=j+1

}

write.table(bootcomp, file="bootcomp_2tip_theta5_tau02_10seq")

#sampling 3 tips per delimited entity
#analyze each simulated tree

bootcomp<-c()

for (j in 1:100) {

 simtree<-simtrees[[j]]

 # force ultrametric tree
 simtree<-phytools::force.ultrametric(simtree, method="extend")

 #do gmyc
 fullgmyc<-splits::gmyc(simtree)

 tryCatch(
 expr = {

 assignments<-splits::spec.list(fullgmyc)

 bootgmyc<-c()

 for (i in 1:100) {

	toretain<-c()
	for (k in 1:max(assignments[,1])) {
		if (dim(assignments[assignments[,1]==k,])[1] > 3) {
			toretain<-c(toretain, as.numeric(sample(assignments[assignments[,1]==k,2],3, replace=FALSE)))
			} else {
			toretain<-c(toretain, as.numeric(assignments[assignments[,1]==k,][,2]))
		}	
		k=k+1
	}
	
	#remove branches
	todrop<-setdiff(c(1:length(simtree$tip.label)), toretain)
	bootree<-ape::drop.tip(simtree, todrop)

	#do gmyc
	tryCatch(
     expr = {newassign<-splits::spec.list(splits::gmyc(bootree))
             bootgmyc<-c(bootgmyc, max(newassign[,1]))},
	 error = function(e) {print("There was an error message.")})

 i=i+1
 }

#compare number of delimited entities

#plot(hist(bootgmyc), xlim=c(0,25), main="8 spp (green), 5 seq/sp, 2 tips/entity (gray), 
#	full gymc with 40 tips, 22 entities (red)")
#arrows(x0=8, y0=20,y1=0, col="green")
#arrows(x0=max(assignments[,1]), y0=20,y1=0, col="red")

bootcomp<-rbind(bootcomp, c(max(assignments[,1]), min(bootgmyc), max(bootgmyc), mean(bootgmyc), sd(bootgmyc)))

},
error = function(e) {print("There was an error message.")})

print(j)
j=j+1

}

write.table(bootcomp, file="bootcomp_3tip_theta5_tau02_10seq")

#sampling 4 tips per delimited entity
#analyze each simulated tree

bootcomp<-c()

for (j in 1:100) {

 simtree<-simtrees[[j]]

 # force ultrametric tree
 simtree<-phytools::force.ultrametric(simtree, method="extend")

 #do gmyc
 fullgmyc<-splits::gmyc(simtree)

 tryCatch(
 expr = {

 assignments<-splits::spec.list(fullgmyc)

 bootgmyc<-c()

 for (i in 1:100) {

	toretain<-c()
	for (k in 1:max(assignments[,1])) {
		if (dim(assignments[assignments[,1]==k,])[1] > 4) {
			toretain<-c(toretain, as.numeric(sample(assignments[assignments[,1]==k,2],4, replace=FALSE)))
			} else {
			toretain<-c(toretain, as.numeric(assignments[assignments[,1]==k,][,2]))
		}	
		k=k+1
	}
	
	#remove branches
	todrop<-setdiff(c(1:length(simtree$tip.label)), toretain)
	bootree<-ape::drop.tip(simtree, todrop)

	#do gmyc
	tryCatch(
     expr = {newassign<-splits::spec.list(splits::gmyc(bootree))
             bootgmyc<-c(bootgmyc, max(newassign[,1]))},
	 error = function(e) {print("There was an error message.")})

 i=i+1
 }


#compare number of delimited entities

#plot(hist(bootgmyc), xlim=c(0,25), main="8 spp (green), 5 seq/sp, 2 tips/entity (gray), 
#	full gymc with 40 tips, 22 entities (red)")
#arrows(x0=8, y0=20,y1=0, col="green")
#arrows(x0=max(assignments[,1]), y0=20,y1=0, col="red")

bootcomp<-rbind(bootcomp, c(max(assignments[,1]), min(bootgmyc), max(bootgmyc), mean(bootgmyc), sd(bootgmyc)))

},
error = function(e) {print("There was an error message.")})

print(j)
j=j+1

}

write.table(bootcomp, file="bootcomp_4tip_theta5_tau02_10seq")

#sampling 6 tips per delimited entity
#analyze each simulated tree

bootcomp<-c()

for (j in 1:100) {

 simtree<-simtrees[[j]]

 # force ultrametric tree
 simtree<-phytools::force.ultrametric(simtree, method="extend")

 #do gmyc
 fullgmyc<-splits::gmyc(simtree)

 tryCatch(
 expr = {

 assignments<-splits::spec.list(fullgmyc)

 bootgmyc<-c()

 for (i in 1:100) {

	toretain<-c()
	for (k in 1:max(assignments[,1])) {
		if (dim(assignments[assignments[,1]==k,])[1] > 4) {
			toretain<-c(toretain, as.numeric(sample(assignments[assignments[,1]==k,2],6, replace=FALSE)))
			} else {
			toretain<-c(toretain, as.numeric(assignments[assignments[,1]==k,][,2]))
		}	
		k=k+1
	}
	
	#remove branches
	todrop<-setdiff(c(1:length(simtree$tip.label)), toretain)
	bootree<-ape::drop.tip(simtree, todrop)

	#do gmyc
	tryCatch(
     expr = {newassign<-splits::spec.list(splits::gmyc(bootree))
             bootgmyc<-c(bootgmyc, max(newassign[,1]))},
	 error = function(e) {print("There was an error message.")})

 i=i+1
 }


#compare number of delimited entities

#plot(hist(bootgmyc), xlim=c(0,25), main="8 spp (green), 5 seq/sp, 2 tips/entity (gray), 
#	full gymc with 40 tips, 22 entities (red)")
#arrows(x0=8, y0=20,y1=0, col="green")
#arrows(x0=max(assignments[,1]), y0=20,y1=0, col="red")

bootcomp<-rbind(bootcomp, c(max(assignments[,1]), min(bootgmyc), max(bootgmyc), mean(bootgmyc), sd(bootgmyc)))

},
error = function(e) {print("There was an error message.")})

print(j)
j=j+1

}

write.table(bootcomp, file="bootcomp_6tip_theta5_tau02_10seq")

#sampling 8 tips per delimited entity
#analyze each simulated tree

bootcomp<-c()

for (j in 1:100) {

 simtree<-simtrees[[j]]

 # force ultrametric tree
 simtree<-phytools::force.ultrametric(simtree, method="extend")

 #do gmyc
 fullgmyc<-splits::gmyc(simtree)

 tryCatch(
 expr = {

 assignments<-splits::spec.list(fullgmyc)

 bootgmyc<-c()

 for (i in 1:100) {

	toretain<-c()
	for (k in 1:max(assignments[,1])) {
		if (dim(assignments[assignments[,1]==k,])[1] > 4) {
			toretain<-c(toretain, as.numeric(sample(assignments[assignments[,1]==k,2],8, replace=FALSE)))
			} else {
			toretain<-c(toretain, as.numeric(assignments[assignments[,1]==k,][,2]))
		}	
		k=k+1
	}
	
	#remove branches
	todrop<-setdiff(c(1:length(simtree$tip.label)), toretain)
	bootree<-ape::drop.tip(simtree, todrop)

	#do gmyc
	tryCatch(
     expr = {newassign<-splits::spec.list(splits::gmyc(bootree))
             bootgmyc<-c(bootgmyc, max(newassign[,1]))},
	 error = function(e) {print("There was an error message.")})

 i=i+1
 }


#compare number of delimited entities

#plot(hist(bootgmyc), xlim=c(0,25), main="8 spp (green), 5 seq/sp, 2 tips/entity (gray), 
#	full gymc with 40 tips, 22 entities (red)")
#arrows(x0=8, y0=20,y1=0, col="green")
#arrows(x0=max(assignments[,1]), y0=20,y1=0, col="red")

bootcomp<-rbind(bootcomp, c(max(assignments[,1]), min(bootgmyc), max(bootgmyc), mean(bootgmyc), sd(bootgmyc)))

},
error = function(e) {print("There was an error message.")})

print(j)
j=j+1

}

write.table(bootcomp, file="bootcomp_8tip_theta5_tau02_10seq")


############################################################

#simulated tree
simtrees<-ape::read.nexus(file="simtrees_theta50_tau2_10seq.nex")

#analyze each simulated tree

bootcomp<-c()

for (j in 1:1000) {

 simtree<-simtrees[[j]]

 # force ultrametric tree
 simtree<-phytools::force.ultrametric(simtree, method="extend")

 #do gmyc
 fullgmyc<-splits::gmyc(simtree)
 
 tryCatch(
 expr = {
 assignments<-splits::spec.list(fullgmyc)

 #sampling one tip per delimited entity
 bootgmyc<-c()

# for (i in 1:1) {

	toretain<-c()
	for (k in 1:max(assignments[,1])) {
	 toretain<-c(toretain, as.numeric(sample(assignments[assignments[,1]==k,2],1)))
	 k=k+1
	}
	
	#remove branches
	todrop<-setdiff(c(1:length(simtree$tip.label)), toretain)
	bootree<-ape::drop.tip(simtree, as.character(todrop))

	#do gmyc
	
	tryCatch(
     expr = {newassign<-splits::spec.list(splits::gmyc(bootree))
             bootgmyc<-c(bootgmyc, max(newassign[,1]))
             bootcomp<-rbind(bootcomp, c(max(assignments[,1]), min(bootgmyc), max(bootgmyc), mean(bootgmyc), sd(bootgmyc)))},
	 error = function(e) {print("There was an error message.")})

  },
  error = function(e) {print("There was an error message.")})

#	print(i)
#	i<-i+1
	
#
#compare number of delimited entities

#plot(hist(bootgmyc), xlim=c(8,22), main="8 spp (green), 5 seq/sp, 1 tip/entity (gray), 
#	full gymc with 40 tips, 22 entities (red)")
#arrows(x0=8, y0=5, y1=0, col="green")
#arrows(x0=max(assignments[,1]), y0=5, y1=0, col="red")

#save max, min and mode of bootstrapped tree

print(j)
j=j+1

}

write.table(bootcomp, file="bootcomp_1tip_theta50_tau2_10seq")

#sampling 2 tips per delimited entity

#simulated tree
#simtrees<-ape::read.nexus(file="simtrees.nex")

#analyze each simulated tree

bootcomp<-c()

for (j in 1:100) {

 simtree<-simtrees[[j]]

 # force ultrametric tree
 simtree<-phytools::force.ultrametric(simtree, method="extend")

 #do gmyc
 fullgmyc<-splits::gmyc(simtree)
 
 tryCatch(
 expr = {

 assignments<-splits::spec.list(fullgmyc)

 bootgmyc<-c()

 for (i in 1:100) {

	toretain<-c()
	for (k in 1:max(assignments[,1])) {
		if (dim(assignments[assignments[,1]==k,])[1] > 2) {
			toretain<-c(toretain, as.numeric(sample(assignments[assignments[,1]==k,2],2, replace=FALSE)))
			} else {
			toretain<-c(toretain, as.numeric(assignments[assignments[,1]==k,][,2]))
		}	
		k=k+1
	}
	
	#remove branches
	todrop<-setdiff(c(1:length(simtree$tip.label)), toretain)
	bootree<-ape::drop.tip(simtree, as.character(todrop))

	#do gmyc
	tryCatch(
     expr = {newassign<-splits::spec.list(splits::gmyc(bootree))
             bootgmyc<-c(bootgmyc, max(newassign[,1]))},
	 error = function(e) {print("There was an error message.")})

 i=i+1
 }

#compare number of delimited entities

#plot(hist(bootgmyc), xlim=c(0,25), main="8 spp (green), 5 seq/sp, 2 tips/entity (gray), 
#	full gymc with 40 tips, 22 entities (red)")
#arrows(x0=8, y0=20,y1=0, col="green")
#arrows(x0=max(assignments[,1]), y0=20,y1=0, col="red")

bootcomp<-rbind(bootcomp, c(max(assignments[,1]), min(bootgmyc), max(bootgmyc), mean(bootgmyc), sd(bootgmyc)))

},
error = function(e) {print("There was an error message.")})

print(j)
j=j+1

}

write.table(bootcomp, file="bootcomp_2tip_theta50_tau2_10seq")

#sampling 3 tips per delimited entity
#analyze each simulated tree

bootcomp<-c()

for (j in 1:100) {

 simtree<-simtrees[[j]]

 # force ultrametric tree
 simtree<-phytools::force.ultrametric(simtree, method="extend")

 #do gmyc
 fullgmyc<-splits::gmyc(simtree)

 tryCatch(
 expr = {

 assignments<-splits::spec.list(fullgmyc)

 bootgmyc<-c()

 for (i in 1:100) {

	toretain<-c()
	for (k in 1:max(assignments[,1])) {
		if (dim(assignments[assignments[,1]==k,])[1] > 3) {
			toretain<-c(toretain, as.numeric(sample(assignments[assignments[,1]==k,2],3, replace=FALSE)))
			} else {
			toretain<-c(toretain, as.numeric(assignments[assignments[,1]==k,][,2]))
		}	
		k=k+1
	}
	
	#remove branches
	todrop<-setdiff(c(1:length(simtree$tip.label)), toretain)
	bootree<-ape::drop.tip(simtree, todrop)

	#do gmyc
	tryCatch(
     expr = {newassign<-splits::spec.list(splits::gmyc(bootree))
             bootgmyc<-c(bootgmyc, max(newassign[,1]))},
	 error = function(e) {print("There was an error message.")})

 i=i+1
 }

#compare number of delimited entities

#plot(hist(bootgmyc), xlim=c(0,25), main="8 spp (green), 5 seq/sp, 2 tips/entity (gray), 
#	full gymc with 40 tips, 22 entities (red)")
#arrows(x0=8, y0=20,y1=0, col="green")
#arrows(x0=max(assignments[,1]), y0=20,y1=0, col="red")

bootcomp<-rbind(bootcomp, c(max(assignments[,1]), min(bootgmyc), max(bootgmyc), mean(bootgmyc), sd(bootgmyc)))

},
error = function(e) {print("There was an error message.")})

print(j)
j=j+1

}

write.table(bootcomp, file="bootcomp_3tip_theta50_tau2_10seq")

#sampling 4 tips per delimited entity
#analyze each simulated tree

bootcomp<-c()

for (j in 1:100) {

 simtree<-simtrees[[j]]

 # force ultrametric tree
 simtree<-phytools::force.ultrametric(simtree, method="extend")

 #do gmyc
 fullgmyc<-splits::gmyc(simtree)

 tryCatch(
 expr = {

 assignments<-splits::spec.list(fullgmyc)

 bootgmyc<-c()

 for (i in 1:100) {

	toretain<-c()
	for (k in 1:max(assignments[,1])) {
		if (dim(assignments[assignments[,1]==k,])[1] > 4) {
			toretain<-c(toretain, as.numeric(sample(assignments[assignments[,1]==k,2],4, replace=FALSE)))
			} else {
			toretain<-c(toretain, as.numeric(assignments[assignments[,1]==k,][,2]))
		}	
		k=k+1
	}
	
	#remove branches
	todrop<-setdiff(c(1:length(simtree$tip.label)), toretain)
	bootree<-ape::drop.tip(simtree, todrop)

	#do gmyc
	tryCatch(
     expr = {newassign<-splits::spec.list(splits::gmyc(bootree))
             bootgmyc<-c(bootgmyc, max(newassign[,1]))},
	 error = function(e) {print("There was an error message.")})

 i=i+1
 }


#compare number of delimited entities

#plot(hist(bootgmyc), xlim=c(0,25), main="8 spp (green), 5 seq/sp, 2 tips/entity (gray), 
#	full gymc with 40 tips, 22 entities (red)")
#arrows(x0=8, y0=20,y1=0, col="green")
#arrows(x0=max(assignments[,1]), y0=20,y1=0, col="red")

bootcomp<-rbind(bootcomp, c(max(assignments[,1]), min(bootgmyc), max(bootgmyc), mean(bootgmyc), sd(bootgmyc)))

},
error = function(e) {print("There was an error message.")})

print(j)
j=j+1

}

write.table(bootcomp, file="bootcomp_4tip_theta50_tau2_10seq")

#sampling 6 tips per delimited entity
#analyze each simulated tree

bootcomp<-c()

for (j in 1:100) {

 simtree<-simtrees[[j]]

 # force ultrametric tree
 simtree<-phytools::force.ultrametric(simtree, method="extend")

 #do gmyc
 fullgmyc<-splits::gmyc(simtree)

 tryCatch(
 expr = {

 assignments<-splits::spec.list(fullgmyc)

 bootgmyc<-c()

 for (i in 1:100) {

	toretain<-c()
	for (k in 1:max(assignments[,1])) {
		if (dim(assignments[assignments[,1]==k,])[1] > 4) {
			toretain<-c(toretain, as.numeric(sample(assignments[assignments[,1]==k,2],6, replace=FALSE)))
			} else {
			toretain<-c(toretain, as.numeric(assignments[assignments[,1]==k,][,2]))
		}	
		k=k+1
	}
	
	#remove branches
	todrop<-setdiff(c(1:length(simtree$tip.label)), toretain)
	bootree<-ape::drop.tip(simtree, todrop)

	#do gmyc
	tryCatch(
     expr = {newassign<-splits::spec.list(splits::gmyc(bootree))
             bootgmyc<-c(bootgmyc, max(newassign[,1]))},
	 error = function(e) {print("There was an error message.")})

 i=i+1
 }


#compare number of delimited entities

#plot(hist(bootgmyc), xlim=c(0,25), main="8 spp (green), 5 seq/sp, 2 tips/entity (gray), 
#	full gymc with 40 tips, 22 entities (red)")
#arrows(x0=8, y0=20,y1=0, col="green")
#arrows(x0=max(assignments[,1]), y0=20,y1=0, col="red")

bootcomp<-rbind(bootcomp, c(max(assignments[,1]), min(bootgmyc), max(bootgmyc), mean(bootgmyc), sd(bootgmyc)))

},
error = function(e) {print("There was an error message.")})

print(j)
j=j+1

}

write.table(bootcomp, file="bootcomp_6tip_theta50_tau2_10seq")

#sampling 8 tips per delimited entity
#analyze each simulated tree

bootcomp<-c()

for (j in 1:100) {

 simtree<-simtrees[[j]]

 # force ultrametric tree
 simtree<-phytools::force.ultrametric(simtree, method="extend")

 #do gmyc
 fullgmyc<-splits::gmyc(simtree)

 tryCatch(
 expr = {

 assignments<-splits::spec.list(fullgmyc)

 bootgmyc<-c()

 for (i in 1:100) {

	toretain<-c()
	for (k in 1:max(assignments[,1])) {
		if (dim(assignments[assignments[,1]==k,])[1] > 4) {
			toretain<-c(toretain, as.numeric(sample(assignments[assignments[,1]==k,2],8, replace=FALSE)))
			} else {
			toretain<-c(toretain, as.numeric(assignments[assignments[,1]==k,][,2]))
		}	
		k=k+1
	}
	
	#remove branches
	todrop<-setdiff(c(1:length(simtree$tip.label)), toretain)
	bootree<-ape::drop.tip(simtree, todrop)

	#do gmyc
	tryCatch(
     expr = {newassign<-splits::spec.list(splits::gmyc(bootree))
             bootgmyc<-c(bootgmyc, max(newassign[,1]))},
	 error = function(e) {print("There was an error message.")})

 i=i+1
 }


#compare number of delimited entities

#plot(hist(bootgmyc), xlim=c(0,25), main="8 spp (green), 5 seq/sp, 2 tips/entity (gray), 
#	full gymc with 40 tips, 22 entities (red)")
#arrows(x0=8, y0=20,y1=0, col="green")
#arrows(x0=max(assignments[,1]), y0=20,y1=0, col="red")

bootcomp<-rbind(bootcomp, c(max(assignments[,1]), min(bootgmyc), max(bootgmyc), mean(bootgmyc), sd(bootgmyc)))

},
error = function(e) {print("There was an error message.")})

print(j)
j=j+1

}

write.table(bootcomp, file="bootcomp_8tip_theta50_tau2_10seq")


############################################################

#simulated tree
simtrees<-ape::read.nexus(file="simtrees_theta50_tau02_10seq.nex")

#analyze each simulated tree

bootcomp<-c()

for (j in 1:1000) {

 simtree<-simtrees[[j]]

 # force ultrametric tree
 simtree<-phytools::force.ultrametric(simtree, method="extend")

 #do gmyc
 fullgmyc<-splits::gmyc(simtree)
 
 tryCatch(
 expr = {
 assignments<-splits::spec.list(fullgmyc)

 #sampling one tip per delimited entity
 bootgmyc<-c()

# for (i in 1:1) {

	toretain<-c()
	for (k in 1:max(assignments[,1])) {
	 toretain<-c(toretain, as.numeric(sample(assignments[assignments[,1]==k,2],1)))
	 k=k+1
	}
	
	#remove branches
	todrop<-setdiff(c(1:length(simtree$tip.label)), toretain)
	bootree<-ape::drop.tip(simtree, as.character(todrop))

	#do gmyc
	
	tryCatch(
     expr = {newassign<-splits::spec.list(splits::gmyc(bootree))
             bootgmyc<-c(bootgmyc, max(newassign[,1]))
             bootcomp<-rbind(bootcomp, c(max(assignments[,1]), min(bootgmyc), max(bootgmyc), mean(bootgmyc), sd(bootgmyc)))},
	 error = function(e) {print("There was an error message.")})

  },
  error = function(e) {print("There was an error message.")})

#	print(i)
#	i<-i+1
	
#
#compare number of delimited entities

#plot(hist(bootgmyc), xlim=c(8,22), main="8 spp (green), 5 seq/sp, 1 tip/entity (gray), 
#	full gymc with 40 tips, 22 entities (red)")
#arrows(x0=8, y0=5, y1=0, col="green")
#arrows(x0=max(assignments[,1]), y0=5, y1=0, col="red")

#save max, min and mode of bootstrapped tree

print(j)
j=j+1

}

write.table(bootcomp, file="bootcomp_1tip_theta50_tau02_10seq")

#sampling 2 tips per delimited entity

#simulated tree
#simtrees<-ape::read.nexus(file="simtrees.nex")

#analyze each simulated tree

bootcomp<-c()

for (j in 1:100) {

 simtree<-simtrees[[j]]

 # force ultrametric tree
 simtree<-phytools::force.ultrametric(simtree, method="extend")

 #do gmyc
 fullgmyc<-splits::gmyc(simtree)
 
 tryCatch(
 expr = {

 assignments<-splits::spec.list(fullgmyc)

 bootgmyc<-c()

 for (i in 1:100) {

	toretain<-c()
	for (k in 1:max(assignments[,1])) {
		if (dim(assignments[assignments[,1]==k,])[1] > 2) {
			toretain<-c(toretain, as.numeric(sample(assignments[assignments[,1]==k,2],2, replace=FALSE)))
			} else {
			toretain<-c(toretain, as.numeric(assignments[assignments[,1]==k,][,2]))
		}	
		k=k+1
	}
	
	#remove branches
	todrop<-setdiff(c(1:length(simtree$tip.label)), toretain)
	bootree<-ape::drop.tip(simtree, as.character(todrop))

	#do gmyc
	tryCatch(
     expr = {newassign<-splits::spec.list(splits::gmyc(bootree))
             bootgmyc<-c(bootgmyc, max(newassign[,1]))},
	 error = function(e) {print("There was an error message.")})

 i=i+1
 }

#compare number of delimited entities

#plot(hist(bootgmyc), xlim=c(0,25), main="8 spp (green), 5 seq/sp, 2 tips/entity (gray), 
#	full gymc with 40 tips, 22 entities (red)")
#arrows(x0=8, y0=20,y1=0, col="green")
#arrows(x0=max(assignments[,1]), y0=20,y1=0, col="red")

bootcomp<-rbind(bootcomp, c(max(assignments[,1]), min(bootgmyc), max(bootgmyc), mean(bootgmyc), sd(bootgmyc)))

},
error = function(e) {print("There was an error message.")})

print(j)
j=j+1

}

write.table(bootcomp, file="bootcomp_2tip_theta50_tau02_10seq")

#sampling 3 tips per delimited entity
#analyze each simulated tree

bootcomp<-c()

for (j in 1:100) {

 simtree<-simtrees[[j]]

 # force ultrametric tree
 simtree<-phytools::force.ultrametric(simtree, method="extend")

 #do gmyc
 fullgmyc<-splits::gmyc(simtree)

 tryCatch(
 expr = {

 assignments<-splits::spec.list(fullgmyc)

 bootgmyc<-c()

 for (i in 1:100) {

	toretain<-c()
	for (k in 1:max(assignments[,1])) {
		if (dim(assignments[assignments[,1]==k,])[1] > 3) {
			toretain<-c(toretain, as.numeric(sample(assignments[assignments[,1]==k,2],3, replace=FALSE)))
			} else {
			toretain<-c(toretain, as.numeric(assignments[assignments[,1]==k,][,2]))
		}	
		k=k+1
	}
	
	#remove branches
	todrop<-setdiff(c(1:length(simtree$tip.label)), toretain)
	bootree<-ape::drop.tip(simtree, todrop)

	#do gmyc
	tryCatch(
     expr = {newassign<-splits::spec.list(splits::gmyc(bootree))
             bootgmyc<-c(bootgmyc, max(newassign[,1]))},
	 error = function(e) {print("There was an error message.")})

 i=i+1
 }

#compare number of delimited entities

#plot(hist(bootgmyc), xlim=c(0,25), main="8 spp (green), 5 seq/sp, 2 tips/entity (gray), 
#	full gymc with 40 tips, 22 entities (red)")
#arrows(x0=8, y0=20,y1=0, col="green")
#arrows(x0=max(assignments[,1]), y0=20,y1=0, col="red")

bootcomp<-rbind(bootcomp, c(max(assignments[,1]), min(bootgmyc), max(bootgmyc), mean(bootgmyc), sd(bootgmyc)))

},
error = function(e) {print("There was an error message.")})

print(j)
j=j+1

}

write.table(bootcomp, file="bootcomp_3tip_theta50_tau02_10seq")

#sampling 4 tips per delimited entity
#analyze each simulated tree

bootcomp<-c()

for (j in 1:100) {

 simtree<-simtrees[[j]]

 # force ultrametric tree
 simtree<-phytools::force.ultrametric(simtree, method="extend")

 #do gmyc
 fullgmyc<-splits::gmyc(simtree)

 tryCatch(
 expr = {

 assignments<-splits::spec.list(fullgmyc)

 bootgmyc<-c()

 for (i in 1:100) {

	toretain<-c()
	for (k in 1:max(assignments[,1])) {
		if (dim(assignments[assignments[,1]==k,])[1] > 4) {
			toretain<-c(toretain, as.numeric(sample(assignments[assignments[,1]==k,2],4, replace=FALSE)))
			} else {
			toretain<-c(toretain, as.numeric(assignments[assignments[,1]==k,][,2]))
		}	
		k=k+1
	}
	
	#remove branches
	todrop<-setdiff(c(1:length(simtree$tip.label)), toretain)
	bootree<-ape::drop.tip(simtree, todrop)

	#do gmyc
	tryCatch(
     expr = {newassign<-splits::spec.list(splits::gmyc(bootree))
             bootgmyc<-c(bootgmyc, max(newassign[,1]))},
	 error = function(e) {print("There was an error message.")})

 i=i+1
 }


#compare number of delimited entities

#plot(hist(bootgmyc), xlim=c(0,25), main="8 spp (green), 5 seq/sp, 2 tips/entity (gray), 
#	full gymc with 40 tips, 22 entities (red)")
#arrows(x0=8, y0=20,y1=0, col="green")
#arrows(x0=max(assignments[,1]), y0=20,y1=0, col="red")

bootcomp<-rbind(bootcomp, c(max(assignments[,1]), min(bootgmyc), max(bootgmyc), mean(bootgmyc), sd(bootgmyc)))

},
error = function(e) {print("There was an error message.")})

print(j)
j=j+1

}

write.table(bootcomp, file="bootcomp_4tip_theta50_tau02_10seq")

#sampling 6 tips per delimited entity
#analyze each simulated tree

bootcomp<-c()

for (j in 1:100) {

 simtree<-simtrees[[j]]

 # force ultrametric tree
 simtree<-phytools::force.ultrametric(simtree, method="extend")

 #do gmyc
 fullgmyc<-splits::gmyc(simtree)

 tryCatch(
 expr = {

 assignments<-splits::spec.list(fullgmyc)

 bootgmyc<-c()

 for (i in 1:100) {

	toretain<-c()
	for (k in 1:max(assignments[,1])) {
		if (dim(assignments[assignments[,1]==k,])[1] > 4) {
			toretain<-c(toretain, as.numeric(sample(assignments[assignments[,1]==k,2],6, replace=FALSE)))
			} else {
			toretain<-c(toretain, as.numeric(assignments[assignments[,1]==k,][,2]))
		}	
		k=k+1
	}
	
	#remove branches
	todrop<-setdiff(c(1:length(simtree$tip.label)), toretain)
	bootree<-ape::drop.tip(simtree, todrop)

	#do gmyc
	tryCatch(
     expr = {newassign<-splits::spec.list(splits::gmyc(bootree))
             bootgmyc<-c(bootgmyc, max(newassign[,1]))},
	 error = function(e) {print("There was an error message.")})

 i=i+1
 }


#compare number of delimited entities

#plot(hist(bootgmyc), xlim=c(0,25), main="8 spp (green), 5 seq/sp, 2 tips/entity (gray), 
#	full gymc with 40 tips, 22 entities (red)")
#arrows(x0=8, y0=20,y1=0, col="green")
#arrows(x0=max(assignments[,1]), y0=20,y1=0, col="red")

bootcomp<-rbind(bootcomp, c(max(assignments[,1]), min(bootgmyc), max(bootgmyc), mean(bootgmyc), sd(bootgmyc)))

},
error = function(e) {print("There was an error message.")})

print(j)
j=j+1

}

write.table(bootcomp, file="bootcomp_6tip_theta50_tau02_10seq")

#sampling 8 tips per delimited entity
#analyze each simulated tree

bootcomp<-c()

for (j in 1:100) {

 simtree<-simtrees[[j]]

 # force ultrametric tree
 simtree<-phytools::force.ultrametric(simtree, method="extend")

 #do gmyc
 fullgmyc<-splits::gmyc(simtree)

 tryCatch(
 expr = {

 assignments<-splits::spec.list(fullgmyc)

 bootgmyc<-c()

 for (i in 1:100) {

	toretain<-c()
	for (k in 1:max(assignments[,1])) {
		if (dim(assignments[assignments[,1]==k,])[1] > 4) {
			toretain<-c(toretain, as.numeric(sample(assignments[assignments[,1]==k,2],8, replace=FALSE)))
			} else {
			toretain<-c(toretain, as.numeric(assignments[assignments[,1]==k,][,2]))
		}	
		k=k+1
	}
	
	#remove branches
	todrop<-setdiff(c(1:length(simtree$tip.label)), toretain)
	bootree<-ape::drop.tip(simtree, todrop)

	#do gmyc
	tryCatch(
     expr = {newassign<-splits::spec.list(splits::gmyc(bootree))
             bootgmyc<-c(bootgmyc, max(newassign[,1]))},
	 error = function(e) {print("There was an error message.")})

 i=i+1
 }


#compare number of delimited entities

#plot(hist(bootgmyc), xlim=c(0,25), main="8 spp (green), 5 seq/sp, 2 tips/entity (gray), 
#	full gymc with 40 tips, 22 entities (red)")
#arrows(x0=8, y0=20,y1=0, col="green")
#arrows(x0=max(assignments[,1]), y0=20,y1=0, col="red")

bootcomp<-rbind(bootcomp, c(max(assignments[,1]), min(bootgmyc), max(bootgmyc), mean(bootgmyc), sd(bootgmyc)))

},
error = function(e) {print("There was an error message.")})

print(j)
j=j+1

}

write.table(bootcomp, file="bootcomp_8tip_theta50_tau02_10seq")



