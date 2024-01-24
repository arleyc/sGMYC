#simulated tree
simtrees<-ape::read.nexus(file="simtrees.nex")

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
             bootcomp<-rbind(bootcomp, c(max(assignments[,1]), min(bootgmyc), max(bootgmyc), mean(bootgmyc)))},
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
		if (dim(assignments[assignments[,1]==k,])[1] > 4) {
			toretain<-c(toretain, as.numeric(sample(assignments[assignments[,1]==k,2],4, replace=FALSE)))
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

bootcomp<-rbind(bootcomp, c(max(assignments[,1]), min(bootgmyc), max(bootgmyc), mean(bootgmyc)))

},
error = function(e) {print("There was an error message.")})

print(j)
j=j+1

}

write.table(bootcomp, file="bootcomp_4tips_theta5_tau2_10seq")

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
		if (dim(assignments[assignments[,1]==k,])[1] > 6) {
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

bootcomp<-rbind(bootcomp, c(max(assignments[,1]), min(bootgmyc), max(bootgmyc), mean(bootgmyc)))

},
error = function(e) {print("There was an error message.")})

print(j)
j=j+1

}

write.table(bootcomp, file="bootcomp_6tips_theta5_tau2_10seq")

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
		if (dim(assignments[assignments[,1]==k,])[1] > 8) {
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

bootcomp<-rbind(bootcomp, c(max(assignments[,1]), min(bootgmyc), max(bootgmyc), mean(bootgmyc)))

},
error = function(e) {print("There was an error message.")})

print(j)
j=j+1

}

write.table(bootcomp, file="bootcomp_8tips_theta5_tau2_10seq")