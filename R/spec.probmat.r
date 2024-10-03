#' @export

spec.probmat <- function(res) {

numtip <- length(res[[1]]$tree$tip.label)
probmat<-matrix(data=( rep(0, times=(numtip*numtip))), ncol=numtip, nrow=numtip)
rownames(probmat)<-res[[1]]$tree$tip.label
colnames(probmat)<-res[[1]]$tree$tip.label
ntrees<-length(res)
totalsamp<-ntrees*length(res[[1]]$par[,1])

for (q in 1:ntrees){
		for (j in 1:length(res[[q]]$par[,3])){
			assignlists<-res[[q]]$reassignment[[j]]
			assignlists<-assignlists[order(assignlists[,1]),]
			levs<-levels(as.factor(assignlists[,1]))
			for (i in 1:length(levs)){
				probmat[as.character(assignlists[which(assignlists[,1]==levs[i]),2]),as.character(assignlists[which(assignlists[,1]==levs[i]),2])]<-probmat[as.character(assignlists[which(assignlists[,1]==levs[i]),2]),as.character(assignlists[which(assignlists[,1]==levs[i]),2])]+1
				}
			}
		}

probmat<-probmat/totalsamp
class(probmat)<-"bgmycprobmat"
return(probmat)
}
