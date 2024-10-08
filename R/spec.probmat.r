#' spec.probmat
#'
#' @description
#'   spec.probmat calculates a matrix of probabilities of conspecificity based on the
#'   original spec.promat function of the bGMYC package. When running sGMYC, the
#'   output list contains an table named 'reassignment', which assigns all samples to the
#'   the species delimited with a subsample of sequences.
#' 
#' @export
#' @param res A table containing the reassignment of sequences to the species delimited 
#' with a subsample of sequences when using sGMYC. This table is the element named
#' 'reassignment' in the list returned by sGMYC.
#'
#' @return
#' An object of class “matrix” with the probabilities of conspecificity between all samples
#' 
#' @author
#' Rafael F. Magalhães, MTT Santos & Arley Camargo
#'
#' @references
#' \itemize{
#'   \item Fujisawa T & Barraclough TG. 2013. Delimiting species using single-locus data and the Generalized Mixed Yule Coalescent approach: A revised method and evaluation on simulated data sets. Systematic Biology 62:707–724, doi: 10.1093/sysbio/syt033
#'   \item Pons J, Barraclough TG, Gomez-Zurita J, Cardoso A, Duran DP, Hazell S, Kamoun S, Sumlin WD & Vogler AP. 2006. Sequence-based species delimitation for the DNA taxonomy of undescribed insects. Systematic Biology 55:595-609, doi: 10.1080/10635150600852011.
#'   \item Reid NM & Carstens BC. 2012. Phylogenetic estimation error can decrease the accuracy of species delimitation: a Bayesian implementation of the general mixed Yule-coalescent model. BMC Evolutionary Biology 12:196, doi: 10.1186/1471-2148-12-196.
#'}
#'
#' @examples
#' data(hypotree)
#' myres<-sGMYC(hypotree, subsamp=2, nreps=100)
#' myprobmat<-spec.probmat(myres)
#' library(bGMYC)
#' plot.bgmycprobmat(myprobmat,hypotree)

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
