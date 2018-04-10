#R code to retrieve the 95% highest posterior distribution of of node ages for a single node given a posterior set of chronograms.  
#J. Schenk 9 April 2018 

library(ape)
library(coda)

#Function to estimate the age of a node given two members of a clade.
AgeDensity <- function(phy, species1, species2){
	NodeNumber <- mrca(phy)[species1, species2]
	ages <- branching.times(phy)[as.character(NodeNumber)]
	return(as.numeric(ages))
}

#read in tree
tree <- read.nexus(file="")

#Apply function across the posterior distribution to obtain node ages. This will take some time to run, depending on the number of trees.  It took me about 20 minutes.
NodeAges <- unlist(lapply(.uncompressTipLabel(tree), AgeDensity, "species1", "species2"))

#estimate the HPD interval for the node ages
HPD <- HPDinterval(as.mcmc(NodeAges), prob = 0.95)

#Density plot for 95% HPD
plot(density(as.numeric(NodeAges[which(NodeAges > HPD[1, 1] & NodeAges < HPD[1, 2])])), main="", las=1, xlab = "Million of years before present", lwd = 2)







