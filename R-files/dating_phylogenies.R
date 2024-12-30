library(ape)
# read in your tree (or trees)
# my file example "sample1 .nex" is actually a list of many trees, so I have an
# example of running the smoothing on one tree
setwd("~/dscam/beast_run")
dscam<-read.nexus("dscam.nex")
# Switch this (^) out with mitochondria as well  
# rerooting the tree onto to Odonata
plot(dscam, cex=0.6)
#here is the range of lambda values we'll use
l <- 10^(-1:6)
#here we are making a blank numeric list
cv <- numeric(length(l))

#we will use a loop to run the smoothing test (to see what lambda smoothing parameter we want)

i<-1
# this does the cross validation step for each parameter value
# here we are putting the value for each smoothing parameter into the blank list called cv we made
# this may take a while, it's reliant on the tips we have
# mine runs for a while since I have 152 tips.  In my chronopl I'm using a specific tree in a list of trees
# hence the test1$gen.#####.  It's listing a specific tree
# so the script is summing all the smoothing info with cross validation and putting it into the "cv" list
# these correspond to each lambda smoothing parameter we are using ('l' list)
for (i in 1:length(l)){
  print(i)
  cv[i] <- sum(attr(chronopl(dscam, lambda = l[i], CV=TRUE), "D2"))
}

# we then plot it to see the results of the CV for each lambda.  The lower value is better
# this is the result (low is better)
plot(l, cv)
cv
l

# now that we know the best lambda value, we can use the "chronos" function to date our tree.
# in my example, lambda=100 is best, so that is what I smoothed it with
# new tree
dscam<- chronos(dscam, lambda=1, model="correlated")
plot(dscam, cex=0.6)

# scaling to root age
# now that I have the three smoothed, I want to date it.  I know my root is 72 million years
# based on data from timetree.org
# so I am just taking my old tree and multiplying the edge.lengths by the root
# then I save it as a newick file
root <- 401
dscam$edge.length <- dscam$edge.length * root
plot(dscam, cex=0.6)
write.tree(new.tree1, file="dscam1.new")
