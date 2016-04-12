#G3tYourMa5t3r!!! Yaarctopus

#library statements
require(igraph)
require(lpbrim)
require(bipartite)

#funtions that used later
#function #1 that create links
link.create<-function(ind.vec,comm.ind,p.in,p.out,n){
  link.list<-NULL
  gp<-comm.ind[ind.vec==1]
  pa<-comm.ind[ind.vec==0]
  ngp<-sum(ind.vec)
  npa<-length(ind.vec)-ngp
  in.vec<-c(outer(gp,pa,"=="))
  prob.vec<-rep(0,ngp*npa)
  prob.vec[in.vec]<-p.in
  prob.vec[!in.vec]<-p.out
  for(i in 1:n){
    #randomly link
    temp<-rbinom(ngp*npa,1,prob.vec)
    temp<-matrix(temp,nrow=ngp,ncol=npa)
    link.list[[length(link.list)+1]]<-temp
  }
  link.list
}

#function #2 compute the normalized mutual information between two GP S-matrix
norm.mut.infor<-function(s.mat.gp.true,s.mat.gp){
    p<-t(s.mat.gp.true)%*%s.mat.gp/nrow(s.mat.gp.true)
    px<-apply(p,1,"sum")
    py<-apply(p,2,"sum")
    logpp<-log(p/(px%*%t(py)))
    logpp[p==0]<-0
    Ixy<-sum(logpp*p)
    logpx<-log(px)
    logpx[px==0]<-0
    Hx<--sum(logpx*px)
    logpy<-log(py)
    logpy[py==0]<-0
    Hy<--sum(logpy*py)
    norm.info<-2*Ixy/(Hx+Hy)
    norm.info
}

#function #3 using the BRIM algorithm to compute the result
brim.mut.info<-function(link.data,n.net){
  norm.mut.info.brim.vec<-NULL
  for(i in 1:n.net){
    #It needs the input matrix to content no all-zero row/column, need removal of all-zero column.
    #What if there are all-zero rows? No, all-zero row hardly exists.
    zero.ind<-!apply(link.data[[i]],2,"sum")==0
    input<-link.data[[i]][,zero.ind]
    #BRIM
    res.brim<-findModules(input,spar=FALSE)
    #extract S-matrix list for GPs only
    s.mat.brim.gp<-res.brim$S[1:ngp,]
    #Compute the normalized mutual information for all n.net networks
    norm.mut.info.brim.vec<-c(norm.mut.info.brim.vec,norm.mut.infor(s.mat.gp.true,s.mat.brim.gp))
  }
  norm.mut.info.brim.vec
}

#function #4 using the walktrap algorithm on the projection of the bipartite graph to get the result
proj.mut.info<-function(link.data,ind.vec,n.net){
  norm.mut.info.proj.vec<-NULL
  ngp=sum(ind.vec)
  for(i in 1:n.net){
    #make adjacency matrix
    link.data.adj<-rbind(cbind(matrix(rep(0,ngp^2),ncol=ngp),link.data[[i]]),cbind(t(link.data[[i]]),matrix(rep(0,npa^2),ncol=npa)))
    #produce the bipartite graph object
    temp.bipa.graph<-graph_from_adjacency_matrix(link.data.adj)
    V(temp.bipa.graph)$name<-c(1:(npa+ngp))
    V(temp.bipa.graph)$type<-ind.vec==0
    temp.proj.graph<-bipartite_projection(temp.bipa.graph)$proj1
    wc<-membership(cluster_walktrap(temp.proj.graph))
    s.mat.proj.gp<-matrix(0,ncol=max(wc),nrow=ngp)
    for(j in 1:ngp)s.mat.proj.gp[j,wc[j]]<-1
    norm.mut.info.proj.vec<-c(norm.mut.info.proj.vec,norm.mut.infor(s.mat.gp.true,s.mat.proj.gp))
  }
  norm.mut.info.proj.vec
}



#control parameters
#setseed
set.seed(2389)
#number of GP
ngp<-20
#number of patients
npa<-50
#number of communitys #TRUE#
ncomm.true<-5
#Probability of the link construction
p.in=0.6
p.out=0.3
#numbers of networks generated
n.net=1000
#debugging setting
#n.net=3

#simulation before apply algorithms
#specify the #TRUE# index vector&matrix
comm.ind.true<-ceiling(runif((ngp+npa),0,ncomm.true))
# true S-martix for GPs only
s.mat.gp.true<-matrix(0,ncol=ncomm.true,nrow=ngp)
for(i in 1:ngp){
  s.mat.gp.true[i,comm.ind.true[i]]<-1
}
#index of whether i-th vertex is gp or not
ind.vec<-c(rep(1,ngp),rep(0,npa))
#simulate n.net network links adjacency matrix
link.data<-link.create(ind.vec,comm.ind.true,p.in,p.out,n.net)




#Bipartite-BRIM
norm.mut.info.brim.vec<-brim.mut.info(link.data,n.net)



#Projected Unipartite-walktrap
norm.mut.info.proj.vec<-proj.mut.info(link.data,ind.vec,n.net)

