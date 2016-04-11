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




#control parameters
#setseed
set.seed(1132)
#number of GP
ngp<-20
#number of patients
npa<-50
#number of communitys #TRUE#
ncomm.true<-5
#Probability of the link construction
p.in=0.8
p.out=0.5
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


#Bipartite-BRIM#oping
res.brim.list<-NULL
s.mat.brim.gp.list<-NULL
norm.mut.info.brim.vec<-NULL
for(i in 1:n.net){
  #BRIM and save the result
  res.brim.list[[length(res.brim.list)+1]]<-findModules(link.data[[i]],spar=FALSE)
  #extract S-matrix list for GPs only
  s.mat.brim.gp.list[[length(s.mat.brim.gp.list)+1]]<-res.brim.list[[i]]$S[1:ngp,]
  #Compute the normalized mutual information for all n.net networks
  norm.mut.info.brim.vec<-c(norm.mut.info.brim.vec,norm.mut.infor(s.mat.gp.true,s.mat.brim.gp.list[[i]]))
}


#Projected Unipartite-walktrap
#create adjacency matrix for each network
link.data.adj<-NULL
bipa.graph.list<-NULL
for(i in 1:n.net){
  #make adjacency matrix
  link.data.adj[[length(link.data.adj)+1]]=rbind(cbind(matrix(rep(0,ngp^2),ncol=ngp),link.data[[i]]),cbind(t(link.data[[i]]),matrix(rep(0,npa^2),ncol=npa)))
  #produce the bipartite graph object
  temp<-graph_from_adjacency_matrix(link.data.adj[[i]])
  V(temp)$name<-c(1:(npa+ngp))
  V(temp)$type<-bipartite_mapping(temp)$type
  bipa.graph.list[[length(bipa.graph.list)+1]]<-temp
  
}

