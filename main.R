#G3tYourMa5t3r!!! Yaarctopus

#library statements
require(igraph)
require(lpbrim)
require(bipartite)
require(Matrix)

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
brim.mut.info<-function(link.data,ind.vec,n.net,s.mat.gp.true){
  norm.mut.info.brim.vec<-NULL
  ngp<-sum(ind.vec)
  for(i in 1:n.net){
    #It needs the input matrix to content no all-zero row/column, need removal of all-zero columns&rows.
    #zero-columns, easier to be treated
    nozero.ind.col<-!apply(link.data[[i]],2,"sum")==0
    input<-link.data[[i]][,nozero.ind.col]
    #zero-rows, a bit harder
    nozero.ind.row<-!apply(input,1,"sum")==0
    post.treat=FALSE
    if(sum(nozero.ind.row)!=ngp){
      post.treat=TRUE
      input<-input[nozero.ind.row,]
    }
    #BRIM
    res.brim<-findModules(input,spar=FALSE)
    #extract S-matrix list for GPs only
    if(post.treat){
      s.mat.brim.gp<-matrix(0,nrow=ngp,ncol=ncol(res.brim$S))
      s.mat.brim.gp[nozero.ind.row,]<-res.brim$S[1:sum(nozero.ind.row),]
      addition<-matrix(0,nrow=ngp,ncol=sum(!nozero.ind.row))
      addition[!nozero.ind.row,]<-diag(1,sum(!nozero.ind.row),sum(!nozero.ind.row))
      s.mat.brim.gp<-cbind(s.mat.brim.gp,addition)      
    }else{
      s.mat.brim.gp<-res.brim$S[1:ngp,]
    }
    #Compute the normalized mutual information for all n.net networks
    norm.mut.info.brim.vec<-c(norm.mut.info.brim.vec,norm.mut.infor(s.mat.gp.true,s.mat.brim.gp))
  }
  norm.mut.info.brim.vec
}

#function #4 using the optimal? algorithm on the projection of the bipartite graph to get the result
proj.mut.info<-function(link.data,ind.vec,n.net,s.mat.gp.true){
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
    wc<-membership(cluster_optimal(temp.proj.graph))
    s.mat.proj.gp<-matrix(0,ncol=max(wc),nrow=ngp)
    for(j in 1:ngp)s.mat.proj.gp[j,wc[j]]<-1
    norm.mut.info.proj.vec<-c(norm.mut.info.proj.vec,norm.mut.infor(s.mat.gp.true,s.mat.proj.gp))
  }
  norm.mut.info.proj.vec
}


#Part 1 #for a set p.in/p.out
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
#simulate n.net network links matrix
link.data<-link.create(ind.vec,comm.ind.true,p.in,p.out,n.net)


#Bipartite-BRIM
norm.mut.info.brim.vec<-brim.mut.info(link.data,ind.vec,n.net,s.mat.gp.true)


#Projected Unipartite-optimal?
norm.mut.info.proj.vec<-proj.mut.info(link.data,ind.vec,n.net,s.mat.gp.true)


#Part 2 #for different set p.in/p.out
#I write it as a function
#a huge one

#control parameters
#setseed
set.seed(2389)

#Probability of the link construction
# Say p.in have a upper bound 1, lower bound 0.5
#p.out is always less than p.in, lower bound 0.1
#takeing step=0.01
p.in.upper=1
p.in.lower=0.5
p.out.lower=0.1
step=0.01

#number of GP
ngp<-20
#number of patients
npa<-50
#number of communitys #TRUE#
ncomm.true<-5

#numbers of networks generated
n.net=1000
#debugging setting
#n.net=1



#function #5 crazy function study for the possible influences for the p.in/p.out setting
crazy.function<-function(p.in.lower,p.in.upper,step,p.out.lower,n.net,ngp,npa,ncomm.true){
  result.list<-NULL
  p.in.vec<-seq(p.in.lower,p.in.upper,step)
  for(p.in in p.in.vec){
    p.out.vec<-seq(p.out.lower,p.in,step)
    temp.result.list.sub<-NULL
    for(p.out in p.out.vec){
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
      #simulate n.net network links matrix
      link.data<-link.create(ind.vec,comm.ind.true,p.in,p.out,n.net)
      
      #Bipartite-BRIM
      norm.mut.info.brim.vec<-brim.mut.info(link.data,ind.vec,n.net,s.mat.gp.true)
            
      #Projected Unipartite-optimal?
      norm.mut.info.proj.vec<-proj.mut.info(link.data,ind.vec,n.net,s.mat.gp.true)
      
      temp.result.list<-list(p.in=p.in,p.out=p.out,BRIM=norm.mut.info.brim.vec,Projection=norm.mut.info.proj.vec)
      
      temp.result.list.sub[[length(temp.result.list.sub)+1]]<-temp.result.list
    }
    result.list[[length(result.list)+1]]<-temp.result.list.sub
  }
  result.list
}



cross.fingers<-crazy.function(p.in.lower,p.in.upper,step,p.out.lower,n.net,ngp,npa,ncomm.true)






















#real GP-patients dataset
#the scrambled_ppn have a initial shift of 99(1->100), therefore create a new raw.csv dataset file.
#Convert the raw data into a huge matrix and store it as "Medicare_link_data"
#raw<-read.csv("raw.csv",header=TRUE)
#medicare.link.data<-Matrix(0,nrow=max(raw$scrambled_prov),ncol=max(raw$scrambled_ppn))
#for(i in 1:nrow(raw)){
#  temp<-raw[i,]
#  medicare.link.data[temp$scrambled_prov,temp$scrambled_ppn]<-medicare.link.data[temp$scrambled_prov,temp$scrambled_ppn]+1
#}
#writeMM(medcare.link.data,file="Medicare_link_data")
medicare.link.data<-readMM("Medicare_link_data")
