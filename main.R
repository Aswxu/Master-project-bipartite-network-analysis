#G3tYourMa5t3r!!! Yaarctopus

#uncomment to install packages
#install.packages("igraph", dependencies = TRUE)
#install.packages("lpbrim", dependencies = TRUE)
#install.packages("bipartite", dependencies = TRUE)
#install.packages("Matrix", dependencies = TRUE)
#install.packages("ggplot2", dependencies = TRUE)
#install.packages("wesanderson")

#library statements
require(igraph)
require(lpbrim)
require(bipartite)
require(Matrix)
require(ggplot2)
require(boot)
require(wesanderson)

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

#function #4 using the cluster_louvain algorithm on the projection of the bipartite graph to get the result
proj.mut.info.louvain<-function(link.data,ind.vec,n.net,s.mat.gp.true){
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
    lv<-membership(cluster_louvain(temp.proj.graph))
    s.mat.proj.gp<-matrix(0,ncol=max(lv),nrow=ngp)
    for(j in 1:ngp)s.mat.proj.gp[j,lv[j]]<-1
    norm.mut.info.proj.vec<-c(norm.mut.info.proj.vec,norm.mut.infor(s.mat.gp.true,s.mat.proj.gp))
  }
  norm.mut.info.proj.vec
}


#function #5 using the cluster_walktrap algorithm on the projection of the bipartite graph to get the result
proj.mut.info.walktrap<-function(link.data,ind.vec,n.net,s.mat.gp.true){
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
    
    #wt <- cluster_walktrap(temp.proj.graph, modularity=TRUE)
    #wc <- membership(cluster_walktrap(temp.proj.graph, 
    #                    steps=which.max(wt$modularity)-1))
    
    #This might not be right
    wc<-membership(cluster_walktrap(temp.proj.graph,steps=4))
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
n.net=10000
#debugging setting
#n.net=1000

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

#Projected Unipartite-louvain
norm.mut.info.louvain.vec<-proj.mut.info.louvain(link.data,ind.vec,n.net,s.mat.gp.true)

#Projected Unipartite-walktrap
norm.mut.info.walktrap.vec<-proj.mut.info.walktrap(link.data,ind.vec,n.net,s.mat.gp.true)


#part 2 #ggplot
#create the data.frame
#random assignment
random.assign<-NULL
for(j in 1:10000){
  random.ncomm<-ceiling(runif(1,0,ngp))
  comm.ind.temp<-as.matrix(ceiling(runif(ngp,0,ncomm.true)))
  s.mat.gp.temp<-matrix(0,ncol=ncomm.true,nrow=ngp)
  for(i in 1:ngp){
    s.mat.gp.temp[i,comm.ind.temp[i]]<-1
  }
  random.assign<-c(random.assign,norm.mut.infor(s.mat.gp.true,s.mat.gp.temp))
}
nmi.dataframe<-data.frame(nmi=c(norm.mut.info.walktrap.vec,norm.mut.info.louvain.vec,norm.mut.info.brim.vec,random.assign),
                          Method=c(rep("Walktrap",10000),rep("Louvain",10000),rep("LP-BRIM",10000),rep("Random",10000)))
graph3method <- ggplot(nmi.dataframe, aes(nmi))
graph3method + geom_histogram(aes(fill=Method,y=..density..),alpha=0.7,size=0.2,colour='black',position='dodge',bins=15)+ 
  stat_density(geom='line',position='identity',size=0.7, aes(colour=Method))+
  #scale_fill_manual(values=wes_palette(n=4, name="FantasticFox"))+
  #scale_color_manual(values=wes_palette(n=4, name="FantasticFox"))+
  scale_fill_brewer(palette="Set1")+
  scale_color_brewer(palette="Set1")+
  scale_y_continuous("Density")+
  scale_x_continuous("Normalized Mutual Information") 


#Part 3 #for different set p.in/p.out
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
step=0.02

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



#function #6 crazy function study for the possible influences for the p.in/p.out setting
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
      
      #Projected Unipartite-walktrap
      norm.mut.info.wt.vec<-proj.mut.info.walktrap(link.data,ind.vec,n.net,s.mat.gp.true)
      
      norm.mut.info.lv.vec<-proj.mut.info.louvain(link.data,ind.vec,n.net,s.mat.gp.true)
      
      temp.result.list<-list(p.in=p.in,p.out=p.out,BRIM=norm.mut.info.brim.vec,Walktrap=norm.mut.info.wt.vec,Louvain=norm.mut.info.lv.vec)
      
      temp.result.list.sub[[length(temp.result.list.sub)+1]]<-temp.result.list
    }
    result.list[[length(result.list)+1]]<-temp.result.list.sub
  }
  result.list
}


#it is working now!
cross.fingers<-crazy.function(p.in.lower,p.in.upper,step,p.out.lower,n.net,ngp,npa,ncomm.true)








#part 4 #ploting stuffs
#projection
adj.mat<-matrix(c(rep(0,3),1,rep(0,3),1,rep(0,3),1,rep(1,3),0),nrow=4,ncol=4)
temp.bipa.graph<-graph_from_adjacency_matrix(adj.mat,mode="undirected")
V(temp.bipa.graph)$name<-c(1:3,"A")
V(temp.bipa.graph)$type<-c(rep(1,3),0)
lay.bi<-layout.bipartite(temp.bipa.graph)
plot(temp.bipa.graph, layout = lay.bi[, c(1,2)],vertex.size=30,vertex.color="darkblue",vertex.label.color="yellow",vertex.label.font=2,edge.color="black")
temp.proj.graph<-bipartite_projection(temp.bipa.graph,multiplicity = TRUE,which="TRUE")
plot(temp.proj.graph,vertex.size=30,vertex.color="darkblue",vertex.label.color="yellow",vertex.label.font=2,edge.color="black")
adj.mat<-matrix(c(rep(0,3),1,1,0,rep(0,3),1,0,1,rep(0,3),0,1,1,1,1,0,rep(0,3),1,0,1,rep(0,3),0,1,1,rep(0,3)),nrow=6,ncol=6)
temp.bipa.graph<-graph_from_adjacency_matrix(adj.mat,mode="undirected")
V(temp.bipa.graph)$name<-c(1:3,"A","B","C")
V(temp.bipa.graph)$type<-c(rep(1,3),0,0,0)
lay.bi<-layout.bipartite(temp.bipa.graph)
plot(temp.bipa.graph, layout = lay.bi[, c(1,2)],vertex.size=30,vertex.color="darkblue",vertex.label.color="yellow",vertex.label.font=2,edge.color="black")







#part 5# weird peak in walktrap of nmi==0.6869212 which is actually the maxium
#norm.mut.info.walktrap.vec[norm.mut.info.walktrap.vec>0.66]==0.6869212
#weird.ind<-norm.mut.info.walktrap.vec>0.66
#weird.network.samples<-link.data[weird.ind]
#  #make adjacency matrix
#  link.data.adj<-rbind(cbind(matrix(rep(0,ngp^2),ncol=ngp),weird.network.samples[[1]]),cbind(t(weird.network.samples[[1]]),matrix(rep(0,npa^2),ncol=npa)))
#  #produce the bipartite graph object
#  temp.bipa.graph<-graph_from_adjacency_matrix(link.data.adj)
#  V(temp.bipa.graph)$name<-c(1:(npa+ngp))
#  V(temp.bipa.graph)$type<-ind.vec==0
#  temp.proj.graph<-bipartite_projection(temp.bipa.graph)$proj1
#  
#  wc<-walktrap.community(temp.proj.graph,steps=5)
#  plot(wc, temp.proj.graph, mark.col=heat.colors(length(wc)),
#       col=rainbow(length(wc))[membership(wc)])
#the conclusion is that each vertex are assigned into separete communities





#part 6 #real GP-patients dataset
#the scrambled_ppn have a initial shift of 99(1->100), therefore create a new raw.csv dataset file.
raw<-read.csv("raw.csv",header=TRUE)
#max(table(raw$scrambled_pcode))
##[1] 11692
#c(1000:1541)[table(raw$scrambled_pcode)==11692]
##[1] 1446
#mean(table(raw$scrambled_pcode))
##[1] 938.9207
#c(1000:1541)[table(raw$scrambled_pcode)==938]
##[1] 1153

#function #6 for rearrange the number
reflect<-function(ref,data){
  neo.ref<-c(1:length(ref))
  neo<-rep(0,length(data))
  for(i in 1:length(ref)){
    temp<-data==ref[i]
    neo[temp]<-neo.ref[i]
  }
  neo
}


#pick pcode=1446 & pcode=1153
raw.1446<-subset(raw,raw$scrambled_pcode==1446)
raw.1153<-subset(raw,raw$scrambled_pcode==1153)

#rearrange the network tag
raw.1446$scrambled_prov<-reflect(as.numeric(names(table(raw.1446$scrambled_prov))),raw.1446$scrambled_prov)
raw.1446$scrambled_ppn<-reflect(as.numeric(names(table(raw.1446$scrambled_ppn))),raw.1446$scrambled_ppn)
raw.1446$scrambled_ppn<-raw.1446$scrambled_ppn+max(raw.1446$scrambled_prov)

raw.1153$scrambled_prov<-reflect(as.numeric(names(table(raw.1153$scrambled_prov))),raw.1153$scrambled_prov)
raw.1153$scrambled_ppn<-reflect(as.numeric(names(table(raw.1153$scrambled_ppn))),raw.1153$scrambled_ppn)
raw.1153$scrambled_ppn<-raw.1153$scrambled_ppn+max(raw.1153$scrambled_prov)

#turn into adj'mat
graph.1446<-graph_from_edgelist(as.matrix(raw.1446[,c(1,3)]),directed=FALSE)
graph.1153<-graph_from_edgelist(as.matrix(raw.1153[,c(1,3)]),directed=FALSE)
ind.1446<-c(rep(1,max(raw.1446$scrambled_prov)),
            rep(0,max(raw.1446$scrambled_ppn)-max(raw.1446$scrambled_prov)))
ind.1153<-c(rep(1,max(raw.1153$scrambled_prov)),
            rep(0,max(raw.1153$scrambled_ppn)-max(raw.1153$scrambled_prov)))
adj.mat.1446<-as_adj(graph.1446)
adj.mat.1153<-as_adj(graph.1153)
adj.prime.mat.1446<-adj.mat.1446[c(1:sum(ind.1446)),c((sum(ind.1446)+1):length(ind.1446))]
adj.prime.mat.1153<-adj.mat.1153[c(1:sum(ind.1153)),c((sum(ind.1153)+1):length(ind.1153))]


#pcode 1446
V(graph.1446)$name<-c(1:length(ind.1446))
V(graph.1446)$type<-ind.1446==0
proj.graph.1446<-bipartite_projection(graph.1446)$proj1
lv.1446<-cluster_louvain(proj.graph.1446)
wc.1446<-cluster_walktrap(proj.graph.1446)
res.brim.1446<-findModules(as.matrix(adj.prime.mat.1446),spar=FALSE)
res.1446<-res.brim.1446$S[c(1:sum(ind.1446)),]
brim.1446<-res.1446%*%as.numeric(row.names(t(res.1446)))
write.csv(data.frame(Label=c(1:sum(ind.1446)),
                     lv.membership=lv.1446$membership,
                     wc.membership=wc.1446$membership,
                     brim.membership=brim.1446),
          file="nodename_pcode_1446_proj.csv")
edgelist.1446.proj<-as_edgelist(proj.graph.1446)
write.csv(data.frame(Source=edgelist.1446.proj[,1]-1,
                     Target=edgelist.1446.proj[,2]-1,
                     Type="Undirected",Weight=E(proj.graph.1446)$weight),
          file="edge_pcode_1446_proj.csv")


#pcode 1153
V(graph.1153)$name<-c(1:length(ind.1153))
V(graph.1153)$type<-ind.1153==0
proj.graph.1153<-bipartite_projection(graph.1153)$proj1
lv.1153<-cluster_louvain(proj.graph.1153)
wc.1153<-cluster_walktrap(proj.graph.1153)
res.brim.1153<-findModules(as.matrix(adj.prime.mat.1153),spar=FALSE)
res.1153<-res.brim.1153$S[c(1:sum(ind.1153)),]
brim.1153<-res.1153%*%as.numeric(row.names(t(res.1153)))
write.csv(data.frame(Label=c(1:sum(ind.1153)),
                     lv.membership=lv.1153$membership,
                     wc.membership=wc.1153$membership,
                     brim.membership=brim.1153),
          file="nodename_pcode_1153_proj.csv")
edgelist.1153.proj<-as_edgelist(proj.graph.1153)
write.csv(data.frame(Source=edgelist.1153.proj[,1]-1,
                     Target=edgelist.1153.proj[,2]-1,
                     Type="Undirected",
                     Weight=E(proj.graph.1153)$weight),
          file="edge_pcode_1153_proj.csv")



#part 7 # look deeper for the sample which have highest differeces between lpbrim and louvain
#max(norm.mut.info.brim.vec)
#[1] 0.8101614
#c(1:10000)[norm.mut.info.brim.vec==max(norm.mut.info.brim.vec)]
#[1] 2299 7952
#make adjacency matrix
dl.ind=2292
deeplook<-rbind(cbind(matrix(rep(0,ngp^2),ncol=ngp),link.data[[dl.ind]]),cbind(t(link.data[[dl.ind]]),matrix(rep(0,npa^2),ncol=npa)))
#produce the bipartite graph object
deeplook.bipa.graph<-graph_from_adjacency_matrix(deeplook)
V(deeplook.bipa.graph)$name<-c(1:(npa+ngp))
V(deeplook.bipa.graph)$type<-ind.vec==0

deeplook.proj.graph<-bipartite_projection(deeplook.bipa.graph)$proj1
deeplook.lv<-membership(cluster_louvain(deeplook.proj.graph))
deeplook.lv<-c(deeplook.lv+1,rep(1,npa))

res.brim.deeplook<-findModules(link.data[[dl.ind]],spar=FALSE)
res.deeplook<-res.brim.deeplook$S[c(1:sum(ind.vec)),]
deeplook.brim<-res.deeplook%*%as.numeric(row.names(t(res.deeplook)))
deeplook.brim<-c(deeplook.brim+1,rep(1,npa))

deeplook.true<-comm.ind.true[c(1:20)]
deeplook.true<-c(deeplook.true+1,rep(1,npa))

write.csv(data.frame(Label=c(1:70),
                     deeplook.lv=deeplook.lv,
                     deeplook.brim=deeplook.brim,
                     deeplook.true=deeplook.true,
                     isgp=ind.vec),
          file="nodename_deeplook.csv")
edgelist.deeplook<-as_edgelist(deeplook.bipa.graph)
write.csv(data.frame(Source=edgelist.deeplook[,1]-1,
                     Target=edgelist.deeplook[,2]-1,
                     Type="Directed"),
          file="edge__deeplook.csv")
#not making sense


#part 8 # t-test and bootstrap
t.fun<-function(d,i){
  boot.sample<-d[i]
  m<-mean(boot.sample)
  s<-sd(boot.sample)
  t<-(m)/s*sqrt(length(boot.sample))
  t
}
#A
sa<-norm.mut.info.brim.vec-norm.mut.info.walktrap.vec
#t-statistic function for bootstrap
find.t<-function(t){
  a<-exp(t*sa)
  p<-a/sum(a)
  out<-sum(p*sa)
  out
}
t.root<-uniroot(find.t,c(-7,7))$root
temp<-exp(t.root*sa)
p.a<-temp/sum(temp)
B=5000
boot.obj.a<-boot(sa,t.fun,B,stype="i",weights=p.a)
p.val.a<-sum(boot.obj.a$t>=boot.obj.a$t0)/B
#B
sb<-norm.mut.info.brim.vec-norm.mut.info.louvain.vec
find.t<-function(t){
  a<-exp(t*sb)
  p<-a/sum(a)
  out<-sum(p*sb)
  out
}
t.root<-uniroot(find.t,c(-7,7))$root
temp<-exp(t.root*sb)
p.b<-temp/sum(temp)
B=5000
boot.obj.b<-boot(sb,t.fun,B,stype="i",weights=p.b)
p.val.b<-sum(boot.obj.b$t>=boot.obj.b$t0)/B

#describe expect lpbrim
mean(norm.mut.info.brim.vec)
var(norm.mut.info.brim.vec)
w=sort(norm.mut.info.brim.vec)
w[0.025*10000+1]
w[0.975*10000]