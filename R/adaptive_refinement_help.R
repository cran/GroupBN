#Assisting functions in alphabetical order

#calculate weighted average cross entropy
#pred: predictions, obs: observations (0,1)
#weights are the group proportions
#man kann beide if abfragen zusammenfassen indem man nur den zweiten teil nimmt?
cross.en<-function(pred, obs, sdpred=NULL,weighted=T){
  if(nlevels(factor(obs))<2|nlevels(factor(obs))>4){
    stop('implemented only for a categorical target variable with 2-4 levels')
  }else if(nlevels(factor(obs))>2){
    k<-nlevels(factor(obs))
    pred[pred==0]<-0.000001
    pred[pred==1]<-0.999999
    weight<-pred[!is.na(pred[,1]),1]
    if(weighted==T){
      for (i in 1:length(pred[!is.na(pred[,1]),1])){
        weight[i]<-length(which(obs[!is.na(pred[,1])]!=obs[!is.na(pred[,1])][i]))/length(obs[!is.na(pred[,1])])
      }
    }else{
      weight<-rep(1, length(weight))
    }
    loss<-0
    for (i in 1:length(pred[!is.na(pred[,1]),1])){
      for (j in 1:k){
        if(obs[i]==j&!is.na(pred[i,1])){
          loss<-loss-weight[i]*log2(pred[i,j])
        }
      }
    }
  }else{
    pred[pred==0]<-0.000001
    pred[pred==1]<-0.999999
    if(weighted==T){
      weight<-pred[!is.na(pred)]
      weight[as.numeric(obs[!is.na(pred)])==2]<-length(which(as.numeric(obs[!is.na(pred)])==1))/length(obs[!is.na(pred)])
      weight[as.numeric(obs[!is.na(pred)])==1]<-length(which(as.numeric(obs[!is.na(pred)])==2))/length(obs[!is.na(pred)])
    } else {
      weight<-rep(1, length(pred[!is.na(pred)]))#/length(pred[!is.na(pred)])
    }
    if(!is.null(sdpred)){
      loss<-(-sum(weight*(as.numeric(obs[!is.na(pred)])-1)*log2(1-pred[!is.na(pred)]+sdpred[!is.na(pred)]))-sum(weight*(2-as.numeric(obs[!is.na(pred)]))*log2(pred[!is.na(pred)]+sdpred[!is.na(pred)])))
    } else {
      loss<-(-sum(weight*(as.numeric(obs[!is.na(pred)])-1)*log2(1-pred[!is.na(pred)]))-sum(weight*(2-as.numeric(obs[!is.na(pred)]))*log2(pred[!is.na(pred)])))
    }
  }
  return(loss)
}

#discretize and factorize scores of a groupbn object
disc.scores<-function(res, seed=NULL){
  k<-max(res$grouping)
  param<-res$disc.param
  data.scores<-res$group.data
  colnames(data.scores)<-paste("cl", 1:k, sep="")
  data.scores<-as.data.frame(data.scores)
  data.scores.dist<-data.scores
  for (i in 1:dim(data.scores.dist)[2]){
    title<-colnames(data.scores)[i]
    #print(title)
    disc<-discretize.dens(data.scores.dist[,i], cluster=T, title=title, seed=seed)
    data.scores.dist[,i]<-as.factor(disc$discretized)
    param[[i]]<-disc$levels

  }
  res$group.data<-data.scores.dist
  res$disc.param<-param
  return(res)
}

#density approximative discretization
#graph: Boolean, if TRUE plots are produced
#rename.level: Boolean, if TRUE levels are renamed to integers (1,2,3,...)
#cluster: Boolean, if TRUE it is checked if a cluster is already a discrete variable
discretize.dens<-function(data, graph=F, title="Density-approxmative Discretization", rename.level=F, return.all=T, cluster=F, seed=NULL){
  if (!is.null(seed)) {set.seed(seed)}
  if (is.factor(data)){
    if(graph){
      plot<-plot(data,
            col="lightblue4", # column color
            border="black",
            xlab="Value",
            ylab="Frequency",
            main = title)
    }
    message("Data are already discrete.")
    result<-list(plot, discretized=data, levels=NULL)
    if(return.all){
      invisible(result)
    } else {
      invisible(data)
    }
  }
  #check if a cluster already has a limited number of states (<10)
  if (cluster==T&length(table(as.factor(data))[table(as.factor(data))>1])<10&length(table(as.factor(data))[table(as.factor(data))>1])>1){
    m<-nlevels(as.factor(data))
    #all levels exist more than once in the data (then keep all levels)
    if(all(table(as.factor(data))>1)){
      d<-as.factor(data)
      lev<-attr(d, "discretized:breaks")
      if(rename.level){
        levels(d)<-1:m
      }
      if(graph){
        plot<-plot(stats::density(data, na.rm=T), col="black", lwd=3,main = paste("Density Approxmative Discretization: ", title, sep=""))
        plot<-graphics::hist(data, # histogram
                   col="lightblue4", # column color
                   border="black",
                   breaks=lev,
                   xlab="Value",
                   prob = TRUE, # show densities instead of frequencies
                   main = title,
                   add=TRUE)
        plot<-graphics::lines(stats::density(data, na.rm=T), col="black", lwd=3, main = paste("Density Approxmative Discretization: ", title, sep=""))}
      #plot<-abline(v=lev, col="lightblue4", lty="dashed", lwd=3)
      #plot<-hist(data, col="red", breaks=lev)

      #set outer boundaries to Inf
      lev[1]<--Inf
      lev[length(lev)]<-Inf
      result<-list(plot, discretized=d, levels=lev)

      if(return.all){
        invisible(result)
      } else {
        invisible(d)
      }
    }
    #otherwise use cluster-based discretization policy
    else{
      cent<-as.numeric(names(table(as.factor(data))[table(as.factor(data))>1]))
      d<-arules::discretize(data, method="cluster", centers=cent)
      lev<-attr(d, "discretized:breaks")
      if(rename.level){
        levels(d)<-1:length(cent)
      }
      result<-list(discretized=d)
      xy=0
      max<-length(cent)
      while(min(table(d))<20&xy<max-2){
        xy=xy+1
        cent<-as.numeric(names(sort(table(as.factor(data))[table(as.factor(data))>1], decreasing=T)[1:(max-xy)]))
        d<-arules::discretize(data, method="cluster", centers=cent)
        lev<-attr(d, "discretized:breaks")
        if(rename.level){
          levels(d)<-1:length(cent)
        }
      }
      if(graph){
        plot<-plot(stats::density(data, na.rm=T), col="black", lwd=3,main = paste("Density Approxmative Discretization: ", title, sep=""))
        plot<-graphics::hist(data, # histogram
                   col="lightblue4", # column color
                   border="black",
                   breaks=lev,
                   xlab="Value",
                   prob = TRUE, # show densities instead of frequencies
                   main = title,
                   add=TRUE)
        plot<-graphics::lines(stats::density(data, na.rm=T), col="black", lwd=3, main = paste("Density Approxmative Discretization: ", title, sep=""))}
      #plot<-abline(v=lev, col="lightblue4", lty="dashed", lwd=3)
      #plot<-hist(data, col="red", breaks=lev)

      #set outer boundaries to Inf
      lev[1]<--Inf
      lev[length(lev)]<-Inf

      result<-list(plot, discretized=d, levels=lev)

      if(return.all){
        invisible(result)
      } else {
        invisible(d)
      }
    }
  }
  else{
    #get local maxima of the density
    max<-get.maxima(data)
    max.nr<-dim(max)[2]
    #if more than one local maximum: discretize around cluster means
    if (max.nr>1){
      int<-max.nr
      dist.max<-max(max[1,])-min(max[1,])
      range<-abs(stats::quantile(data, prob=0.05, na.rm=T)-stats::quantile(data, prob=0.95, na.rm=T))
      if (dist.max<range/5){
        int<-int+1
        d<-arules::discretize(data, method="cluster", centers=c(max[1,],mean(data, na.rm=T)+0.3))
      }
      else{
        d<-arules::discretize(data, method="cluster", centers=max[1,])
      }
      lev<-attr(d, "discretized:breaks")
      if(rename.level){
        levels(d)<-1:int
      }
    }
    #if only one local maximum and quartiles are unique: discretize to quartiles
    else if (max.nr==1){
      #unique quartiles?
      if (length(unique(stats::quantile(data, probs=c(0,0.25,0.5,0.75,1), na.rm=TRUE))) == 5){
        #rem<-which.min(abs(quantile(data, na.rm = T)[2:4]-max[1,]))
        #if (rem==2){
        int=4
        br<-stats::quantile(data, na.rm=T)
        d<-arules::discretize(data, method="fixed", breaks=br)
        lev<-attr(d, "discretized:breaks")
        if(rename.level){
          levels(d)<-1:int
        }
        #}
        #else if (rem==1){
        #  int=3
        #  br<-quantile(data, na.rm=T)[c(1,3,4,5)]
        #  d<-arules::discretize(data, method="fixed", breaks=br)
        #  lev<-attr(d, "discretized:breaks")
        #   levels(d)<-1:int
        #}
        # else {
        #  int=3
        #  br<-quantile(data, na.rm=T)[c(1,2,3,5)]
        #  d<-arules::discretize(data, method="fixed", breaks=br)
        # lev<-attr(d, "discretized:breaks")
        #  levels(d)<-1:int
        #}
      }
      else {
        x=0.05
        lev<-unique(stats::quantile(data, prob=c(0,x,1-x,1), na.rm = T))
        while(length(lev)<=2&x>=0){
          x=x-0.01
          lev<-unique(stats::quantile(data, prob=c(0,x,1-x,1), na.rm = T))
        }
        d<-arules::discretize(data, method="fixed", breaks=lev)
        if(rename.level){
          levels(d)<-1:length(unique(stats::quantile(data, na.rm=TRUE, prob=c(0,0.05,0.95,1))))
        }
        #int<-3
        #d<-arules::discretize(data, method="cluster", breaks=int)
        #lev<-attr(d, "discretized:breaks")
        #levels(d)<-1:int
      }
    }

    #plot histograms with intervals to check
    if(graph){
      plot<-plot(stats::density(data, na.rm=T), col="black", lwd=3,main = paste("Density Approxmative Discretization: ", title, sep=""))
      plot<-graphics::hist(data, # histogram
                 col="lightblue4", # column color
                 border="black",
                 breaks=lev,
                 xlab="Value",
                 prob = TRUE, # show densities instead of frequencies
                 main = title,
                 add=TRUE)
      plot<-graphics::lines(stats::density(data, na.rm=T), col="black", lwd=3, main = paste("Density Approxmative Discretization: ", title, sep=""))
      plot<-graphics::lines(max[1,], max[2,], type="p",pch=8, col="lightblue")}
    #plot<-abline(v=lev, col="lightblue4", lty="dashed", lwd=3)
    #plot<-hist(data, col="red", breaks=lev)
    #set outer boundaries to Inf
    lev[1]<--Inf
    lev[length(lev)]<-Inf
    result<-list(plot, discretized=d, levels=lev, optima=max)
    if(return.all){
      invisible(result)
    } else {
      invisible(d)
    }
  }
}

#get all variables belonging to specific cluster from merge matrix
get.cluster<-function(cl, ME){
  idx<-cl
  while(any(idx>0)){
    new<-c()
    for (i in 1:length(idx)){
      if (idx[i]>0){
        new<-c(new, ME[idx[i],])
      }
    }
    idx<-c(idx[which(idx<0)], new)
  }
  return(-1*idx)
}

#find local optima in data (used for discretization function)
get.maxima<-function(data1){
  dens<-stats::density(data1, na.rm = T)
  vec<-dens$y
  #ordered observations
  vec.zoo<- zoo::as.zoo(vec)
  #wrap zeros around to take start and end also into account
  vec.zoo<-c(rep(0,10), vec.zoo, rep(0,10))
  #find maximum within a range of 21 data points
  top <- zoo::rollapply(vec.zoo, 21, function(x) which.max(x)==11)
  #get indices
  top.ind<-zoo::index(top)[zoo::coredata(top)]

  #only keep if peak higher than 0.01
  #if (length(which(dens$y[top.ind]>mean(vec)))>1){
  #  top.ind<-top.ind[which(dens$y[top.ind]>mean(vec))]
  # }
  #else{
  top.y<-dens$y[top.ind]
  sub<-sort(top.y, decreasing = TRUE)[which(sort(top.y, decreasing = T)>0.05*sort(top.y, decreasing = T)[1])]
  if(length(sub)>4){
    sub<-sort(top.y, decreasing = T)[1:4]
  }
  top.ind<-top.ind[which(top.ind%in%which(dens$y%in%sub))]
  #}
  top.x<-dens$x[top.ind]
  top.y<-dens$y[top.ind]

  return(rbind(top.x, top.y))
}

#get row/step in which a variable or cluster is merged
#returns vector of merging step (row of merge matrix) and column of the merged cluster/variable
get.merge.step<-function(j.var=0, j.cl=0, ME){
  if(j.var>0){
    idx<-which(ME[,]==-j.var)
    if (idx>=dim(ME)[1]+1){
      idx<-idx-dim(ME)[1]
      idx.2<-1
    }
    else{
      idx.2<-2
    }
    return(c(idx,idx.2))
  }
  if (j.cl>0){
    if (j.cl==dim(ME)[1] | j.cl==2*dim(ME)[1]){
      idx=402
    }
    else{
      idx<-which(ME==j.cl)
    }
    if (idx>=dim(ME)[1]+1){
      idx<-idx-dim(ME)[1]
      idx.2<-1
    }
    else{
      idx.2<-2
    }
    return(c(idx,idx.2))
  }
}

#get cluster number from merge matrix of a set of variables
get.number<-function(vars, ME){
  list<-sapply(vars, function(x){get.merge.step(j.var=x, ME=ME)})[1,]
  vars.n<-unique(list)
  while(length(vars.n)>1){
    list<-c(list[-which(list == min(list, na.rm = TRUE))], get.merge.step(j.cl=min(vars.n), ME=ME)[1])
    #list<-sapply(vars.n, function(x){get.merge.step(j.cl=x, ME=ME)})
    vars.n<-unique(list)
    length(vars.n)
  }
  return(vars.n)
}

#plot *layer*th-order neighbourhood of a variable in a large network
#return the neighbourhood subgraph and the neighbour nodes as a vector
graphviz.plot.neighbourhood<-function(net, target=NULL, layer=1){
  if(is.groupbn(net)){
    target<-net$target
    net<-net$bn
  }
  if(layer=="all"){
    x<-bnlearn::graphviz.plot(cpdag(net),shape="ellipse", highlight = list(nodes=target, fill="gray", col="black"))
    invisible(list(plot=x, subnetwork=net, neighbours=nodes(net)))
  }else{
    nodes1=target
    for (i in 1:layer){
      nodes1=c(nodes1, unique(unlist(lapply(nodes1,FUN=function(x) mb(net, x)))))
    }
    if (!target%in%nodes1){
      nodes1<-c(nodes1, target)
    }
    nodes1<-unique(nodes1)
    sub=bnlearn::subgraph(net, nodes1)
    x<-bnlearn::graphviz.plot(bnlearn::cpdag(sub),shape="ellipse", highlight = list(nodes=target, fill="gray", col="black"))
    invisible(list(plot=x, subnetwork=sub, neighbours=nodes1))
  }
}


#group data, discretize and factorize and separate target etc.
#res: groupbn object with data, target, separate, grouping
group.data.preproc<-function(res, seed=NULL){
  if (!is.null(seed)) {set.seed(seed)}
  #combine data
  X.quanti=res$X.quanti
  X.quali=res$X.quali
  target=res$target
  separate=res$separate

  data<-cbind(X.quanti, X.quali)
  if(dim(X.quanti)[1]!=dim(X.quali)[1]){
    stop("wrong data dimensions")
  }

  #discretize aggregation variables
  res<-disc.scores(res, seed=seed)
  cluster<-res$grouping
  k<-max(cluster)
  n<-max(cluster)
  data.scores<-res$group.data

  #separate target and other variables that need to be separated
  if(!is.null(target)){
    if(!(target%in%colnames(data))){
      stop("Could not find target variable")
    } else {
      if (target%in%colnames(X.quanti)){
        target.data<-X.quanti[target]
      } else {
        target.data<-X.quali[target]
      }
        target.number<-which(colnames(data)==target)
        cl.target<-cluster[target]
        cl.target.names<-names(which(cluster==cl.target))
        cl.target.names<-cl.target.names[-which(cl.target.names==target)]
        #do pca
        pc<-PCAmix.groupbn(X.quanti, X.quali, cl.target.names, seed=seed)
        res$pca.param[[cl.target]]<-pc
        #discretize score
        disc<-discretize.dens(pc$scores[,1], cluster=T, seed=seed)
        data.scores[,cl.target]<-as.factor(disc$discretized)
        res$disc.param[[cl.target]]<-disc$levels
        data.scores<-cbind(data.scores, target.data)
        colnames(data.scores)[length(colnames(data.scores))]<-target
        cluster[target]<-n+1
        res$disc.param<-lappend(res$disc.param, NULL)
        res$pca.param<-lappend(res$pca.param, NULL)
        n<-n+1
    }
  }
  if(!is.null(separate)){
    for (sep in 1:length(separate)){
      variable<-separate[sep]
      if(!(variable%in%colnames(data))){
        message("Could not find variable")
      } else {

        res$disc.param<-lappend(res$disc.param, NULL)
        res$pca.param<-lappend(res$pca.param, NULL)

        if (variable%in%colnames(X.quanti)){
          #discretize data
          disc<-discretize.dens(X.quanti[,variable], seed=seed)
          variable.data<-as.factor(disc$discretized)
          res$pca.param[[n+1]]<-NULL
          res$disc.param[[n+1]]<-disc$levels
        } else if (variable%in%colnames(X.quali)){
          variable.data<-X.quali[variable]
        }

        variable.number<-which(colnames(data)==variable)
        cl.variable<-cluster[variable]
        cl.variable.names<-names(which(cluster==cl.variable))
        cl.variable.names<-cl.variable.names[-which(cl.variable.names==variable)]
        #do pca
        pc<-PCAmix.groupbn(X.quanti, X.quali, cl.variable.names, seed=seed)
        res$pca.param[[cl.variable]]<-pc
        #discretize score
        disc<-discretize.dens(pc$scores[,1], cluster=T, seed=seed)
        data.scores[,cl.variable]<-as.factor(disc$discretized)
        res$disc.param[[cl.variable]]<-disc$levels

        data.scores<-cbind(data.scores, variable.data)
        colnames(data.scores)[length(colnames(data.scores))]<-variable
        cluster[variable]<-n+1
        n<-n+1
      }
    }
  }

  for (i in 1:k){
    if(is.null(res$pca.param[[i]])){
      cl.variable.names<-names(which(cluster==i))
      pc<-PCAmix.groupbn(X.quanti, X.quali, cl.variable.names, seed=seed)
      res$pca.param[[i]]<-pc
      #discretize score
      disc<-discretize.dens(pc$scores[,1], cluster=T, seed=seed)
      #data.scores[,cl.variable]<-as.factor(disc$discretized)
      res$disc.param[[i]]<-disc$levels
    }
  }
  res$group.data<-data.scores
  res$grouping<-cluster
  res$k<-k
  return(res)
}

#add lists to each other
lappend <- function (lst, x){
  lst <- c(lst, list(x))
  return(lst)
}

#add elements of list of lists into separate lists
#rearrange.listoflists <- function (lst){
#  outcome<-vector(mode = "list", length = length(lst[[1]]))
#for (i in 1:length(lst[[1]])){
#    outcome[[i]]<-vector(mode = "list", length = length(lst))
#    for(j in 1:length(lst)){
#      outcome[[i]][[j]]<-lst[[j]][[i]]
#    }
#  }
#  return(outcome)
#}

#get nth element of several lists as a list
comb.output.scores<-function(x,...,n=13){
  sapply(seq_along(x), FUN=function(i) c(rlist::list.last(x[[i]][[n]][[1]])))
}

#Learn Network
#data: Dataframe without missing values
network<-function(data, start=NULL, R=200, restart=5, perturb=max(1,round(0.1*dim(data)[2])), blacklist=NULL, debug=F, seed=NULL){
  if(debug){message("learning arcs")}
  if(!is.null(seed)){set.seed(seed)}
  arcs<-suppressWarnings(bnlearn::boot.strength(data, R=R, algorithm = "hc", algorithm.args = list(restart=restart, perturb=perturb, start=start, blacklist=blacklist)))
  net<-bnlearn::averaged.network(arcs)
  res <- try(bnlearn::cextend(net), silent = T)
  #try if extension is possible
  while(inherits(res, "try-error")){
    message("try-error")
    arcs<-bnlearn::boot.strength(data, R=R, algorithm = "hc", algorithm.args = list(restart=restart, perturb=perturb))
    net<-bnlearn::averaged.network(arcs)
    res <- try(bnlearn::cextend(net))
  }
  return(list(net=net, arc.confid<-arcs))
}

#Create an output table with clusters and included variables with similarity scores
#res: groupbn object
#save.name: filename for saving
#pdf: Boolean, should the file be saved as pdf
groupbn.output.table<-function(res){
  net<-res$bn
  X.quali<-res$X.quali
  X.quanti<-res$X.quanti
  data.scores<-res$group.data
  cluster<-res$grouping

  data<-cbind(X.quanti, X.quali)
  #Calculate similarity scores of synthetic variables and original variables
  #calculate similarity score to cluster representant
  scores<-cluster
  names(scores)<-colnames(data)
  for (i in 1:length(scores)){
    if(nlevels(factor( data[stats::complete.cases(data.scores)&stats::complete.cases(data),i]))==1|nlevels(factor(data.scores[stats::complete.cases(data.scores)&stats::complete.cases(data),cluster[i]]))==1){
      scores[i]<-NA
    } else if(is.factor(data[stats::complete.cases(data.scores)&stats::complete.cases(data),i])){
      scores[i]<-ClustOfVar::mixedVarSim(factor(data.scores[stats::complete.cases(data.scores)&stats::complete.cases(data),cluster[i]]), factor(data[stats::complete.cases(data.scores)&stats::complete.cases(data),i]))
    } else {
      scores[i]<-ClustOfVar::mixedVarSim(factor(data.scores[stats::complete.cases(data.scores)&stats::complete.cases(data),cluster[i]]), data[stats::complete.cases(data.scores)&stats::complete.cases(data),i])
    }
  }

  #create data frame with information about cluster and similarity
  names<-paste("cluster",1:max(cluster), sep="")
  max.len<-max(summary(factor(cluster)))
  #first column
  df<-data.frame(c(names(sort(scores[names(which(cluster==1))], decreasing=T, na.last=TRUE)), rep("", max.len - length(names(which(cluster==1))))),stringsAsFactors =FALSE)
  colnames(df)<-names[1]
  #other columns
  for (i in 2:max(cluster)){
    x<-names(sort(scores[names(which(cluster==i))], decreasing=T, na.last=TRUE))
    df[,i]<-c(x, rep("", max.len - length(x)))
    colnames(df)[i]<-names[i]
  }

  #add scores
  for (i in 1:max.len){
    for (j in 1:max(cluster)){
      if (df[i,j]!=""){
        if(is.na(scores[df[i,j]])){
          df[i,j]<-paste(df[i,j], ": NA", sep="")
        }
        else if (round(scores[df[i,j]],2)==0){
          df[i,j]<-paste(df[i,j], ": <0.01", sep="")
        } else{
          df[i,j]<-paste(df[i,j], ": ",round(scores[df[i,j]],2), sep="")
        }
      }
    }
  }
  invisible(df)
}

#do PCA for a given vector of variables (splitting, PCAmix call with right configuration)
#X.quanti, X.quali: data
#names: variable names for which a pca should be done
#graph: if PCA plots should be produced
PCAmix.groupbn<-function(X.quanti, X.quali, names, graph=FALSE, seed=NULL){
  #split to quali and quantitative variables
  if (!is.null(seed)) {set.seed(seed)}
  quali.names<-names[which(names%in%colnames(X.quali))]
  quanti.names<-names[which(names%in%colnames(X.quanti))]
  if (length(quanti.names)==1&length(quali.names)>0){
    #1 Quanti, many Quali: Error: discretize
    quan2qual<-discretize.dens(X.quanti[,quanti.names], seed=seed)$discretized
    quan2qual<-as.data.frame(quan2qual)
    colnames(quan2qual)<-quanti.names
    pc<-PCAmixdata::PCAmix(X.quali=cbind(X.quali[,quali.names, drop=F], quan2qual), rename.level = T, graph=graph, ndim=2)
  } else if (length(quanti.names)>0&length(quali.names)==0){
    #Only Quanti
    pc<-PCAmixdata::PCAmix(X.quanti=X.quanti[,quanti.names, drop=F], graph=graph, rename.level=TRUE)
  } else if (length(quanti.names)==0&length(quali.names)>0){
    #Only Quali
    pc<-PCAmixdata::PCAmix(X.quali=X.quali[,quali.names, drop=F], graph=graph, rename.level=TRUE)
  } else if (length(quanti.names)>0&length(quali.names)>0){
    #both
    df1<-X.quanti[,quanti.names, drop=F]
    df2<-X.quali[,quali.names, drop=F]
    pc<-PCAmixdata::PCAmix(X.quanti=df1, X.quali=df2, graph=graph, rename.level=TRUE)
  }
  return(pc)
}

plot.groupbn<-function(x,...){
  plot(bnlearn::cpdag(x$bn), highlight=c(x$target), color="coral1",...)
  print(x)
}

#remove variables where all the categories are identical
rem.id.col<-function(X.quali){
  t<-sapply(X.quali, table)
  vars.rem<-c()
  for (i in 1:length(t)){
    if(sum(t[[i]]==0)>=length(t[[i]])-1){
      vars.rem<-c(vars.rem, names(t)[i])
    }
  }
  return(vars.rem)
}

#get division of one cluster to the next two clusters, vars=indices of variables
split.cluster<-function(vars, ME){
  nr<-get.number(vars, ME)
  split.nr<-ME[nr,]
  cl1<-get.cluster(split.nr[1], ME)
  cl2<-get.cluster(split.nr[2], ME)
  return(list(cl1, cl2))
}

#Parameter fitting and several validation scores
#if do.fit=T, fitting is done within the function
#if do.fit=F, net must be an already fitted object of class bn.fit , method %in% c("parents", "bayes-lw)
#n bootstapping for probability estimation
validation<-function(net, data, target, thresh=0.5, do.fit=TRUE, n=2000 ,method="bayes-lw", debug=F, seed=seed){
  if (!is.null(seed)) {set.seed(seed)}
  if(do.fit){
    if(debug){message("fitting")}
      if (is.factor(data[[target]])){
        fit<-bnlearn::bn.fit(cextend(net), data, method="bayes")
      } else {
        fit<-bnlearn::bn.fit(cextend(net), data)
     }
 } else {
    fit=net
    for (clus in 1:bnlearn::nnodes(fit)){
      if (any(levels(data[,clus])!=rownames(fit[[clus]][[4]]))){
      message("Attention: Level names changed of cluster", clus)
      levels(data[,clus])<-rownames(fit[[clus]][[4]])
      }
    }
 }
  if(debug){message("validation")}
  #save actual values
  obs<-data[[target]]
  if(!is.factor(obs)){
    temp<-sapply(1:dim(data)[1], function(j){
      if(sum(is.na(data[j, -which(colnames(data)==target)]))==0){
        prob<-replicate(10, predict(object=fit, node=target, data=data[j, -which(colnames(data)==target)], method = method, n=n))
        return(c(mean(prob), stats::sd(prob)))
      }
      else{
        return(c(NA, NA))
      }
    })
    p<-temp[1,]
    sdprob<-temp[2,]
    rm(temp)

    #transform probabilities to predictions
    result<-MLmetrics::MSE(p, obs)
    attr(result, "pscores")<-p
    psd<-p
    psd[p>obs]<-psd[p>obs]-sdprob[p>obs]
    psd[p<=obs]<-psd[p<=obs]+sdprob[p<=obs]
    result.error<-MLmetrics::MSE(psd, obs)
    attr(result, "error.th")<-result.error+(abs(result-result.error))/2
    x<-paste("cor: ",round(cor(p, obs),2),"; R-sq: ", round(cor(p, obs)^2,2))
    if(debug){message(x)}
    attr(result, "scores") <- x
    return(result)
    #stop('implemented only for a categorical target variable with 2-4 levels. Haha')
  } else if(nlevels(factor(obs))<2|nlevels(factor(obs))>4){
    stop('implemented only for a categorical target variable with 2-4 levels')
  } else if(nlevels(factor(obs))>2){
    k<-nlevels(obs)
    levels(obs)<-0:(k-1)
    p<-matrix(0, dim(data)[1], k)
    for (j in 1:dim(data)[1]){
      if(sum(is.na(data[j, -which(colnames(data)==target)]))==0){
        p[j,]<-attr(predict(object=fit, node=target, data=data[j, -which(colnames(data)==target)], method = method, n=n, prob = TRUE), "prob")
      }
      else{
        p[j,]<-rep(NA,k)
      }
    }
    #transform probabilities to predictions
    predictions<-apply(p,1, which.max)-1
    predictions<-factor(predictions, levels=c(0:k))

    #calculate results
    result<-cross.en(p, obs, weighted=T)
    attr(result, "scores") <- NA
    attr(result, "auc")<- NA
    attr(result, "sd")<- NA
    attr(result, "confusion") <- table(obs, predictions)
    attr(result, "pscores")<-p
    if(debug){print(table(predictions,obs))}
    return(result)
  } else {
  levels(obs)<-c(0,1)
  temp<-sapply(1:dim(data)[1], function(j){
    if(sum(is.na(data[j, -which(colnames(data)==target)]))==0){
      prob<-replicate(10, attr(predict(object=fit, node=target, data=data[j, -which(colnames(data)==target)], method = method, n=n, prob = TRUE), "prob")[1])
      return(c(mean(prob), stats::sd(prob)))
    }
    else{
      return(c(NA, NA))
    }
  })
  p<-temp[1,]
  sdprob<-temp[2,]
  rm(temp)

  #transform probabilities to predictions
  predictions<-p
  predictions[p>=thresh]<-0
  predictions[p<thresh]<-1
  predictions<-factor(predictions, levels=c(0,1))

  #calculate different scores
  f11<-MLmetrics::F1_Score(y_true=obs[!is.na(predictions)], y_pred = predictions[!is.na(predictions)], positive=1)
  pre1<-MLmetrics::Precision(y_true=obs[!is.na(predictions)], y_pred = predictions[!is.na(predictions)], positive=1)
  re1<-MLmetrics::Recall(y_true=obs[!is.na(predictions)], y_pred = predictions[!is.na(predictions)], positive=1)
  pr.cu<-PRROC::pr.curve(scores.class0 = 1-p[!is.na(p)], weights.class0 = as.numeric(levels(obs[!is.na(p)]))[obs[!is.na(p)]], curve=T, rand.compute = T, min.compute = T, max.compute = T)
  #plot(pr.cu,max.plot = TRUE, min.plot = TRUE, rand.plot = TRUE, fill.area = TRUE, auc.main=TRUE)
  pr1<-pr.cu$auc.integral
  roc.cu<-PRROC::roc.curve(scores.class0 = 1-p[!is.na(p)], weights.class0 = as.numeric(levels(obs[!is.na(p)]))[obs[!is.na(p)]], curve=T, rand.compute = T, min.compute = T, max.compute = T)
  #plot(roc.cu,max.plot = TRUE, min.plot = TRUE, rand.plot = TRUE, fill.area = TRUE, auc.main=TRUE)
  pr2<-roc.cu$auc
  result<-cross.en(p, obs)
  result.error<-cross.en(p, obs, weighted=T,sdpred=sdprob)
  x<-paste("F1: ",round(f11,2),"; Precision: ", round(pre1,2), "; Recall: ", round(re1,2), "; AUC-PR: ", round(pr1,3), "; AUC-ROC: ", round(pr2,3), "; cross-entr.: ", round(result,2))
  if(debug){message(x)}
  attr(result, "scores") <- x
  attr(result, "auc")<-roc.cu$auc
  attr(result, "error.th")<-result.error+(result-result.error)/2
  attr(result, "confusion")<-table(obs, predictions)
  attr(result, "pscores")<-p
  if(debug){print(table(predictions,obs))}
  return(result)
  }
}

#Create an interactive html network object with visNet (displaying similarity scores and number of variables in a score)
#res: groupbn object
#df: data frame from groupbn.output.table
#save.name: Name for file
groupbn.vis.html.plot<-function(res, df=NULL, save.file=TRUE, save.name=NULL, hierarchical=FALSE, nodecolor.all="#E0F3F8", nodecolor.special="black", main=NULL){
  if (is.null(df)){
    df<-groupbn.output.table(res)
  }
  if(is.null(main)){
    main=paste("Network model of ",res$target)
  }
  net<-res$bn
  cluster<-res$grouping
  arrows<-rep("to", dim(bnlearn::arcs(net))[1])
  undir<-bnlearn::undirected.arcs(bnlearn::cpdag(net))
  for (i in 1:dim(bnlearn::arcs(net))[1]){
    x<-bnlearn::arcs(net)[i,]
    set<-intersect(which(x[1]==undir[,1]),which(x[2]==undir[,2]))
    if(length(set)>0){
      arrows[i]<-"enabled"
    }
  }
  #node titles
  title<-paste("<b>Cluster", 1:max(cluster), "</b>")
  for (i in 1:dim(df)[2]){
    for (j in 1:dim(df)[1])
      if(df[j,i]!=""){
        title[i]<-paste(title[i], "<br>", df[j,i])
      }
  }
  #use most representative variables as names
  names_P2<-df[1,]
  for (i in 1:length(names_P2)){
    names_P2[i]<-paste("[",i, "] ", stringr::str_split(names_P2[i], stringr::fixed(":"))[[1]][1], " (",length(which(cluster==i)),")", sep="")
  }
  #shorter cluster names with number of variables
  names_P3<-names_P2
  for (i in 1:length(names_P2)){
    names_P3[i]<-paste("Cl", i, " (", length(which(cluster==i)), ")", sep="")
  }
  names_P3[which(bnlearn::nodes(net)==res$target)]<-toupper(res$target)
  #nodecolors
  #all blue
  nodecolors<-rep(nodecolor.all, length(bnlearn::nodes(net)))
  #special nodes: other blue
  nodecolors[which(bnlearn::nodes(net)%in%c(res$target, res$separate))]<-nodecolor.special
  #node settings
  nodes <- data.frame(id = bnlearn::nodes(net),
                      shape = c(rep("dot", length(nodes(net)))),#shape of nodes
                      title = as.character(title),                #node title
                      color = nodecolors,                         #colors
                      label = as.character(names_P2),             #node labels
                      font=list(color="black"))                   #font color
                       #shadow = c(FALSE, TRUE)                   #shadow
  strength<-round(res$arc.confid$strength[as.numeric(rownames(plyr::match_df(as.data.frame(res$arc.confid[,1:2]), data.frame(sapply(as.data.frame(arcs(res$bn)), as.character), stringsAsFactors=FALSE), on=c("from", "to"))))],2)
  #edge settings
  edges <- data.frame(from = bnlearn::arcs(net)[,1],
                      to = bnlearn::arcs(net)[,2],
                      arrows=arrows,
                      color=rep("darkgray", dim(bnlearn::arcs(net))[1]),
                      font.color = "gray",
                      title=strength,
                      value=strength,
                      scaling=list(min=5, max=10))
  #build network
  if(hierarchical){
    graph<-visNetwork::visNetwork(nodes, edges, width = "100%", height="1000px", main = main,
                      footer=paste("Group Bayesian Network.\n The thickness of an edge represents its confidence.", sep="")) %>%
      visNetwork::visHierarchicalLayout(blockShifting=TRUE, edgeMinimization=FALSE)%>%
      visNetwork::visOptions(highlightNearest = list(enabled = T, hover = T),nodesIdSelection = T)%>%
      visNetwork::visInteraction(navigationButtons = TRUE)

  } else {
    graph<-visNetwork::visNetwork(nodes, edges, width = "100%", height="1000px", main = main,
                      footer=paste("Group Bayesian Network.\n The thickness of an edge represents its confidence. Clusters are named by its most central variable. The number is brackets denotes the number of variables in the cluster.", sep="")) %>%
      visNetwork::visIgraphLayout(type = "square")%>%
      visNetwork::visOptions(highlightNearest = list(enabled = T, hover = T),nodesIdSelection = T)%>%
      visNetwork::visInteraction(navigationButtons = TRUE)
  }
  if(save.file){
    #save as html object
    if(is.null(save.name)){
      save.name=paste(format(Sys.Date(), format="%Y%m%d"),"_GroupBN", sep="")
    }
    visNetwork::visSave(graph, file=paste(save.name, ".html", sep=""), selfcontained = TRUE, background = "white")
    message("html written. (", save.name, ".html)")
    invisible(graph)
   } else {
     #optionally: output in RStudio Viewer
     print(graph)
     invisible(graph)
   }
}

#is.groupbn
is.groupbn <- function(x){
  inherits(x, "groupbn")
}

#print.groupbn
print.groupbn<-function(x,...){
  if (!inherits(x, "groupbn"))
    stop("use only with \"groupbn\" objects")
    gnr <- dim(x$group.data)[2]
    separate<-x$separate
    cat("group Bayesian network (class 'groupbn') \n\n")
    cat("name of target variable: ", x$target, "\n",sep = "")
    if(!is.null(separate)){
      cat("separated: ", paste(separate, collapse=" & "), "\n",sep = "")
    }
    if (is.numeric(gnr)){
      cat("number of groups: ", gnr,"\n", sep = " ")
    }
    cat("achieved scoring: ", attr(x$score, "scores"), "; BIC (netw.): ", round(stats::BIC(x$fit, stats::na.omit(x$group.data)),2), "\n",sep = "")
    cat("\n")
  res <- matrix("", 13, 2)
  colnames(res) <- c("name", "description")
  res[1, ] <- c("$bn", "Bayesian network structure")
  res[2, ] <- c("$fit", "fitted Bayesian network (multinomial)")
  res[3, ] <- c("$arc.confid", "arc confidence")
  res[4, ] <- c("$X.quali", "qualitative variables in a data.frame")
  res[5, ] <- c("$X.quanti", "quantitative variables in a data.frame")
  res[6, ] <- c("$grouping", "group memberships")
  res[7, ] <- c("$k", "number of groups of initial grouping")
  res[8, ] <- c("$group.data", "group representatives used for network inference")
  res[9, ] <- c("$target", "name of target variable")
  res[10,] <- c("$separate", "name of any other separated variables")
  res[11, ] <- c("$pca.param", "pca parameters of each group")
  res[12, ] <- c("$disc.param", "discretization intervals of each group")
  res[13, ] <- c("$score", "cross entropy and additional scoring information")
  row.names(res) <- rep("", 13)
  print(res,...)
}


