
#create groupbn object (group network) (determines initial clustering with target, separate separated, learns network from group variables and saves disc and pca params)
#hierarchy: clustofvar clustering
#k: number of clusters (where to cut the dendrogram)
#target, separate
#X.quanti, X.quali: data
groupbn<-function(hierarchy, k, target, separate=NULL, separate.as.roots=FALSE, X.quanti=NULL, X.quali=NULL, struct.alg="hc", boot=TRUE, discretize=TRUE, arc.thresh=NULL, debug=FALSE, R=100, seed=NULL){
  start_time <- Sys.time()

  if(!is.null(X.quanti)&!is.null(X.quali)&dim(X.quanti)[1]!=dim(X.quali)[1]){
    stop("The number of samples of X.quali and X.quanti is different")
  }
  if(!is.null(X.quanti)){
    names.quanti<-colnames(unique(as.matrix(X.quanti),MARGIN=2))
    if(length(names.quanti)!=dim(X.quanti)[2]){
      warning(paste("Identical columns in X.quanti:", paste(colnames(X.quanti)[which(!colnames(X.quanti)%in%names.quanti)], collapse=" ")))
    }
  }
  if(!is.null(X.quali)){
    names.quali<-colnames(unique(as.matrix(X.quali),MARGIN=2))
    if(length(names.quali)!=dim(X.quali)[2]){
      warning(paste("Identical columns in X.quali:", paste(colnames(X.quali)[which(!colnames(X.quali)%in%names.quali)], collapse=" ")))
    }
  }

  data<-cbind(X.quanti, X.quali)

  #create gn object
  res <- list(bn = NULL,
              fit = NULL,
              arc.confid=NULL,
              X.quali = X.quali,
              X.quanti=X.quanti,
              grouping=NULL,
              k=k,
              group.data=NULL,
              target=target,
              separate=separate,
              pca.param=vector(mode="list", length=k),
              disc.param=vector(mode="list", length=k),
              score=NULL)
  class(res) <- "groupbn"
  if(!is.null(separate)){
  attributes(res$separate)$"separate.as.roots"<-separate.as.roots
  }else{
    separate.as.roots<-FALSE
  }

  if (!(target%in%colnames(X.quanti)|target%in%colnames(X.quali))){
    stop("target variable not found in data")
  }
  if(debug){message("grouping")}
  #determine initial clustering, cut dendrogram
  P<-ClustOfVar::cutreevar(hierarchy, k=k)
  res$grouping<-P$cluster
    if(sum(is.na(names(res$grouping)))>0){
      names(res$grouping)<-c(names.quanti, names.quali)
    }
  res$group.data<-P$scores
  for (i in 1:length(res$pca.param)){
      pc<-PCAmix.groupbn(X.quanti, X.quali, names(which(res$grouping==i)), seed=seed)
      res$pca.param[[i]]<-pc
      res$group.data[,i]<-pc$scores[,1]
  }
  #data preprocessing
  res<-group.data.preproc(res, discretize=discretize, seed=seed)
  #determine variables
  data.scores<-res$group.data
  cluster<-res$grouping
  n<-max(cluster)
  pca.param<-res$pca.param
  names(res$pca.param)<-colnames(res$group.data)
  disc.param<-res$disc.param

  target.number<-which(colnames(data)==target)
  variable.number<-which(colnames(data)%in%separate)


  ##############################
  ####Learn initial network#####
  ##############################

  #learn a network from complete cases
  blacklist<-NULL
  if(!is.null(separate)&any(separate.as.roots)){
    blacklist<-groupbn.build.blacklist(data.scores, separate[separate.as.roots])
  }
  net<-network(data.scores[complete.cases(data.scores),], R=R, seed=seed, debug=debug, blacklist=blacklist, arc.thresh=arc.thresh, struct.alg=struct.alg, boot=boot)
  res$bn<-net[[1]]
  res$arc.confid<-net[[2]]

  #Prediction
  if (!is.null(seed)) {set.seed(seed)}
  if(debug){message("fitting")}
  if (target%in%colnames(X.quanti)|discretize==FALSE){
    fit<-bnlearn::bn.fit(bnlearn::cextend(res$bn), data.scores[complete.cases(data.scores),])
  } else if (target%in%colnames(X.quali)){
    fit<-bnlearn::bn.fit(bnlearn::cextend(res$bn), data.scores[complete.cases(data.scores),], method="bayes")
  }

  result<-validation(fit, data.scores[complete.cases(data.scores),], target, do.fit=FALSE, seed=seed, debug=debug)
  res$fit<-fit
  res$score<-result

  attributes(res$bn)$struct.alg<-struct.alg
  attributes(res$bn)$boot<-boot

  #timing
  end_time <- Sys.time()
  duration=end_time-start_time
  if(debug){message("Calculation Time: ", round(duration[[1]], 2),  units(duration))}

  return(res)
}


groupbn_refinement<-function(res, hierarchy, refinement.part="mb", restart=0, perturb=1, max.step=10, max.min=Inf, R=100, return.all=FALSE, arc.thresh=NULL, debug=FALSE, seed=NULL){
  #hierarchy: cluster object from ClustOfVar
  #res: groupbn object
  #perturb: perturbations per restart
  #restart: number of restarts after a local optima
  #max.step: maximal number of optimizing steps per restart
  #refinement.part: "mb", "mb2" or "all" to refine only within markov blanket or all clusters

  #timing
  start_time <- Sys.time()

  ######Initialization#####
  if (!is.null(seed)) {set.seed(seed)}
  k<-res$k
  target<-res$target
  separate<-res$separate
  X.quali<-res$X.quali
  X.quanti<-res$X.quanti
  data<-cbind(X.quanti, X.quali)
  net<-res$bn
  arc.confid<-res$arc.confid
  cluster<-res$grouping
  data.scores<-res$group.data
  result<-res$score
  target.number<-which(colnames(data)==target)
  variable.number<-which(colnames(data)%in%separate)
  pca.param<-res$pca.param
  disc.param<-res$disc.param
  n<-max(cluster)
  struct.alg<-attributes(res$bn)$struct.alg
  boot<-attributes(res$bn)$boot
  data<-cbind(X.quanti, X.quali)
  if(!is.null(separate)){
    separate.as.roots<-attr(res$separate, "separate.as.roots")
  } else {
    separate.as.roots<-FALSE
  }
  if(all(vapply(res$disc.param, is.null, TRUE))){
    discretize<-FALSE
  } else {discretize<-TRUE}

  ##############################
  ####Adaptive Refinement#######
  ##############################

  #get merge matrix from cluster object
  ME<-hierarchy$merge

  #save result for every restart in a list
  restart.list<-vector(mode = "list", length = restart+2)
  restart.list[[1]]<-res
  #loop for every restart
  print("here")
  for (re in 0:restart){
    #initialize result list
    result.list<-list()

    #After first run: add random perturbation
    if (re>0){
      #start with best model and split *perturb* random cluster of the markov blanket
      if(re==1){
        minimum=2
      }
      else{
        last<-comb.output.scores(restart.list)
        if(length(which(sapply(last, is.null)))>0){
          last[which(sapply(last, is.null))]<-Inf
        }
        minimum<-which.min(unlist(last))
      }

      #get optimal model
      message("minimum:", minimum)
      net<-restart.list[[minimum]]$bn
      arc.confid<-restart.list[[minimum]]$arc.confid
      data.scores<-restart.list[[minimum]]$group.data
      cluster<-restart.list[[minimum]]$grouping
      pca.param<-restart.list[[minimum]]$pca.param
      disc.param<-restart.list[[minimum]]$disc.param
      #debugging
      if(debug){message("restart", re, ":")}

      #Determine clusters in markov blanket
      if(refinement.part=="mb"){
        listnod<-bnlearn::mb(net, target)
        listnod<-c(listnod, bnlearn::nodes(net)[-which(bnlearn::nodes(net)%in%c(target, separate))])
        listnod<-unique(listnod)
        if(length(listnod)==0){
          stop("You chose the Markov Blanket for refinement, but it is empty. Use 'refinement.part='all'' instead.")
        }
      }
      if(refinement.part=="mb2"){
        listnod<-bnlearn::mb(net, target)
        for(g in listnod){
          listnod<-c(listnod, bnlearn::mb(net, g))
        }
        listnod<-c(listnod, bnlearn::nodes(net)[-which(bnlearn::nodes(net)%in%c(target, separate))])
        listnod<-unique(listnod)
        if(length(listnod)==0){
          stop("You chose the Markov Blanket for refinement, but it is empty. Use 'refinement.part='all'' instead.")
        }
      }
      if(refinement.part=="all"){
        listnod<-bnlearn::nodes(net)[-which(bnlearn::nodes(net)%in%c(target, separate))]
        listnod<-unique(listnod)
      }
      if(refinement.part=="arc.confid"){      #based on arc confidence > 0
        listnod<-arc.confid[arc.confid$from==target&arc.confid$strength>0,2]
      }

      #remove age and sex from the list (cannot be split)
      if(any(listnod%in%c(separate,target))){
        listnod<-listnod[-which(listnod%in%c(separate,target))]
      }
      cl.nrs<-substr(listnod, rep(3, length(listnod)), nchar(listnod))
      remove<-c()
      for (nod in 1:length(cl.nrs)){
        if(length(which(cluster==cl.nrs[nod]))==1){
          remove<-c(remove, nod)
        }
      }
      if(length(remove)){
        listnod<-listnod[-remove]
      }
      #get the number of perturbations
      nr.perturbations<-min(length(listnod), perturb)

      #sample which clusters are split
      if(refinement.part=="mb"){
        if(nr.perturbations<=length(bnlearn::mb(net, target))){
          pert<-sample(listnod[1:length(bnlearn::mb(net, target))], nr.perturbations)
        } else {
          pert<-sample(listnod[1:length(bnlearn::mb(net, target))], length(bnlearn::mb(net, target)))
          pert<-c(pert, sample(listnod[(length(bnlearn::mb(net, target))+1):length(listnod)], nr.perturbations-length(bnlearn::mb(net, target))))
        }
      } else {
        pert<-sample(listnod, nr.perturbations)
      }
      if(debug){message(nr.perturbations, " perturbation(s) are introduced: ", paste(pert , collapse=" "))}
      #loop through number of perturbations and introduce them
      for (pert.idx in 1:nr.perturbations){
        #which cluster is split
        cl<-substr(pert[pert.idx], 3, nchar(pert[pert.idx]))

        #set number of clusters
        n<-max(cluster)

        pca.param<-lappend(pca.param, NULL)
        disc.param<-lappend(disc.param, NULL)
        #variable names
        vars<-names(cluster[which(cluster==cl)])
        if(debug){message("Variables: ", paste(vars, collapse=" "))}

        #transform to numbers
        vars<-which(colnames(data)%in%vars)
        #only one variable: skip, not splittable
        if (length(vars)==1){
          if(debug){message("not splittable, skip")}
          next
        }

        #two variables, one is sex/age/target: skip
        if(length(vars)==2&any(vars%in%c(separate,target))){
          if(debug){message("not splittable, skip")}
          next
        }

        #otherwise determine splitting from merge matrix
        splits<-split.cluster(vars=vars, ME=ME)
        #first split
        j=1
        if(length(which(colnames(data)[splits[[j]]]%in%c(separate,target)))){
          splits[[j]]<-splits[[j]][-which(colnames(data)[splits[[j]]]%in%c(separate,target))]
        }
        #if only one variable in cluster
        if(length(splits[[j]])==1){
          if(is.factor(data[,splits[[j]]])){
            #if it is already a factor
            data.scores[,as.numeric(cl)]<-data[,splits[[j]]]
            pca.param[[as.numeric(cl)]]<-NULL
            disc.param[[as.numeric(cl)]]<-NULL
          }
          else{
            #if it is continuous: discretize
            if(discretize){
              disc<-discretize.dens(data[,splits[[j]]], cluster=F, seed=seed)
              data.scores[,as.numeric(cl)]<-as.factor(disc$discretized)
              pca.param[[as.numeric(cl)]]<-NULL
              disc.param[[as.numeric(cl)]]<-disc$levels
            } else {
              data.scores[,as.numeric(cl)]<-data[,splits[[j]]]
            }
          }
        }else{
          names<-colnames(data)[splits[[j]]]
          pc<-PCAmix.groupbn(X.quanti, X.quali, names, seed=seed)
          pca.param[[as.numeric(cl)]]<-pc
          if(discretize){
          disc<-discretize.dens(pc$scores[,1], cluster=T, seed=seed)
          disc.param[[as.numeric(cl)]]<-disc$levels
          #substitute column from before
          data.scores[,as.numeric(cl)]<-as.factor(disc$discretized)
          } else {
            data.scores[,as.numeric(cl)]<- pc$scores[,1]
          }
        }
        #second split
        j=2
        if(length(which(colnames(data)[splits[[j]]]%in%c(separate,target)))){
          splits[[j]]<-splits[[j]][-which(colnames(data)[splits[[j]]]%in%c(separate,target))]
        }
        if(length(splits[[j]])==1){
          if(is.factor(data[,splits[[j]]])){
            #if it is already a factor
            data.scores[,n+1]<-data[,splits[[j]]]
          }
          else{
            if(discretize){
              disc<-discretize.dens(data[,splits[[j]]], cluster=F, seed=seed)
              data.scores[,n+1]<-as.factor(disc$discretized)
              disc.param[[n+1]]<-disc$levels
            } else {
              data.scores[,n+1]<-data[,splits[[j]]]
            }
          }
        }
        else{
          names<-colnames(data)[splits[[j]]]
          pc<-PCAmix.groupbn(X.quanti, X.quali, names, seed=seed)
          pca.param[[n+1]]<-pc
          #add column
          if(discretize){
            disc<-discretize.dens(pc$scores[,1], cluster=T, seed=seed)
            data.scores[,n+1]<-as.factor(disc$discretized)
            disc.param[[n+1]]<-disc$levels
          } else {
            data.scores[,n+1]<-pc$scores[,1]
          }
        }
        colnames(data.scores)[dim(data.scores)[2]]<-paste("cl", n+1, sep="")
        cluster[splits[[j]]]<-n+1
        n<-max(cluster)
      }
      #if(cluster==cluster.old){
      # if(debug){message("clusters unchanged")}
      #  next
      #}
      #learn network for perturbed cluster structure
      if(!is.null(separate)&any(separate.as.roots)){
        blacklist<-groupbn.build.blacklist(data.scores[complete.cases(data.scores),], separate[separate.as.roots])
      }
      netw<-network(data.scores[complete.cases(data.scores),], R=R, seed=seed, debug=F, blacklist=blacklist, arc.thresh=arc.thresh, struct.alg=struct.alg, boot=boot)
      net<-netw[[1]]
      arc.confid<-netw[[2]]
      result<-validation(net, data.scores, target, seed=seed, debug=F)
      result.list<-lappend(result.list, result)
      if(debug) {message("Score: ", round(result,2), ", threshold: ", round(attr(result, "error.th"),2))}
    }


    #Do steps until max.step
    for (step in 1:(max.step+1)){
      if(Sys.time()>start_time+max.min*60){
        message("Did not converge in ", max.min, " minutes. Stopped.")
        break
      }
      if(debug) {message(".......................................................")}
      if(debug) {message("Iteration ", step)}
      if(step==max.step+1){
        message("Stopped after the maximum of ", max.step, " iterations. Set max.step higher for full optimization.")
        break
      }

      #Which clusters to try??
      #only markov blanket
      if(refinement.part=="mb"){
      listnod<-bnlearn::mb(net, target)
      if(length(listnod)==0){
        stop("You chose the Markov Blanket for refinement, but it is empty. Use 'refinement.part='all'' instead.")
      }
      } else if(refinement.part=="mb2"){#Second Layer Markov Blanket
        listnod<-bnlearn::mb(net, target)
        for(g in listnod){
          listnod<-c(listnod, bnlearn::mb(net, g))
          }
        listnod<-unique(listnod)
        if(length(listnod)==0){
          stop("You chose the Markov Blanket for refinement, but it is empty. Use 'refinement.part='all'' instead.")
        }
      } else if(refinement.part=="all"){      #all
      listnod<-bnlearn::nodes(net)[-which(bnlearn::nodes(net)%in%c(target, separate))]
      listnod<-unique(listnod)
      } else if(refinement.part=="arc.confid"){      #based on arc confidence > 0
        listnod<-arc.confid[arc.confid$from==target&arc.confid$strength>0,2]
      }
      #add markov blanket of markov blanket
      #for(g in listnod){
      #  listnod<-c(listnod, mb(net, g))
      #}
      #listnod<-unique(listnod)

      #remove target from this list
      if(any(listnod%in%c(separate, target))){
        listnod<-listnod[-which(listnod%in%c(target, separate))]
      }
      if(length(listnod)==0){
        stop("All variables in the refinement part are already single variables and cannot be split.")
      }
      l<-length(listnod)
      if(debug) {message("Splits to be tested: ", paste(listnod, collapse =" "))}

      eval.neighbourhood<-function(node){
        if(debug){
          message("\n")
          message("Testing to split ", node, ":")
        }
        #get cluster number
        cl<-substr(node, 3, nchar(node))

        #save in preliminary variables
        cluster[target]<-k+1
        #for (i in 1:length(separate)){
        #  cluster[separate[i]]<-k+1+i
        #}
        cluster1<-cluster
        n<-max(cluster)
        data.scores1<-data.scores
        pca.param1<-lappend(pca.param, NULL)
        disc.param1<-lappend(disc.param, NULL)

        #get variable names
        vars<-names(cluster1[which(cluster1==cl)])
        if(debug){message("Variables: ", paste(vars, collapse=" "))}
        #transform to numbers
        vars<-which(colnames(data)%in%vars)

        #only one variable: next
        if (length(vars)==1){
          if(debug){message("not splittable, skip")}
          return(list(NA, NA, NA, NA, NA, NA, NA))
        } else {
          #do splitting
          splits<-split.cluster(vars=vars, ME=ME)
          #check if target is contained
          if (length(splits[[1]])==1&any(splits[[1]]==target.number)){
            splits<-split.cluster(vars=vars[-which(vars==target.number)], ME=ME)
          }
          if (length(splits[[2]])==1&any(splits[[2]]==target.number)){
            splits<-split.cluster(vars=vars[-which(vars==target.number)], ME=ME)
          }
        }
        if(any(splits[[1]]==target.number)){
          splits[[1]]<-splits[[1]][-which(splits[[1]]==target.number)]
        }
        if(any(splits[[2]]==target.number)){
          splits[[2]]<-splits[[2]][-which(splits[[2]]==target.number)]
        }
        if(any(variable.number%in%splits[[1]])){
          splits[[1]]<-splits[[1]][-which(splits[[1]]%in%variable.number)]
        }
        if(any(variable.number%in%splits[[2]])){
          splits[[2]]<-splits[[2]][-which(splits[[2]]%in%variable.number)]
        }
        #split 1
        j=1
        if(length(splits[[j]])==1){
          if(is.factor(data[,splits[[j]]])){
            #if it is already a factor
            data.scores1[,as.numeric(cl)]<-data[,splits[[j]]]
            disc.param1[[as.numeric(cl)]]<-NULL
            pca.param1[[as.numeric(cl)]]<-NULL
          }else{
            #if it is continuous: discretize
            if(discretize){
              disc<-discretize.dens(data[,splits[[j]]], cluster=F, seed=seed)
              data.scores1[,as.numeric(cl)]<-as.factor(disc$discretized)
              disc.param1[[as.numeric(cl)]]<-disc$levels
              pca.param1[[as.numeric(cl)]]<-NULL
            } else {
              data.scores1[,as.numeric(cl)]<-data[,splits[[j]]]
            }
          }
        }else{
          names<-colnames(data)[splits[[j]]]
          if(any(c(target, separate)%in%names)){
            names<-names[-which(names%in%c(target, separate))]
          }
          pc<-PCAmix.groupbn(X.quanti, X.quali, names, seed=seed)
          pca.param1[[as.numeric(cl)]]<-pc
          #keep column from before
          if(discretize){
            disc<-discretize.dens(pc$scores[,1], cluster=T, seed=seed)
            data.scores1[,as.numeric(cl)]<-as.factor(disc$discretized)
            disc.param1[[as.numeric(cl)]]<-disc$levels
          } else {
            data.scores1[,as.numeric(cl)]<-pc$scores[,1]
          }
        }

        #second split
        j=2
        if(length(splits[[j]])==1){
          if(is.factor(data[,splits[[j]]])){
            #if it is already a factor
            data.scores1[,n+1]<-data[,splits[[j]]]
            disc.param1[[n+1]]<-NULL
            pca.param1[[n+1]]<-NULL
            names(disc.param1)[n+1]<-paste("cl", n+1, sep="")
            names(pca.param1)[n+1]<-paste("cl", n+1, sep="")
          }else{
            #if it is continuous: discretize
            if(discretize){
              disc<-discretize.dens(data[,splits[[j]]], cluster=F, seed=seed)
              data.scores1[,n+1]<-as.factor(disc$discretized)
              disc.param1[[n+1]]<-disc$levels
              pca.param1[[n+1]]<-NULL
              names(disc.param1)[n+1]<-paste("cl", n+1, sep="")
              names(pca.param1)[n+1]<-paste("cl", n+1, sep="")
            } else {
              data.scores1[,n+1]<-data[,splits[[j]]]
            }
          }
        }else{
          names<-colnames(data)[splits[[j]]]
          if(any(c(target, separate)%in%names)){
            names<-names[-which(names%in%c(target, separate))]
          }
          pc<-PCAmix.groupbn(X.quanti, X.quali, names, seed=seed)
          pca.param1[[n+1]]<-pc
          #add column
          if(discretize){
            disc<-discretize.dens(pc$scores[,1], cluster=T, seed=seed)
            data.scores1[,n+1]<-as.factor(disc$discretized)
            disc.param1[[n+1]]<-disc$levels
            names(disc.param1)[n+1]<-paste("cl", n+1, sep="")
            names(pca.param1)[n+1]<-paste("cl", n+1, sep="")
          } else{
            data.scores1[,n+1]<-pc$scores[,1]
          }
        }
        cluster1[splits[[j]]]<-n+1
        colnames(data.scores1)[dim(data.scores1)[2]]<-paste("cl", n+1, sep="")

        #learn network
        start=bnlearn::empty.graph(nodes=colnames(data.scores1))
        arcs(start)<-bnlearn::arcs(net)[-which(arcs(net)[,1]==node|arcs(net)[,2]==node),]
        start<-bnlearn::cextend(start)
        blacklist<-NULL
        if(!is.null(separate)&any(separate.as.roots)){
          blacklist<-groupbn.build.blacklist(data.scores1[complete.cases(data.scores1),], separate[separate.as.roots])
        }
        netw1<-network(data.scores1[complete.cases(data.scores1),], R=R, start=start, seed=seed, debug=F, blacklist=blacklist, arc.thresh=arc.thresh, struct.alg=struct.alg, boot=boot)
        net1<-netw1[[1]]
        arc.confid1<-netw1[[2]]
        #Scoring
        result1<-validation(net1, data.scores1[complete.cases(data.scores1),], target, seed=seed, debug=F)

        return(list(data.scores1, cluster1, net1, arc.confid1, result1, pca.param1, disc.param1))
      }
      nbh<-lapply(listnod, eval.neighbourhood)
      data.scores.help<- sapply(nbh, `[`, 1)
      cluster.help<-sapply(nbh, `[`, 2)
      net.help<-sapply(nbh, `[`, 3)
      arc.confid.help<-sapply(nbh, `[`, 4)
      result.help<-sapply(nbh, `[`, 5)
      pca.param.help<-sapply(nbh, `[`, 6)
      disc.param.help<-sapply(nbh, `[`, 7)
      #compare prediction of each network and keep the best
      if(debug){message("Scores: ", paste(round(unlist(lapply(result.help, `[`, 1)),5), collapse=" "))}
      if(debug) {message("Score before: ", round(result,5), ", threshold: ", round(attr(result, "error.th"),5))}

      #minimum better than before?
      if(suppressWarnings(min(unlist(lapply(result.help, `[`, 1)), na.rm=T))==Inf|(suppressWarnings(min(unlist(lapply(result.help, `[`, 1)), na.rm=T))>attr(result, "error.th"))){
        if(debug){message("no further improvement. Stop.")}
        cluster[target]<-k+1
        break
      }

      if (suppressWarnings(min(unlist(lapply(result.help, `[`, 1)), na.rm=T))<=attr(result, "error.th")&l==1){
        choose<-1
        if(debug){message("Chosen to split: ",listnod[choose]," (",choose, ")")}else{message("Split: ", listnod[choose])}

        data.scores<-data.scores.help[[1]]
        cluster<-cluster.help[[1]]
        net<-net.help[[1]]
        arc.confid<-arc.confid.help[[1]]
        disc.param<-disc.param.help[[1]]
        pca.param<-pca.param.help[[1]]
        result.list<-lappend(result.list, result.help[[choose]])
        result<-result.help[[1]]
        n<-n+1
      }else{
        choose<-which.min(unlist(lapply(result.help, `[`, 1)))
        if(debug){message("Chosen to split: ",listnod[choose]," (",choose, ")")}else{message("Split: ", listnod[choose])}
        data.scores<-data.scores.help[[choose]]
        cluster<-cluster.help[[choose]]
        net<-net.help[[choose]]
        arc.confid<-arc.confid.help[[choose]]
        result.list<-lappend(result.list, result.help[[choose]])
        result<-result.help[[choose]]
        disc.param<-disc.param.help[[choose]]
        pca.param<-pca.param.help[[choose]]
        n<-n+1
      }
      #close step loop
    }
    cluster[target]<-k+1
    #for (i in 1:length(separate)){
    #  cluster[separate[i]]<-k+1+i
    #}
    res$bn<-net
    res$arc.confid<-arc.confid
    res$group.data<-data.scores
    res$grouping<-cluster
    res$score<-result
    res$pca.param<-pca.param
    res$disc.param<-disc.param
    if (target%in%colnames(X.quanti)|discretize==FALSE){
      res$fit<-bnlearn::bn.fit(bnlearn::cextend(net), data.scores[stats::complete.cases(data.scores),])
    } else if (target%in%colnames(X.quali)){
      res$fit<-bnlearn::bn.fit(bnlearn::cextend(net), data.scores[stats::complete.cases(data.scores),], method="bayes")
    }
    res$score<-result

    restart.list[[re+2]]<-res
    message("Run ", re+1, " of ",restart+1," done.")
    message(step-1, " refinement step(s) done.")
    message(".......................................................")
    }
  end_time <- Sys.time()
  duration=end_time-start_time
  if(debug){message("Calculation Time: ", round(duration[[1]], 2),  units(duration))}
  if(return.all){
    return(restart.list)
  } else {
    idx<-which.min(unlist(comb.output.scores(restart.list)))
    return(restart.list[[idx]])
  }
}


groupbn_refine_manually<-function(res, hierarchy, refine, arc.thresh=NULL, R=100, debug=FALSE, seed=NULL){
  #hierarchy: cluster object from ClustOfVar
  #res: groupbn object

  #timing
  start_time <- Sys.time()

  ######Initialization#####
  if (!is.null(seed)) {set.seed(seed)}
  k<-res$k
  target<-res$target
  separate<-res$separate
  X.quali<-res$X.quali
  X.quanti<-res$X.quanti
  net<-res$bn
  arc.confid<-res$arc.confid
  cluster<-res$grouping
  data.scores<-res$group.data
  result<-res$score
  target.number<-which(colnames(data)==target)
  variable.number<-which(colnames(data)%in%separate)
  pca.param<-res$pca.param
  disc.param<-res$disc.param
  n<-max(cluster)
  struct.alg<-attributes(res$bn)$struct.alg
  boot<-attributes(res$bn)$boot
  data<-cbind(X.quanti, X.quali)
  if(!is.null(separate)){
    separate.as.roots<-attr(res$separate, "separate.as.roots")
  } else {
    separate.as.roots<-FALSE
  }

  ##############################
  #########Refinement###########
  ##############################

  #get merge matrix from cluster object
  ME<-hierarchy$merge
  listnod<-refine

  if(any(listnod%in%c(separate,target))){
    stop("Cluster cannot be split")
  }
  cl.nrs<-substr(listnod, rep(3, length(listnod)), nchar(listnod))
  remove<-c()
  for (nod in 1:length(cl.nrs)){
    if(length(which(cluster==cl.nrs[nod]))==1){
      remove<-c(remove, nod)
    }
  }
  if(length(remove)){
    listnod<-listnod[-remove]
  }
  pert<-listnod
  pert.idx<-1

  if(debug){message("Split is introduced")}
  #loop through number of perturbations and introduce them
    cl<-substr(pert[pert.idx], 3, nchar(pert[pert.idx]))

    #set number of clusters
    n<-max(cluster)

    pca.param<-lappend(pca.param, NULL)
    disc.param<-lappend(disc.param, NULL)
    #variable names
    vars<-names(cluster[which(cluster==cl)])
    if(debug){message("Variables: ", paste(vars, collapse=" "))}

    #transform to numbers
    vars<-which(colnames(data)%in%vars)
    #only one variable: skip, not splittable
    if (length(vars)==1){
      if(debug){message("not splittable, skip")}
    }


    else if(length(vars)==2&any(vars%in%c(separate,target))){
      #two variables, one is sex/age/target: skip
      if(debug){message("not splittable, skip")}
    }

    #otherwise determine splitting from merge matrix
    splits<-split.cluster(vars=vars, ME=ME)
    #first split
    j=1
    if(length(which(colnames(data)[splits[[j]]]%in%c(separate,target)))){
      splits[[j]]<-splits[[j]][-which(colnames(data)[splits[[j]]]%in%c(separate,target))]
    }
    #if only one variable in cluster
    if(length(splits[[j]])==1){
      if(is.factor(data[,splits[[j]]])){
        #if it is already a factor
        data.scores[,as.numeric(cl)]<-data[,splits[[j]]]
        pca.param[[as.numeric(cl)]]<-NULL
        disc.param[[as.numeric(cl)]]<-NULL
      }
      else{
        #if it is continuous: discretize
        disc<-discretize.dens(data[,splits[[j]]], cluster=F, seed=seed)
        data.scores[,as.numeric(cl)]<-as.factor(disc$discretized)
        pca.param[[as.numeric(cl)]]<-NULL
        disc.param[[as.numeric(cl)]]<-disc$levels
      }
    }else{
      names<-colnames(data)[splits[[j]]]
      pc<-PCAmix.groupbn(X.quanti, X.quali, names, seed=seed)
      pca.param[[as.numeric(cl)]]<-pc
      disc<-discretize.dens(pc$scores[,1], cluster=T, seed=seed)
      disc.param[[as.numeric(cl)]]<-disc$levels
      #substitute column from before
      data.scores[,as.numeric(cl)]<-as.factor(disc$discretized)
    }
    #second split
    j=2
    if(length(which(colnames(data)[splits[[j]]]%in%c(separate,target)))){
      splits[[j]]<-splits[[j]][-which(colnames(data)[splits[[j]]]%in%c(separate,target))]
    }
    if(length(splits[[j]])==1){
      if(is.factor(data[,splits[[j]]])){
        #if it is already a factor
        data.scores[,n+1]<-data[,splits[[j]]]
      }
      else{
        disc<-discretize.dens(data[,splits[[j]]], cluster=F, seed=seed)
        data.scores[,n+1]<-as.factor(disc$discretized)
        disc.param[[n+1]]<-disc$levels
      }
    }
    else{
      names<-colnames(data)[splits[[j]]]
      pc<-PCAmix.groupbn(X.quanti, X.quali, names, seed=seed)
      pca.param[[n+1]]<-pc
      #add column
      disc<-discretize.dens(pc$scores[,1], cluster=T, seed=seed)
      data.scores[,n+1]<-as.factor(disc$discretized)
      disc.param[[n+1]]<-disc$levels
    }
    colnames(data.scores)[dim(data.scores)[2]]<-paste("cl", n+1, sep="")
    cluster[splits[[j]]]<-n+1
    n<-max(cluster)
  #learn network for perturbed cluster structure
  blacklist<-NULL
  if(!is.null(separate)&any(separate.as.roots)){
    blacklist<-groupbn.build.blacklist(data.scores[complete.cases(data.scores),], separate[separate.as.roots])
  }
  netw<-network(data.scores[complete.cases(data.scores),], R=R, seed=seed, debug=F, blacklist=blacklist , arc.thresh=arc.thresh, struct.alg=struct.alg, boot=boot)
  net<-netw[[1]]
  arc.confid<-netw[[2]]
  result<-validation(net, data.scores, target, seed=seed, debug=F)
  if(debug) {message("Score: ", round(result,2), ", threshold: ", round(attr(result, "error.th"),2))}
  res$bn<-net
  res$arc.confid<-arc.confid
  res$group.data<-data.scores
  res$grouping<-cluster
  res$score<-result
  res$pca.param<-pca.param
  res$disc.param<-disc.param
  if (target%in%colnames(X.quanti)){
    res$fit<-bnlearn::bn.fit(bnlearn::cextend(net), data.scores[stats::complete.cases(data.scores),])
  } else if (target%in%colnames(X.quali)){
    res$fit<-bnlearn::bn.fit(bnlearn::cextend(net), data.scores[stats::complete.cases(data.scores),], method="bayes")
  }
  res$score<-result
  end_time <- Sys.time()
  duration=end_time-start_time
  if(debug){message("Calculation Time: ", round(duration[[1]], 2),  units(duration))}
  return(res)
}
