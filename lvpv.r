# Code for pvclust approximately unbiased p-value realisation wrt Newman's leading eigenvector algorithm for detection of integration and modularity in complex systems

# an implementation of a multiscale bootstrap procedure designed to estimate the p-values in the modular structure analysis with the leading eigenvector algorithm

# by Oksana Vertsimakha and Igor Dzeverin, 2020-09-10

# References
# Csardi, G., Nepusz, T. (2006). The igraph software package for complex network research. InterJournal, Complex Systems 1695. http://igraph.org
# Melo D., Garcia G., Hubbe A., Assis A. P., Marroig G. (2016). EvolQG – An R package for evolutionary quantitative genetics [version 3; peer review: 2 approved, 1 approved with reservations]. F1000Research, 4(925), doi: 10.12688/f1000research.7082.3
# Newman, M. E. J. (2006). Modularity and community structure in networks. Proceedings of the National Academy of Sciences of the USA, 103, 8577-8582, doi: 10.1073/pnas.0601602103
# Shimodaira, H. (2004) Approximately unbiased tests of regions using multistep-multiscale bootstrap resampling. The Annals of Statistics, 32, 2616-2641. doi: 10.1214/009053604000000823


# creates graph from data's cor matrix, returns clusters, merges & modularity;
# (function requires igraph package)

# method and use same as cor() function options;
# gauss allows for prior gaussian tranformation of the correlation matrix, default = FALSE;
# q.cut cuts off correlations below specified quantile (counted w/o diagonal), default = 0.5;
# diag allows for loops in the graph, default = FALSE;
# square allows for squared cor matrix usage in case of negative correlations, default = TRUE;
# warnings notifies of squared cor matrix usage in case of negative correlations, default = FALSE


cor.graph = function(data, 
                     method = "pearson",
                     use = "everything",
                     gauss = FALSE,
                     q.cut = 0.5,
                     diag = FALSE,
                     square = TRUE,
                     warnings = FALSE){
  suppressMessages(require(igraph))
  C = cor(data, method = method, use = use)
  if(any(C<0)&square){
    C = C^2
    if(warnings) warning("Some correlations are negative. Using squared correlations.")
  }
  if(gauss|!square) C = exp(-(1-C)^2/(2*sd(C)^2))
  if((0<q.cut)&(q.cut<1)) diag(C) = 0; C[C^2 < quantile(C^2, q.cut)] = 0
  
  g = graph.adjacency(C, weighted = TRUE, diag = diag, mode = 'undirected')
  res = cluster_leading_eigen(g, options = list(maxiter = 1000000))
  orig.cl = res$membership; orig.merges = res$merges
  names(orig.cl)<-colnames(data)
  orig.mod = res$modularity
  
  return(list(orig.cl, orig.merges, orig.mod))
}


# (inner function)
# inputs: membership and merges as obtained by cluster_leading_eigen igraph function.
# returns list corresponding all obtained clusters:
# [[1]] pattern - 0/1 strings (e.g. "110100" = cluster of elements 1,2 and 4 out of 1-6);
# [[2]] member - list of indexes of cluster elements (e.g. c(1,2,4) for the previous example)

SplitMerges = function(c.membership, c.merges){
  n <- length(c.membership)
  m <- nrow(c.merges)
  K <- max(c.membership)
  B <- list()
  #print(m)
  
  for(i in 1:K) B[[i]] <- which(c.membership == i)
  
  if(m > 0){
    for(j in 1:m){
      m1 = c.merges[j,1]; m2 = c.merges[j,2]
      B[[j+K]] = sort(c(B[[m1]], B[[m2]]))
    }
  }
  
  CC <- matrix(rep(0,n*(K+m)), ncol=n)
  for(i in 1:length(B)){
    bi <- B[[i]]
    mm <- length(bi)
    for(j in 1:mm)
      CC[i,bi[j]] <- 1
  }
  
  split <- list(pattern=apply(CC, 1, paste, collapse=""), member=B)
  
  return(split)
}

# bootstrap for clusters and modularity
# r = bootstrap scale parameter: nboot - number of bootstrap samples;
# store = store sample correlations;
# quiet hides progress messages
# mod.alpha sets levels of modularity bootstrap sample quanties (alpha/2, 1-alpha/2)
# returns list of counted clusters, parameters e.g. nboot, r, ...


PV.boot<-function(data, r, nboot = 1000,
                  method = "pearson", use = "everything", 
                  warnings = F, gauss = F, q.cut = 0, diag = F, square = square,
                  store = F, quiet = F,
                  mod.alpha = .05){
  
  n <- nrow(data)
  size<- round(n*r, digits=0)
  if(size == 0)
    stop("invalid scale parameter(r)")
  r <- size/n
  
  g.cor = cor.graph(data, method = method, use = use, 
                    gauss = gauss, q.cut = q.cut, 
                    diag = diag, square = square, warnings = warnings)
  pattern   <- SplitMerges(g.cor[[1]], g.cor[[2]])$pattern
  edges.cnt <- table(factor(pattern)) - table(factor(pattern))
  st <- list()
  rp <- as.character(round(r, digits=2)); if(r == 1) rp <- paste(rp,".0",sep="")
  mod = list()
  
  if(!quiet)
    cat(paste("Bootstrap (r = ", rp, ")... ", sep=""))
  
  na.flag = 0
  for(i in 1:nboot){
    smpl <- sample(1:n, size, replace=TRUE)
    boot.data = data[smpl,]
    
    x.cor = cor.graph(boot.data, method = method, use = use, gauss = gauss, square = square, q.cut = q.cut, diag = diag)
    mod[[i]] = x.cor[[3]]
    
    if(sum(is.na(x.cor)[[2]])==0){
      
      pattern.i <- SplitMerges(x.cor[[1]], x.cor[[2]])$pattern   # split
      edges.cnt <- edges.cnt + table(factor(pattern.i,  levels=pattern))
      
    } else na.flag <- na.flag + 1
    
    if(store)
      st[[i]] <- x.cor
  }
  if(!quiet)
    
    cat("Done.\n")
  
  if(na.flag == 1)
    warning(paste("unseccessful networks divisions are omitted in computation: r = ", r), call.=FALSE)
  
  
  cmod = unlist(mod)
  
  mod.mean = mean(cmod); mod.sd = sd(cmod)
  mod.q1 = quantile(cmod, mod.alpha/2); mod.q2 = quantile(cmod, 1 - mod.alpha/2)
  if(!square) gauss = TRUE
  
  boot <- list(edges.cnt = edges.cnt,
               mod.mean = mod.mean, mod.sd = mod.sd,
               mod.q1 = mod.q1, mod.q2 = mod.q2,
               mod.alpha = mod.alpha,
               method = method, use = use,
               q.cut = q.cut,
               diag = diag,
               gauss = gauss,
               square = square,
               nboot = nboot, 
               size=size, r=r, 
               na.flag = na.flag, store=st)
  #class(boot) <- "boot.clust"
  
  return(boot)
}


# (inner function)
# calculates BP and AU p-values estimates 
# bp = boostrap counts; r = vector of scaling parameters; nboot = number of samples

PV.au_bp.fit <- function(bp, r, nboot) {
  
  if(length(bp) != length(r))
    stop("bp and r should have the same length")
  nboot <- rep(nboot, length=length(bp))
  ## suitable parameters check
  min.use <- 3 
  eps <- 0.001
  bp.use <- bp > eps & bp < 1-eps
  
  p <- se <- c(0,0)
  names(p) <- names(se) <- c("au", "bp")
  coef <- c(0,0); names(coef) <- c("v", "c")
  
  a <- list(p=p, se=se, coef=coef, df=0, rss=0, pchi=0) 
  #class(a) <- "msfit"
  
  if(sum(bp.use) < min.use) {
    if(mean(unlist(bp)) < .5) a$p[] <- c(0, 0) else a$p[] <- c(1, 1)
    return(a)
  }
  
  
  
  bp <- bp[bp.use]; r <- r[bp.use]; nboot <- nboot[bp.use]
  zz <- -qnorm(bp)
  vv <- ((1-bp)*bp)/(dnorm(zz)^2*nboot)
  a$bp.use <- bp.use; a$r <- r; a$zz <- zz
  
  X   <- cbind(sqrt(r), 1/sqrt(r)); dimnames(X) <- list(NULL, c("v","c"))
  fit <- lsfit(X, zz, 1/vv, intercept=FALSE)
  a$coef <- coef <- fit$coef
  
  h.au <- c(1, -1); h.bp <- c(1, 1)
  z.au <- drop(h.au %*% coef); z.bp <- drop(h.bp %*% coef)
  p.au <- pnorm(-z.au); p.bp <- pnorm(-z.bp)
  
  a$p["au"] <- p.au; a$p["bp"] <- p.bp
  
  
  V <- solve(crossprod(X, X/vv))
  vz.au <- drop(h.au %*% V %*% h.au); vz.bp <- drop(h.bp %*% V %*% h.bp)
  
  a$se["au"] <- dnorm(z.au) * sqrt(vz.au)
  a$se["bp"] <- dnorm(z.bp) * sqrt(vz.bp)
  a$rss <- sum(fit$residual^2/vv)
  
  
  if((a$df <- sum(bp.use) - 2) > 0) {
    a$pchi <- pchisq(a$rss, lower.tail=FALSE, df=a$df)
  }
  else a$pchi <- 1.0
  
  return(a)
  
}


# finds AU and BP p-values estimates obtained from data w/ parameters defined in previous function;
# returns list of obtained p-values for each split & scale parameter, modularity bootstrap stats,... 
# ...input parameters, communities and modularity of original partition etc
# iseed sets seed for the bootstrap 

PV.complete<- function(data, nboot = 1000, 
                       r = seq(.5,1.4,by=.1),
                       method="pearson", use="everything",
                       warnings=FALSE, gauss=FALSE, square=TRUE, q.cut=0.5, diag=FALSE,
                       iseed=NULL, quiet=FALSE, store=FALSE,
                       mod.alpha = .05){
  
  if(!is.null(iseed))
    set.seed(seed = iseed)
  
  
  mboot <- lapply(r, PV.boot, data = data, 
                  method = method, use = use,
                  q.cut = q.cut,
                  gauss = gauss,
                  square = square,
                  diag = diag,
                  nboot = nboot,
                  mod.alpha = mod.alpha)
  
  ### merge results
  g.full = cor.graph(data, method = method, use = use, q.cut = q.cut, 
                     diag = diag, square = square, gauss = gauss, warnings = warnings)
  
  pattern <- SplitMerges(g.full[[1]], g.full[[2]])$pattern
  r <- unlist(lapply(mboot,"[[","r"))
  nboot <- unlist(lapply(mboot,"[[","nboot"))
  mod.mean = unlist(lapply(mboot,"[[","mod.mean"))
  mod.sd = unlist(lapply(mboot,"[[","mod.sd"))
  mod.q1 = unlist(lapply(mboot,"[[","mod.q1"))
  mod.q2 = unlist(lapply(mboot,"[[","mod.q2"))
  store <- lapply(mboot,"[[", "store"); names(store)=paste("r", r, sep = "")
  rl <- length(mboot); ne <- length(pattern)
  if(!square) gauss = TRUE

  
  ### bootstrap bp results
  edges.bp <- edges.cnt <- data.frame(matrix(rep(0,ne*rl), ncol=rl))
  row.names(edges.bp) <- pattern
  
  names(edges.cnt) <- paste("r", 1:rl, sep="")
  for(j in 1:rl) {
    edges.cnt[,j] <- as.vector(mboot[[j]]$edges.cnt) 
    edges.bp[,j]  <- edges.cnt[,j] / nboot[j]
  }
  
  ms.fitted <- lapply(as.list(1:ne),
                      function(x, edges.bp, r, nboot){
                        PV.au_bp.fit(as.vector(t(edges.bp[x,])), r, nboot)},
                      edges.bp, r, nboot)
  
  p    <- lapply(ms.fitted,"[[","p")
  
  se   <- lapply(ms.fitted,"[[","se")
  
  coef <- lapply(ms.fitted,"[[","coef")
  
  au    <- unlist(lapply(p,"[[","au"))
  
  bp    <- unlist(lapply(p,"[[","bp"))
  
  se.au <- unlist(lapply(se,"[[","au"))
  
  se.bp <- unlist(lapply(se,"[[","bp"))
  
  v     <- unlist(lapply(coef,"[[","v"))
  
  cc    <- unlist(lapply(coef,"[[","c"))
  
  pchi  <- unlist(lapply(ms.fitted,"[[","pchi"))
  
  edges.pv <- data.frame(au=au, bp=bp, se.au=se.au, se.bp=se.bp,
                         v=v, c=cc, pchi=pchi)

  
  row.names(edges.pv) <- row.names(edges.cnt) <-pattern
  
  
  ### modularity stats
  
  modularity.bp = t(data.frame(r=r, mean=mod.mean, sd=mod.sd, q1=mod.q1, q2=mod.q2))
  colnames(modularity.bp) = rep("", length(r))
  
  mod.stats = list()

  if(1 %in% r){
    
    mod.stats[[1]] = mod.mean[r==1]
    
    mod.stats[[2]] = mod.sd[r==1]
    
    mod.stats[[3]] = mod.q1[r==1]
    
    mod.stats[[4]] = mod.q2[r==1]
  }
  else{
    nboot1 = nboot(data, r=1, nboot = nboot,
                   method = method, use = use, 
                   warnings = warnings, gauss = gauss,
                   square = square,
                   q.cut = q.cut, diag = diag,
                   store = F, quiet= T,
                   mod.alpha = mod.alpha)
    
    mod.stats[[1]] = nboot1$mod.mean
    
    mod.stats[[2]] = nboot1$mod.sd
    
    mod.stats[[3]] = nboot1$mod.q1
    
    mod.stats[[4]] = nboot1$mod.q2
  }
  
  mod.stats = data.frame(t(unlist(mod.stats))); colnames(mod.stats) = c("mean", "sd", "q1", "q2")

  
  pv.result = list(
    membership = g.full[[1]],
    modularity = g.full[[3]],
    modularity.bp = modularity.bp,
    mod.stats = mod.stats,
    merges = g.full[[2]],
    #hclust = make_dendro(g.full[[1]], g.full[[2]]),
    count = edges.cnt,
    edges.bp = edges.bp,
    nboot = nboot, r=r,
    edges = edges.pv,
    method = list(cor.method = method, use = use, gauss = gauss, q.cut = q.cut, 
                  diag = diag, square = square, mod.alpha = mod.alpha, seed = iseed),
    
    store=store
  )
  return(pv.result)
}



# (inner function)
# transforms merge matrix for the dendrogram parameters purposes
# returns data.frame with x and y coordinates of the nodes

lv.hc2axes <- function(x){
  A <- x$merge 
  m <- nrow(A); K = length(x$labels)
  t = t(x$merge)
  #x.axis <- -t[t<0]
  x.axis = 1:K
  
  y.axis <- c(rep(1:m, each = 2)[as.vector(t<0)], x$height)
  
  x.tmp  <- rep(0,2)
  zz <- match(1:length(x$order), x$order)
  
  for(i in 1:m) {
    ai <- A[i,1]
    
    if(ai < 0)
      x.tmp[1] <- zz[-ai]
    else
      x.tmp[1] <- x.axis[ai+K]
    
    ai <- A[i,2]
    
    if(ai < 0)
      x.tmp[2] <- zz[-ai]
    else
      x.tmp[2] <- x.axis[ai+K]
    x.axis[i+K] <- mean(x.tmp)
    
  }
  return(data.frame(x.axis=x.axis, y.axis=y.axis))
  
}

# (inner function)
# creates dendrogram from the membership vector and merges matrix (of igraph type)
# returns hclust object or empty list if only one cluster was detected

make_dendro = function(xC, xM){
  
  K = length(unique(xC))
  m = nrow(xM)
  
  a = list()
  
  if(m==0|m==1) cat("Less then three clusters, cannot make a dendrogram")
  else if(m > 1) {
    new_xM = apply(xM, 2, function(z){(z>K)*(z-K)+(z<=K)*(-z)})
    a$merge = new_xM
    a$height = 1:m
    
    a$order = -t(new_xM)[t(new_xM)<0]
    a$labels <- sapply(1:K, function(z) paste("C",z,sep = " "))
    class(a) <- "hclust"      
  }
  
  return(a)
}


# creates dendrogram according to leading eigenvector communities
# x = PV.complete function output, plot parameters; 
# additional parameters (...) for the plot function

PV.dendro <-function(x, print.pv=TRUE, float=0.01,
                     col.pv=c(au=2, bp=3), cex.pv=0.8, font.pv=NULL,
                     col=NULL, cex=NULL, font=NULL, lty=NULL, lwd=NULL,
                     main=NULL, sub=NULL, xlab=NULL, ylab = NULL, ...){
  
  xC = x$membership; xM = x$merges
  
  a = make_dendro(xC, xM)
  if(is.null(main))
    main="Cluster dendrogram with p-values (%)"
  
  if(is.null(xlab))
    xlab = ""
  
  if(is.null(ylab))
    ylab = ""
  if(length(xM)>2) {
    plot(a, main=main, sub=sub, 
         xlab=xlab, ylab=ylab, col=col,
         cex=cex, font=font, lty=lty, lwd=lwd, ...)
  }
  
  
  if(isTRUE(print.pv)) print.pv   <- c("au", "bp")
  col.text <- col.pv[print.pv]
  
  #text(x, col=col.text, cex=cex.pv, font=font.pv, float=float, print.num=print.num)
  
}


# adds estimated BP and AU p-values next to the splits to the dendrogram by PV.dendro
# x = PV.complete function output, text parameters
# to hide au or bp values setrespective colorto transparent


PV.text <- function(x, col=c(au=2, bp=3), 
                    print.num=TRUE, float=0.01, cex=NULL, font=NULL, ...)
{
  
  xC = x$membership; xM = x$merges
  a = make_dendro(xC, xM)
  
  if(length(col) == 1) col = rep(col, 2)
  names(col) <- c("au", "bp")
  
  axes <- lv.hc2axes(a)
  usr  <- par()$usr; wid <- usr[4] - usr[3]
  
  
  num_str <- lapply(
    x$edges[seq_len(which(names(x$edges) == "bp"))],
    function(p) round(p * 100))
  
  
  # change the last element to the name of p-value
  for(i in names(num_str)) {
    num_str[[i]][length(num_str[[i]])] <- i
  }
  
  K = length(a$labels)
  m = nrow(a$merge)
  
  
  
  range <- seq_len(min(3, length(col)))
  pos <- c(2, 4, 1)
  y_off <- float * wid
  y_offset<- c(rep(-y_off*3, K), rep(y_off*2, m))
  
  #print(axes[,2]); print(y_offset)
  ind = c(a$order, K+(1:m))
  for(i in range) {
    
    name <- names(col)[i]
    text(x=axes[,1], y=axes[,2] + y_offset, num_str[[name]][ind],
         col=col[name], pos=pos[i], offset=.2, cex=cex, font=font, ...)
    
  }
}



# creates summary: cluster elements and p-value information, used method info
# x = PV.complete output list;
# which specifies cluster numbers according to the dendoram, e.g. for C1 enter which = 1
# if not specified all clusters are shown; 
# modularity shows modularity bootstrap results 
# digits sets round parameter

PV.summary<-function(x, which = NULL, digits = 3, modularity = TRUE){
  d = TRUE
  n = nrow(x$edges.bp); K = length(unique(x$membership->cl))
  
  if(is.null(which)) which = 1:(n-1)
  mtch = intersect(which, 1:K)
  
  if(n == 1) {print("One cluster detected"); d = F}
  
  if(length(mtch)==0) {print(paste("Wrong cluster indexes,", K, "clusters found", sep = " ")); d = F}
  
  if(d){
    names = names(x$membership)
    if(is.null(names)) names = 1:length(x$membership) 
    
    cat("Clusters:\n\n")
    for(j in mtch){
      
      cat(paste("C", j, ": ", sep = ""))
      
      cat(paste(names[cl==j], sep = ", "), labels = NULL)
      
      cat("\n")
    }
    cat("Estimates on edges:\n\n")
    
    print(round(x$edges[mtch,], digits=digits))
    
    cat("\n")
    
    cat(paste("Correlation method: ", x$method[[1]],', use = ', x$method[[2]], "\n", sep = ""))
    if(!x$method[[6]]) cat('Negative correlations used',"\n")
    if(x$method[[3]]) cat('Gauss transformation used',"\n")
    if((is.numeric(x$method[[4]]))&(x$method[[4]]>0)) cat('correlations cut at', 
                                                          paste(round(100*x$method[[4]],2),"%",sep=""), 
                                                          'quantile', "\n")
    if(x$method[[5]]) cat('diag = TRUE, graph built with loops', "\n")
  }
  if(modularity) cat(paste('Modularity: ', round(x$modularity, 4), "\n", 
                           'Bootstrap estimate: ', round(x$mod.stats["mean"],4), "\n",
                           'Bootstrap sd: ', round(x$mod.stats["sd"],4), "\n", 
                           'Bootstrap ', paste(round(100*x$method$mod.alpha/2,2),"%",sep=""), ' quantile: ',
                           round(x$mod.stats["q1"], 4), "\n", 
                           'Bootstrap ', paste(round(100*(1-x$method$mod.alpha/2),2),"%",sep=""), ' quantile: ',
                           round(x$mod.stats["q2"], 4), sep = ""))
}



# creates rectangles around clusters of the dendrogram that satisfy condition on p-value estimates
# x = PV.complete outut, other parameters the same as of pvrect function

PV.rectangle <- function(x, alpha = 0.95, pv = "au", type = "geq", max.only = T,
                          border = NULL,...){
  
  m = nrow(x$merges)
  K = max(x$membership)
  
  order = t(x$merges)[t(x$merges)<=K]
  member = SplitMerges(1:K, x$merges)$member
  #a = make_dendro(x$membership, x$merges)
  
  height = rep(1:m,each=2)[t(x$merges)<=K][order(order)]
  yt_height = c(height, 1:(m-1))
  
  usr = par("usr")
  xwd = usr[2]-usr[1]
  ywd = usr[4]-usr[3]
  #ywd = m
  cin = par()$cin
  
  if(is.null(border)) border = c(au = 2, bp = 3)[pv]
  
  ht <- c()
  j <- 1
  
  if(is.na(pm <- pmatch(type, c("geq", "leq", "gt", "lt"))))
    stop("Invalid type argument: see help(pvrect)")
  
  
  for(i in (K+m-1):1)
  {
    if     (pm==1) wh <- (x$edges[i,pv] >= alpha) # Greater than or EQuals
    else if(pm==2) wh <- (x$edges[i,pv] <= alpha) # Lower than or EQuals
    else if(pm==3) wh <- (x$edges[i,pv] >  alpha) # Greater Than
    else if(pm==4) wh <- (x$edges[i,pv] >  alpha) # Lower Than
    
    if(wh)
    {
      mi <- member[[i]]
      ma <- match(mi, order)
      
      if(max.only == FALSE || (max.only && sum(match(ma, ht, nomatch=0)) == 0))
      {
        xl <- min(ma)
        xr <- max(ma)
        yt <- yt_height[i]
        yb <- usr[3]
        
        mx <- xwd / length(member) / 4
        my <- ywd / 200
        
        rect(xl - mx, yb + 3*my, xr + mx, yt + my, 
             border=border, shade=NULL, ...)
        j <- j + 1
      }
      ht <- c(ht, ma)
      
    }
  }
}


### creates via igraph a binary graph corresponding to detected module structure
#input= PV.complete ouput; members = create nodes for membersof each final cluser, 
# pv = add AU and BP p-values text labels, m.labels = add names of module memebers, 
# which = to which modules does m.labels apply (all by default, else give cluster number vector)

PV.graph = function(x, members = TRUE, 
                    pv = TRUE, 
                    m.labels = TRUE,
                    which = NULL,
                    
                    vertex.shape = "rectangle",
                    vertex.color = "white",
                    vertex.frame.color = "black",
                    vertex.size = 16,
                    vertex.label.color = "black",
                    vertex.label.cex = .9,
                    vertex.label.color2 = NULL,
                    vertex.label.cex2 = .9, ...
                    
)
{
  suppressMessages(require(igraph))
  xC = x$membership; xM = x$merges
  K = max(xC); m = nrow(xM)
  if(m==0) stop("No clusters, no graph")
  
  ##add au bp test here
  
  G<-graph.tree(n = 3, children = 2, mode = "out")
  V(G)$name = c("", paste(xM[m,]))
  if(m>=2){
    for(j in (m-1):1){
      G <- G + vertices(paste(xM[j,]))
      G <- G + edge(paste(K+j), paste(xM[j,1]))
      G <- G + edge(paste(K+j), paste(xM[j,2]))
      
    }
  }
  
  nn = 0
  index = c(1, rep(2,2*m), rep(3, nn))
  v.labels = V(G)$name
  
  if(members){
    
    if(is.null(which)) {which = 1:K}
    include = xC %in% which; xCin = xC[include]
    
    names = names(xCin); nn = sum(include)
    
    for(k in 1:nn){
      G = G+vertices(paste(names[k]))+edge(paste(xCin[k]), paste(names[k]))
    }
    
    index = c(1, rep(2,2*m), rep(3, nn))
    v.labels = V(G)$name
    
    if(m.labels) v.labels[(2*m+1)+1:nn] = names
    else v.labels[(2*m+1)+1:nn] = rep("", nn)
    if(is.null(vertex.label.color2)) vertex.label.color2 = vertex.label.color
    if(is.null(vertex.label.cex2)) vertex.label.cex2 = vertex.label.cex
  }
  
  index[1] = 3
  if(pv)
  {
    index[1] = 1
    o = c(t(x$merges[m:1,]))
    v.labels[1] = "AU|BP"
    v.labels[2:(2*m+1)] = paste(
      round(100*x$edges$au[o], 0),
      round(100*x$edges$bp[o], 0),
      sep = "|"
    )
    
  }
  
  v.color = c("white", vertex.color, "black")[index]
  v.size = c(30, vertex.size, 4)[index]
  v.shape = c("none", vertex.shape, "circle")[index]
  v.frame = c("none", vertex.frame.color, "black")[index]
  
  l.color = c("black", vertex.label.color, vertex.label.color2)[index]
  v.label.cex = c(1, vertex.label.cex, vertex.label.cex2)[index]
  label.dist = c(0,0,1)[index]
  
  co <- layout.reingold.tilford(G, mode = "out", root="") 
  vertex_attr(G) = list(
    shape = v.shape,
    color = v.color,
    size = v.size,
    frame.color = v.frame,
    label = v.labels,
    label.cex = v.label.cex,
    label.dist = label.dist,
    label.color = l.color,
    label.degree = c(0, pi/2)[label.dist+1]
  )
  
  e = length(E(G))
  edge_attr(G) = list(
    size = rep(.8, e),
    arrow.size = rep(.2, e)
  )
  graph_attr(G, "layout") = co
  
  
  return(G)
}



# draws graph and highlights (option highlight = TRUE, default) 
# nodes that satisfy p-values conditions (see pvrect options)
# input: x =PV.complete output, ... = PV.graph visual options
PV.graph.highlight = function(x, highlight = TRUE, 
                              alpha=0.95, with="au", type="geq", col = "pink", ...){
  m = nrow(x$merges); K = max(x$membership)
  ht <- c()
  j <- 1
  g = PV.graph(x, ...)
  
  if(is.na(pm <- pmatch(type, c("geq", "leq", "gt", "lt"))))
    stop("Invalid type argument: see help(pvrect)")
  
  if(with == "au") w = 1 else if(with == "bp") w = 2
  else stop("Invalid with argument: use au or bp")
  
  order = c(t(x$merges[m:1,]))
  for(i in 1:(K+m-1))
  {
    if     (pm==1) wh <- (x$edges[order[i],w] >= alpha) # Greater than or EQuals
    else if(pm==2) wh <- (x$edges[order[i],w] <= alpha) # Lower than or EQuals
    else if(pm==3) wh <- (x$edges[order[i],w] >  alpha) # Greater Than
    else if(pm==4) wh <- (x$edges[order[i],w] >  alpha) # Lower Than
    
    if(wh) ht = c(ht, i)
  }
  
  
  vc = vertex.attributes(g)$color
  if(highlight) vc[ht+1] = col else vc[ht+1] = "white"
  plot(g,
       vertex.color = vc,
       layout = graph_attr(g, "layout"))
  
}



# ### TEST RUN
# source("lvpv.r")
# 
# ## create TestData dataset
# library(evolqg)
# library(MASS)
# 
# set.seed(254)
# m = rnorm(8, mean=10, sd=4)
# s = m*0.1 + rnorm(8, mean = 0, sd = 1)
# V = s^2
# VCV = RandomMatrix(8, variance=V, LKJ=FALSE)
# TestData = mvrnorm(n=240, mu=m, Sigma=VCV)
# colnames(TestData) = c("v1", "v2", "v3", "v4", "v5", "v6", "v7", "v8")
# 
# ## analyze TestData using Newman's leading eigenvector method
# pvTD = PV.complete(TestData, iseed=767)
# PV.dendro(pvTD); PV.text(pvTD) 
# plot(PV.graph(pvTD))
# PV.graph.highlight(pvTD, col="lightblue")
# PV.summary(pvTD)
# 
# ## compare to the results of cluster analysis and Sammon's mapping
# 
# library(pvclust)
# hclTD = pvclust(TestData, method.hclust = "average")
# plot(hclTD)
# pvrect(hclTD)
# 
# sTD = sammon(as.dist(1-cor(TestData)), k = 2)
# plot(sTD$points[,1], sTD$points[,2], bg="springgreen", pch=21, cex=2, xlab="Axis 1", ylab="Axis 2"); text(sTD$points[,1], sTD$points[,2], label=rownames(sTD$points), pos=4)
# 
# # the results of cluster analysis and Sammon's method differ from the results of Newman's leading eigenvector method because v5 negatively correlates with other variables
# 
# # to obtain the results of cluster analysis and Sammon's method similar to the results of Newman's leading eigenvector method, use
# TestData5 = TestData
# TestData5[,5] = TestData[,5]*(-1)
# 
# hclTD5 = pvclust(TestData5, method.hclust = "average")
# plot(hclTD5)
# pvrect(hclTD5)
# 
# sTD5 = sammon(as.dist(1-cor(TestData5)), k = 2)
# plot(sTD5$points[,1], sTD5$points[,2], bg="springgreen", pch=21, cex=2, xlab="Axis 1", ylab="Axis 2"); text(sTD5$points[,1], sTD5$points[,2], label=rownames(sTD5$points), pos=4)
# 
# # to obtain the results of Newman's leading eigenvector method similar to the results of cluster analysis and Sammon's method, use
# pvTD.n = PV.complete(TestData, iseed=767, square=FALSE)
# PV.dendro(pvTD.n); PV.text(pvTD.n)
# PV.graph.highlight(pvTD.n, col="lightblue")
# PV.summary(pvTD.n)
# 
# ## compare to the modularity detection algorithm realized in evolqg
# 
# # use evolqg
# lmod = LModularity(cor(TestData))
# lmod.partition = integer(dim(TestData)[2])
# for (j in 1:dim(TestData)[2]) lmod.partition[j] = which(lmod$Modularity_hypothesis[j,]==1)
# names(lmod.partition) = colnames(TestData)
# lmod$LModularity
# lmod$Modularity_hypothesis
# lmod.partition
# 
# # use lvpv, specify the same calculation algorithm as in evolqg 
# pvTD.diag = PV.complete(TestData, iseed=767, q.cut=0, diag=TRUE)
# PV.dendro(pvTD.diag); PV.text(pvTD.diag) 
# PV.graph.highlight(pvTD.diag, col="lightblue")
# PV.summary(pvTD.diag)
# 
# # compare
# c(pvTD.diag$modularity, lmod$LModularity)
# lmod.partition
# pvTD.diag$membership
# 
# ## use Gaussian transformation before analyzing the data
# pvTD.gauss = PV.complete(TestData, iseed=767, method="kendall", gauss=TRUE)
# PV.dendro(pvTD.gauss); PV.text(pvTD.gauss) 
# PV.graph.highlight(pvTD.gauss, col="lightblue")
# PV.summary(pvTD.gauss)

