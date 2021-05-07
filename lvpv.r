# Code for pvclust approximately unbiased p-value realisation wrt Newman's leading eigenvector algorithm for detection of integration and modularity in complex systems# Code for pvclust approximately unbiased p-value realisation wrt Newman's leading eigenvector algorithm for detection of integration and modularity in complex systems

# an implementation of a multiscale bootstrap procedure designed to estimate the p-values in the modular structure analysis with the leading eigenvector algorithm

# by Oksana Vertsimakha and Igor Dzeverin, 2021-04-27

# References
# Csardi, G., Nepusz, T. (2006). The igraph software package for complex network research. InterJournal, Complex Systems 1695. http://igraph.org
# Melo D., Garcia G., Hubbe A., Assis A. P., Marroig G. (2016). EvolQG вЂ“ An R package for evolutionary quantitative genetics [version 3; peer review: 2 approved, 1 approved with reservations]. F1000Research, 4(925), doi: 10.12688/f1000research.7082.3
# Newman, M. E. J. (2006). Modularity and community structure in networks. Proceedings of the National Academy of Sciences of the USA, 103, 8577-8582, doi: 10.1073/pnas.0601602103
# Shimodaira, H. (2004) Approximately unbiased tests of regions using multistep-multiscale bootstrap resampling. The Annals of Statistics, 32, 2616-2641. doi: 10.1214/009053604000000823


# COR.GRAPH creates graph from data's correlaion matrix, returns list of clusters vector, 
# merges matrix & modularity value;
# Note: function requires igraph package

# method and use are cor() function parameters;
# tuning corresponds to the fine-tuning procedure of the LE algorithm;
# the correlation mattrix is used in power of the beta parameter;
# option ="exp" uses exponential tranformation on the correlation matrix;
# if there are negative correlations after beta parameter is used, option = "exp" is enforced;
# if diag = TRUE the fraph is built with loops, default = TRUE;
# q.cut cuts off correlations below specified quantile (counted w/o diagonal), default = 0.5;


cor.graph = function(data, 
                     method = "pearson",
                     use = "everything",
                     exp = FALSE,
                     beta = 1,
                     q.cut = 0,
                     diag = TRUE,
                     tuning = FALSE,
                     warnings = FALSE){
  suppressMessages(require(igraph))
  C = cor(data, method = method, use = use)
  C = C^beta
  
  if((0<q.cut)&(q.cut<1)){diag(C) = 0; C[abs(C) < quantile(abs(C), q.cut)] = 0}
  
  if(exp) C = exp(-(max(C)-C)/sd(C))
  else if(any(C<0)){
    C = exp(-(max(C)-C)/sd(C))
    if(warnings) warning(paste("Some correlations are negative. Using exponential transformation."))
  }
  
  if(diag) diag(C) = 1
  
  g = graph.adjacency(C, weighted = TRUE, diag = diag, mode = 'undirected')
  if(tuning){
    res = NewmanSpectral(g, weights = TRUE, warnings = FALSE)
  } else{
    res = cluster_leading_eigen(g, options = list(maxiter = 1000000)) 
  }
  
  orig.cl = res$membership; orig.merges = res$merges
  if(length(colnames(data))==0) names(orig.cl) = sapply(1:ncol(data), function(j) paste("col", j))
  else{
    names(orig.cl) = colnames(data)
    for(j in 1:ncol(data)) if(names(orig.cl)[j]=="") names(orig.cl)[j] = paste("col", j)
  }
  
  orig.mod = res$modularity
  
  return(list(orig.cl, orig.merges, orig.mod))
}



SplitMerges = function(c.membership, c.merges){
  n <- length(c.membership)
  m <- nrow(c.merges)
  K <- max(c.membership)
  B <- list()
  
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

# inner function, implementation of the bootstrap procedure


PV.boot<-function(data, g.cor,
                  r, nboot = 1000,
                  method = "pearson", use = "everything", 
                  exp = TRUE,
                  warnings = FALSE, 
                  q.cut = 0, beta = 1, diag = TRUE,
                  store = FALSE, quiet = F, tuning = FALSE,
                  mod.alpha = .05){
  
  n <- nrow(data)
  size<- round(n*r, digits=0)
  if(size == 0)
    stop("Invalid scale parameter(r)")
  r <- size/n
  
  if(!is.na(g.cor[[3]])){
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
      x.cor = tryCatch(cor.graph(boot.data, method = method, use = use, 
                                 exp = exp, q.cut = q.cut, beta = beta, diag = diag, 
                                 tuning = tuning),
                       error=function(cond) {
                         return(as.list(rep(NA,3)))
                       }
      )
      if(sum(is.na(x.cor[[3]]))==0){
        
        pattern.i <- SplitMerges(x.cor[[1]], x.cor[[2]])$pattern   # split
        edges.cnt <- edges.cnt + table(factor(pattern.i,  levels=pattern))
        
      } else na.flag <- 1
      
      if(store) {st[[i]] <- x.cor}
      
      mod[[i]] = x.cor[[3]]
    }
    if((!quiet)&na.flag==0)
      cat("Done.\n")
    if((!quiet)&na.flag==1)
      cat("Error.\n")
    
    if(na.flag == 1)
      warning(paste("Unsuccessful networks divisions are omitted in computation: r = ", r), call.=FALSE)
    
    
    cmod = unlist(mod); cmod = cmod[!is.na(cmod)]
    
    if(length(cmod)>3){
      mod.mean = mean(cmod); mod.sd = sd(cmod)
      mod.q1 = quantile(cmod, mod.alpha/2)
      mod.q2 = quantile(cmod, 1 - mod.alpha/2)
    } else mod.mean<-mod.sd<-mod.q1 <-mod.q2 <-NA
    
    boot <- list(edges.cnt = edges.cnt,
                 
                 method = method, 
                 use = use,
                 q.cut = q.cut,
                 beta = beta,
                 diag = diag,
                 exp = exp,
                 tuning = tuning,
                 nboot = nboot, 
                 size=size, r=r, 
                 na.flag = na.flag, 
                 store=st,
                 
                 mod.mean = mod.mean, 
                 mod.sd = mod.sd,
                 mod.q1 = mod.q1, 
                 mod.q2 = mod.q2,
                 mod.alpha = mod.alpha)
  }
  else boot <- list(edges.cnt = NA,
                    
                    method = method, 
                    use = use,
                    q.cut = q.cut,
                    beta = beta,
                    diag = diag,
                    exp = exp,
                    tuning = tuning,
                    nboot = nboot, 
                    size=size, r=r, 
                    na.flag = na.flag, 
                    store=st,
                    
                    mod.mean = mod.mean, 
                    mod.sd = mod.sd,
                    mod.q1 = mod.q1, 
                    mod.q2 = mod.q2,
                    mod.alpha = mod.alpha)
  
  return(boot)
}

# inner function, AU and naive bootstrap (BP) p-value estimation

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


# Finds community structure of the inpout dataset; input options are usedin the cor.graph function;
# If bootstrap = TRUE, AU and BP p-values estimates are calculated;
# r = scaling parameters vector; 
# nboot = number ofiterations;
# iseed = set.seed option value
# mod.alpha is use in the naive bootstrap estimation of the modularity;

# returns list of the obtained p-values for each split & scale parameter, modularity bootstrap stats, input options etc



PV.complete<- function(data, 
                       bootstrap = TRUE,
                       nboot = 1000, 
                       r = seq(.5,1.4,by=.1),
                       
                       method="pearson", use="everything",
                       exp = TRUE, q.cut=0, 
                       beta = 1, diag = TRUE, tuning = FALSE,
                       mod.alpha = .05,
                       
                       warnings=FALSE,
                       
                       iseed=NULL, quiet=FALSE, store=FALSE){
  
  if(!is.null(iseed))
    set.seed(seed = iseed)
  
  g.full = cor.graph(data, method = method, use = use, q.cut = q.cut, beta = beta, tuning = tuning,
                     diag = diag, exp = exp, warnings = TRUE)
  
  if(bootstrap){
    mboot <- lapply(r, PV.boot, 
                    g.cor = g.full,
                    data = data, 
                    
                    method = method, 
                    use = use,
                    q.cut = q.cut,
                    beta = beta,
                    exp = exp,
                    diag = diag,
                    tuning = tuning,
                    nboot = nboot,
                    store = store,
                    mod.alpha = mod.alpha)
    
    ### merge results
    
    pattern <- SplitMerges(g.full[[1]], g.full[[2]])$pattern
    r <- unlist(lapply(mboot,"[[","r"))
    nboot <- unlist(lapply(mboot,"[[","nboot"))
    mod.mean = unlist(lapply(mboot,"[[","mod.mean"))
    mod.sd = unlist(lapply(mboot,"[[","mod.sd"))
    mod.q1 = unlist(lapply(mboot,"[[","mod.q1"))
    mod.q2 = unlist(lapply(mboot,"[[","mod.q2"))
    store_list <- lapply(mboot,"[[", "store"); names(store_list)=paste("r", r, sep = "")
    rl <- length(mboot); ne <- length(pattern)
    
    
    
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
      nboot1 = PV.boot(data, r=1, g.cor = g.full,
                       nboot = nboot[[1]],
                       method = method, use = use, 
                       exp = exp,
                       tuning = tuning,
                       q.cut = q.cut, beta = beta, diag = diag,
                       mod.alpha = mod.alpha,
                       
                       warnings = warnings,
                       store = F, quiet= T
                       )
      
      mod.stats[[1]] = nboot1$mod.mean
      
      mod.stats[[2]] = nboot1$mod.sd
      
      mod.stats[[3]] = nboot1$mod.q1
      
      mod.stats[[4]] = nboot1$mod.q2
    }
    
    mod.stats = data.frame(t(unlist(mod.stats))); 
    colnames(mod.stats) = c("mean", "sd", "q1", "q2")

    pv.result = list(
      membership = g.full[[1]],
      modularity = g.full[[3]],
      merges = g.full[[2]],
      modularity.bp = modularity.bp,
      mod.stats = mod.stats,
      
      count = edges.cnt,
      edges.bp = edges.bp,
      nboot = nboot, r=r,
      edges = edges.pv,
      method = list(cor.method = method, use = use, exp = exp, 
                    q.cut = q.cut, beta = beta, diag = diag, 
                    tuning = tuning,
                    bootstrap = TRUE, store = store,
                    mod.alpha = mod.alpha, seed = iseed),
      
      store=store_list
    )
  }
  else{
    pv.result = list(
      membership = g.full[[1]],
      modularity = g.full[[3]],
      merges = g.full[[2]],
      
      method = list(cor.method = method, use = use, exp = exp, 
                    q.cut = q.cut, beta = beta, diag = diag, 
                    tuning = tuning, bootstrap = F)
    )
  }
  return(pv.result)
}

PV.complete = compiler::cmpfun(PV.complete)


# Produces summary report on the PV.complete function from the output;
# which gives detaled information on the selected modules only (referenced by the numbers)
# modularity gives information on the modularity bootstrap estimates

PV.summary<-function(x, which = NULL, digits = 3, modularity = TRUE){
  d = TRUE
  n = 2*nrow(x$merges); K = length(unique(x$membership->cl))
  
  if(is.null(which)) which = 1:n
  mtch = intersect(which, 1:K)
  
  if(n == 1) {print("One cluster detected"); d = F}
  
  if(length(mtch)==0) {print(paste("Wrong cluster indexes,", K, "clusters found", sep = " ")); d = F}
  
  if(d){
    names = names(cl)
    if(is.null(names)) names = 1:length(cl) 
    
    cat("Clusters:\n\n")
    for(j in mtch){
      
      cat(paste("C", j, ": ", sep = ""))
      
      cat(paste(names[cl==j], sep = ", "), labels = NULL)
      
      cat("\n")
    }
    if(x$method[[8]]){
      cat("Estimates on edges:\n\n")
      
      print(round(x$edges[mtch,], digits=digits))
      
      cat("\n")
    }

    
    cat(paste("Correlation method: ", x$method[[1]],', use = ', x$method[[2]], "\n", sep = ""))
    
    if(x$method[[3]]=="exp") cat('Exponential transformation used',"\n")
    if(x$method[[5]]!=1) cat('Correlation matrix in ', x$method[[5]], ' power used', "\n")
    if((is.numeric(x$method[[4]]))&(x$method[[4]]>0)) cat('Correlations cut at', 
                                                          paste(round(100*x$method[[4]],2),"%",sep=""), 
                                                          'quantile', "\n")
    if(x$method[[6]]) cat('diag = TRUE, graph built with loops', "\n")
    if(x$method[[7]]) cat('Fine-tuning procedure used in the community detection', "\n")
  }
  
  if(x$method[[8]]&modularity) cat(paste('Modularity: ', round(x$modularity, 4), "\n", 
                           'Bootstrap estimate: ', round(x$mod.stats["mean"],4), "\n",
                           'Bootstrap sd: ', round(x$mod.stats["sd"],4), "\n", 
                           'Bootstrap ', paste(round(100*x$method$mod.alpha/2,2),"%",sep=""), ' quantile: ',
                           round(x$mod.stats["q1"], 4), "\n", 
                           'Bootstrap ', paste(round(100*(1-x$method$mod.alpha/2),2),"%",sep=""), ' quantile: ',
                           round(x$mod.stats["q2"], 4), sep = ""))
}

lv.hc2axes <- function(x, K, m){
  x = as.hclust(x)
  A <- x$merge 
  
  x.axis = 1:K
  y.axis <- c(rep(0.5, 2*m)[as.vector(t(A)<0)], x$height)
  
  x.tmp  <- rep(0,2)
  zz <- match(1:K, x$order)
  
  for(i in 1:m){
    
    ai <- A[i,1]
    x.tmp = ifelse(ai < 0,zz[-ai],x.axis[ai+K])
    ai <- A[i,2]
    x.tmp[2] <-ifelse(ai < 0,zz[-ai], x.axis[ai+K])
    
    x.axis[i+K] <- mean(x.tmp)
  }
  
  rm(x.tmp, zz)
  return(data.frame(x.axis=x.axis, y.axis=y.axis))
}

make_dendro = function(xC, xM, K, m){
  
  a = list()
  
  if(m==0|m==1) a = NULL
  else if(m > 1)
  {
    new_xM = apply(xM, 2, function(z){(z>K)*(z-K)+(z<=K)*(-z)})
    a$merge = new_xM
    a$height = 1:m
    a$order = -t(new_xM)[t(new_xM)<0]
    a$labels <- sapply(1:K, function(z) paste("C", z, sep = ""))
    class(a) <- "hclust"      
  }
  return(a)
}

# Creates a dendrogram based on the PV.complete function output. 
# members gives labels of the modules' (can be referenced by the numbers in the which option) elements;
# print.pv = FALSE hides the estimates; 

PV.dendro <-function(x, members = FALSE, which = NULL,
                     print.pv=TRUE, float=0.01,
                     col.pv=c(au=2, bp=3), cex.pv=0.8, font.pv=NULL,
                     col=NULL, cex=.8, font=NULL, lty=NULL, lwd=NULL,
                     main=NULL, sub=NULL, xlab=NULL, ylab = NULL, ...){
  
  xC = x$membership; xM = x$merges
  K = length(unique(xC)); m = nrow(xM)
  
  a = make_dendro(xC, xM, K, m)
  
  if(is.null(main))
    main ="Cluster dendrogram with p-values (%)"
  
  if(is.null(xlab))
    xlab = ""
  
  if(is.null(ylab))
    ylab = ""
  
  if(is.null(cex))
    cex = .7
  
  if(members){
    names = names(xC)
    if(is.null(which)) which = 1:(2*m)
    mtch = intersect(which, 1:K)
    
    a$labels[mtch] = sapply(mtch, 
                            function(j){paste(paste("C", j, ": ", sep = ""), 
                                              paste(names[xC==j], collapse = ", "), sep = "")})
  }
  if(m<2) cat("Less then three clusters, no dendrogram available", "\n")
  else if(m>=2) {
    plot(as.dendrogram(a), main=main, sub=sub, 
         xlab=xlab, ylab=ylab, col=col,
         cex=cex, font=font, lty=lty, lwd=lwd, ...)
  }
  
  
  if(isTRUE(print.pv)){
    print.pv <- c("au", "bp"); col.text <- col.pv[print.pv]
    PV.text(x, col=col.text, cex=cex.pv, font=font.pv, float=float, print.num=print.num)
  }
}

# Adds text of the estimated BP and AU p-values next to the splits to the dendrogram by PV.dendro;
# PV.text is part ofthe PV.ddemndro function called on by the print.pv option;
# to hide au or bp values set respective color to transparent


PV.text <- function(x, col=c(au=2, bp=3), 
                    print.num=TRUE, float=0.01, cex=NULL, font=NULL, ...)
{
  xC = x$membership; xM = x$merges
  K = length(unique(xC)); m = nrow(xM)
  a = as.dendrogram(make_dendro(xC, xM, K, m))
  a = as.hclust(a)
  xC = a$membership; xM = a$merges
  
  if(length(a)>0){
    if(length(col) == 1) col = rep(col, 2)
    names(col) <- c("au", "bp")
    
    axes <- lv.hc2axes(a, K, m)
    usr  <- par()$usr; wid <- usr[4] - usr[3]
  
    
    range <- seq_len(min(3, length(col)))
    pos <- c(2, 4, 1)
    y_off <- float * wid
    y_offset<- c(rep(-y_off*3, K), rep(y_off*2, m))
    
    ind = c(a$order, K+(1:m))
    if(x$method[[8]]){
      num_str <- lapply(
        x$edges[seq_len(which(names(x$edges) == "bp"))],
        function(p) round(p * 100))
      for(i in names(num_str)) {
        num_str[[i]][length(num_str[[i]])] <- i
      }
      for(i in range) {
        name <- names(col)[i]
        text(x=axes[,1], y=axes[,2] + y_offset, num_str[[name]][ind],
             col=col[name], pos=pos[i], offset=.2, cex=cex, font=font, ...)
      }
    }
    
  }
}


# Highlights selected modules in the PV.dendro output dendrogram;
# pv = "au" or "bp"; type and alpha specify conditions (see pvclust options);
# fill.col and fill.alpha define filled area color and transparence;


PV.rectangle <- function(x, alpha = 0.95, pv = "au", type = "geq", max.only = T,
                         fill.col = "cyan", fill.alpha = .4, border = NULL,...){
  
  m = nrow(x$merges)
  K = max(x$membership)
  
  order = t(x$merges)[t(x$merges)<=K]
  member = SplitMerges(1:K, x$merges)$member
  
  fill.col = col2rgb(fill.col)/max(col2rgb(fill.col))
  col = rgb(fill.col[1], fill.col[2], fill.col[3], alpha = fill.alpha)
  
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
        
        mx <- xwd / length(member) / 5.2
        my <- ywd / 300
        
        rect(xl - mx, yb + my, xr + mx, yt + .8*my, 
             border=border, col=col, ...)
        j <- j + 1
      }
      ht <- c(ht, ma)
      
    }
  }
}



# Creates graphbased on the PV.complete function output;
# members gives labels of the modules' (can be referenced by the numbers in the which option) elements;
# print.pv = FALSE hides the estimates;
# highlight option colors selected modules (see PV.rectangle);


PV.graph = function(x, members = TRUE, 
                    print.pv = TRUE,
                    m.labels = TRUE,
                    which = NULL,
                    highlight = FALSE,
                    alpha=0.95, pv ="au", type="geq", 
                    
                    fill.col = "cyan",
                    fill.alpha = 0.4,
                    
                    vertex.shape = "rectangle",
                    vertex.color = "white",
                    vertex.frame.color = "black",
                    vertex.size = 16,
                    vertex.label.color = "black",
                    vertex.label.cex = .9,
                    vertex.label.color2 = NULL,
                    vertex.label.cex2 = .9, ...
                    
){
  suppressMessages(require(igraph))
  xC = x$membership; xM = x$merges
  K = max(xC); m = nrow(xM)
  if(m==0) stop("No clusters, no graph")
  
  ##add au bp test here
  fill.col = col2rgb(fill.col)/max(col2rgb(fill.col))
  col = rgb(fill.col[1], fill.col[2], fill.col[3], alpha = fill.alpha)
  
  G<-graph.tree(n = 3, children = 2, mode = "out")
  V(G)$name = c("", paste(xM[m,]))
  if(m>=2){
    for(j in (m-1):1){
      G <- G + vertices(paste(xM[j,]))
      G <- G + edge(paste(K+j), paste(xM[j,1]))
      G <- G + edge(paste(K+j), paste(xM[j,2]))
      
    }
  }
  #print(V(G))
  nn = 0
  index = c(1, rep(2,2*m), rep(3, nn))
  v.labels = V(G)$name
  
  
  if(members){
    
    if(is.null(which)) {which = 1:K}
    include = xC %in% which; xCin = xC[include]
    
    names = names(xCin)
    
    nn = sum(include)
  
      
    for(k in 1:nn){
      if(names[k] %in% V(G)) names[k] = paste("v", names[k])
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
  if(x$method[[8]]&print.pv)
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
  v.label.cex = c(vertex.label.cex, vertex.label.cex, vertex.label.cex2)[index]
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
  vc = vertex.attributes(G)$color
  
  e = length(E(G))
  edge_attr(G) = list(
    size = rep(.8, e),
    arrow.size = rep(.2, e)
  )
  graph_attr(G, "layout") = co
  
  if(x$method[[8]]&highlight){
    ht <- c()
    j <- 1
    if(is.na(pm <- pmatch(type, c("geq", "leq", "gt", "lt"))))
      stop("Invalid type argument: see help(pvrect)")
    
    if(pv == "au") w = 1 else if(pv == "bp") w = 2
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
    
    vc[ht+1] = col
    
  }
  plot(G,
       vertex.color = vc,
       layout = graph_attr(G, "layout"))
  #return(G)
}

#######

library(MASS)
library(Matrix)
dims = c(3,4,2,5);
sig = .1
d = length(dims); D = sum(dims); cl = rep(1:d, dims)
X1 = mvrnorm(100, rep(0, d), diag(rep(1, d)))
X2 = mvrnorm(100, rep(0, D), diag(rep(sig, D)))             

corvec = rnorm(D, 0.75, .1)
X = X1[,1]*corvec[1]
for(j in 2:D) X = cbind(X, X1[,cl[j]]*corvec[j])

Xc = cor(X+X2)

XX = mvrnorm(250, rep(0,D), Xc)
xpv1 = PV.complete(XX, nboot = 500)
corrplot::corrplot(cor(XX), method = "color", tl.col = xpv1$membership)



# an implementation of a multiscale bootstrap procedure designed to estimate the p-values in the modular structure analysis with the leading eigenvector algorithm

# by Oksana Vertsimakha and Igor Dzeverin, 2020-11-03

# References
# Csardi, G., Nepusz, T. (2006). The igraph software package for complex network research. InterJournal, Complex Systems 1695. http://igraph.org
# Melo D., Garcia G., Hubbe A., Assis A. P., Marroig G. (2016). EvolQG – An R package for evolutionary quantitative genetics [version 3; peer review: 2 approved, 1 approved with reservations]. F1000Research, 4(925), doi: 10.12688/f1000research.7082.3
# Newman, M. E. J. (2006). Modularity and community structure in networks. Proceedings of the National Academy of Sciences of the USA, 103, 8577-8582, doi: 10.1073/pnas.0601602103
# Shimodaira, H. (2004) Approximately unbiased tests of regions using multistep-multiscale bootstrap resampling. The Annals of Statistics, 32, 2616-2641. doi: 10.1214/009053604000000823


# creates graph from data's cor matrix, returns clusters, merges & modularity;
# (function requires igraph package and rARPACK package for the fine-tuning implementation)


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
                     beta = 1,
                     q.cut = .5,
                     diag = FALSE,
                     square = TRUE,
                     tuning = FALSE,
                     warnings = FALSE){
  suppressMessages(require(igraph))
  C = cor(data, method = method, use = use)
  C = C^beta
  
  if(any(C<0)&square){
    C = C^2
    if(warnings) warning("Some correlations are negative. Using squared correlations.")
  }
  
  if((0<q.cut)&(q.cut<1)){diag(C) = 0; C[abs(C) < quantile(abs(C), q.cut)] = 0}
  
  if(gauss|!square) C = exp(-(1-abs(C))^2/(2*sd(C)^2))
  
  if(diag) diag(C)=1
  
  g = graph.adjacency(C, weighted = TRUE, diag = diag, mode = 'undirected')
  if(tuning){
    res = NewmanSpectral(g, weights = TRUE)
  } else{
    res = cluster_leading_eigen(g, options = list(maxiter = 1000000)) 
  }
  
  orig.cl = res$membership; orig.merges = res$merges
  if(!is.null(colnames(data))){names(orig.cl)<-colnames(data)} 
  else{names(orig.cl) = paste("col", 1:ncol(data), sep = "")}
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
                  warnings = F, gauss = F, q.cut = 0.5, beta = 1, diag = F, square = square,
                  store = F, quiet = F, tuning = F,
                  mod.alpha = .05){
  
  n <- nrow(data)
  size<- round(n*r, digits=0)
  if(size == 0)
    stop("invalid scale parameter(r)")
  r <- size/n
  
  g.cor = cor.graph(data, method = method, use = use, 
                    gauss = gauss, q.cut = q.cut, beta = beta, 
                    diag = diag, square = square,
                    tuning = tuning, warnings = warnings)
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
    
    x.cor = cor.graph(boot.data, method = method, use = use, gauss = gauss, 
                      square = square, q.cut = q.cut, beta = beta, diag = diag, tuning = tuning)
    
    mod[[i]] = x.cor[[3]]
    
    if(sum(is.na(x.cor[[2]]))==0){
      
      pattern.i <- SplitMerges(x.cor[[1]], x.cor[[2]])$pattern   # split
      edges.cnt <- edges.cnt + table(factor(pattern.i,  levels=pattern))
      
    } else na.flag <- 1
    
    if(store)
      st[[i]] <- x.cor
  }
  if(!quiet)
    
    cat("Done.\n")
  
  if(na.flag == 1)
    warning(paste("unsuccessful networks divisions are omitted in computation: r = ", r), call.=FALSE)
  
  
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
               beta = beta,
               diag = diag,
               gauss = gauss,
               square = square,
               tuning = tuning,
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
                       warnings=FALSE, gauss=FALSE, square=TRUE, q.cut=.5, beta = 1, diag=FALSE,
                       iseed=NULL, quiet=FALSE, tuning = FALSE, store=FALSE,
                       mod.alpha = .05){
  
  if(!is.null(iseed))
    set.seed(seed = iseed)
  
  
  mboot <- lapply(r, PV.boot, data = data, 
                  method = method, use = use,
                  q.cut = q.cut,
                  beta = beta,
                  gauss = gauss,
                  square = square,
                  diag = diag,
                  tuning = tuning,
                  nboot = nboot,
                  mod.alpha = mod.alpha)
  
  ### merge results
  g.full = cor.graph(data, method = method, use = use, q.cut = q.cut, beta = beta,  tuning = tuning,
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
                   square = square, tuning = tuning,
                   q.cut = q.cut, beta = beta, diag = diag,
                   store = F, quiet= T,
                   mod.alpha = mod.alpha)
    
    mod.stats[[1]] = nboot1$mod.mean
    
    mod.stats[[2]] = nboot1$mod.sd
    
    mod.stats[[3]] = nboot1$mod.q1
    
    mod.stats[[4]] = nboot1$mod.q2
  }
  
  mod.stats = data.frame(t(unlist(mod.stats))); 
  colnames(mod.stats) = c("mean", "sd", "q1", "q2")
  
  
  pv.result = list(
    membership = g.full[[1]],
    modularity = g.full[[3]],
    modularity.bp = modularity.bp,
    mod.stats = mod.stats,
    merges = g.full[[2]],
  
    count = edges.cnt,
    edges.bp = edges.bp,
    nboot = nboot, r=r,
    edges = edges.pv,
    method = list(cor.method = method, use = use, gauss = gauss, 
                  q.cut = q.cut, beta = beta, diag = diag, square = square, 
                  mod.alpha = mod.alpha, tuning = tuning, seed = iseed),
    
    store=store
  )
  return(pv.result)
}
PV.complete = compiler::cmpfun(PV.complete)


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
    names = names(cl)
    if(is.null(names)) names = 1:length(cl) 
    
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
    if(!x$method[[7]]) cat('Negative correlations used',"\n")
    if(x$method[[3]]) cat('Gauss transformation used',"\n")
    if(x$method[[5]]!=1) cat('Correlation matrix in power', x$method[[5]], 'used', "\n")
    if((is.numeric(x$method[[4]]))&(x$method[[4]]>0)) cat('Correlations cut at', 
                                                          paste(round(100*x$method[[4]],2),"%",sep=""), 
                                                          'quantile', "\n")
    if(x$method[[6]]) cat('diag = TRUE, graph built with loops', "\n")
    if(x$method[[9]]) cat('Fine-tuning procedure used in the community detection', "\n")
  }
  if(modularity) cat(paste('Modularity: ', round(x$modularity, 4), "\n", 
                           'Bootstrap estimate: ', round(x$mod.stats["mean"],4), "\n",
                           'Bootstrap sd: ', round(x$mod.stats["sd"],4), "\n", 
                           'Bootstrap ', paste(round(100*x$method$mod.alpha/2,2),"%",sep=""), ' quantile: ',
                           round(x$mod.stats["q1"], 4), "\n", 
                           'Bootstrap ', paste(round(100*(1-x$method$mod.alpha/2),2),"%",sep=""), ' quantile: ',
                           round(x$mod.stats["q2"], 4), sep = ""))
}





# (inner function)
# transforms merge matrix for the dendrogram parameters purposes
# returns data.frame with x and y coordinates of the nodes

lv.hc2axes <- function(x, K, m){
  
  A <- x$merge 
  
  x.axis = 1:K
  y.axis <- c(rep(1:m, each = 2)[as.vector(t(A)<0)], x$height)
  x.tmp  <- rep(0,2)
  zz <- match(1:K, x$order)
  
  for(i in 1:m){
    
    ai <- A[i,1]
    x.tmp = ifelse(ai < 0,zz[-ai],x.axis[ai+K])
    ai <- A[i,2]
    x.tmp[2] <-ifelse(ai < 0,zz[-ai], x.axis[ai+K])
    
    x.axis[i+K] <- mean(x.tmp)
  }
  
  rm(x.tmp, zz)
  return(data.frame(x.axis=x.axis, y.axis=y.axis))
}

# (inner function)
# creates dendrogram from the membership vector and merges matrix (of igraph type)
# returns hclust object or empty list if only one cluster was detected

make_dendro = function(xC, xM, K, m){
  
  a = list()
  
  if(m==0|m==1) cat("Less then three clusters, no dendrogram available", "\n")
  else if(m > 1)
    {
    new_xM = apply(xM, 2, function(z){(z>K)*(z-K)+(z<=K)*(-z)})
    a$merge = new_xM
    a$height = 1:m
    a$order = -t(new_xM)[t(new_xM)<0]
    a$labels <- sapply(1:K, function(z) paste("C", z, sep = ""))
    class(a) <- "hclust"      
  }
  return(a)
}

# creates dendrogram according to leading eigenvector communities
# x = PV.complete function output, additional parameters (...) for the plot function; 
# print.pv adds AU and BP p-value estimates (default; otherwise, can be done with PV>text function)
# members option adds module members labels; which with allows to specify particular modules by their numbers


PV.dendro <-function(x, members = FALSE, which = NULL,
                     print.pv=TRUE, float=0.01,
                     col.pv=c(au=2, bp=3), cex.pv=0.8, font.pv=NULL,
                     col=NULL, cex=NULL, font=NULL, lty=NULL, lwd=NULL,
                     main=NULL, sub=NULL, xlab=NULL, ylab = NULL, ...){
  
  xC = x$membership; xM = x$merges
  K = length(unique(xC)); m = nrow(xM)
  a = make_dendro(xC, xM, K, m)
  
  if(is.null(main))
    main ="Cluster dendrogram with p-values (%)"
  
  if(is.null(xlab))
    xlab = ""
  
  if(is.null(ylab))
    ylab = ""
  
  if(is.null(cex))
    cex = .7

  if(members){
    names = names(xC)
    if(is.null(which)) which = 1:(2*m)
    mtch = intersect(which, 1:K)
    
    a$labels[mtch] = sapply(mtch, 
                            function(j){paste(paste("C", j, ": ", sep = ""), 
                                              paste(names[xC==j], collapse = ", "), sep = "")})
  }
  if(length(xM)>2) {
    plot(a, main=main, sub=sub, 
         xlab=xlab, ylab=ylab, col=col,
         cex=cex, font=font, lty=lty, lwd=lwd, ...)
  }
  
  
  if(isTRUE(print.pv)){
    print.pv <- c("au", "bp"); col.text <- col.pv[print.pv]
    PV.text(x, col=col.text, cex=cex.pv, font=font.pv, float=float, print.num=print.num)
  }
}

# adds estimated BP and AU p-values next to the splits to the dendrogram by PV.dendro
# x = PV.complete function output, text parameters
# to hide au or bp values setrespective colorto transparent


PV.text <- function(x, col=c(au=2, bp=3), 
                    print.num=TRUE, float=0.01, cex=NULL, font=NULL, ...)
  {
  xC = x$membership; xM = x$merges
  K = length(unique(xC)); m = nrow(xM)
  a = make_dendro(xC, xM, K, m)
  
  if(length(a)>0){
    if(length(col) == 1) col = rep(col, 2)
    names(col) <- c("au", "bp")

    axes <- lv.hc2axes(a, K, m)
    usr  <- par()$usr; wid <- usr[4] - usr[3]
    num_str <- lapply(
      x$edges[seq_len(which(names(x$edges) == "bp"))],
      function(p) round(p * 100))
  
  
  # change the last element to the name of p-value
    for(i in names(num_str)) {
      num_str[[i]][length(num_str[[i]])] <- i
    }
  
    range <- seq_len(min(3, length(col)))
    pos <- c(2, 4, 1)
    y_off <- float * wid
    y_offset<- c(rep(-y_off*3, K), rep(y_off*2, m))

    ind = c(a$order, K+(1:m))
    
    for(i in range) {
      name <- names(col)[i]
      text(x=axes[,1], y=axes[,2] + y_offset, num_str[[name]][ind],
           col=col[name], pos=pos[i], offset=.2, cex=cex, font=font, ...)
    }
    
  }
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
        
        mx <- xwd / length(member) / 5.2
        my <- ywd / 300
        
        rect(xl - mx, yb + my, xr + mx, yt + .8*my, 
             border=border, shade=NULL, ...)
        j <- j + 1
      }
      ht <- c(ht, ma)
      
    }
  }
}


### creates via igraph a binary graph corresponding to detected module structure
# input= PV.complete ouput; members = create nodes for membersof each final cluser, 
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
    
    if(is.null(names(xCin))){names = paste("col", 1:length(xC))} else {
      names = names(xCin); 
      names[names(xCin)==""]=paste("col", 1:length(xC), sep = "")[names(xCin)==""]} 
    
    nn = sum(include)
    
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
  v.label.cex = c(vertex.label.cex, vertex.label.cex, vertex.label.cex2)[index]
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
# input: x = PV.complete output, ... = PV.graph visual options
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



# Fine-tuning procedure for a single split (see Newman (2006) for details)
ft_method<-function(s,B,m){
  eps = length(s); notMv = 1:eps
  dQ_full = 0; 
  s1 = s;
  count = 0; upd = T
  
  while(count < eps){
    if(upd){
      Bs = B%*%s; upd = F
    }
    dQ1 = -1; dQ_id = 0
    
    dQ_vec = (-s*Bs)[notMv]/m; dQ = max(dQ_vec)
    if(dQ > dQ1){
      dQ1 = dQ; i = which.max(dQ_vec); dQ_id = notMv[i]
    }
    s[dQ_id] = -s[dQ_id]; notMv = notMv[-i]
    count = count + 1
    if(dQ1>=0){s1 <- s; upd <- T}
  }
  return(s1)
}


# Split of clusters with fine-tuning procedure correction
SubSplit<-function(B,m,c_g){
  res = list(FALSE)
  if(length(c_g)>1){
    B_g<-B[c_g,c_g]-diag(apply(B[c_g,c_g], 1, sum))
    
    if(nrow(B_g)==2){
      s_g = c(1,-1)
    }
    else{
      u_g<-eigs_sym(B_g, 1, which = "LA", options = list(maxitr = 1000000))$vectors[,1]
      s_g<-(u_g>=0)*2-1
    }
    
    s_g = ft_method(s_g,B_g,m)
    dQ = t(s_g)%*%B_g%*%s_g/4/m
    
    if((dQ>=0)&(abs(sum(s_g))<length(s_g))){
      res=list(T,c_g[s_g>0],c_g[s_g<0])
    }

  }
  return(res)
}


# Leading eigenvector algorithm with fine-tuning procedure implemented
# G = input graph; 
# default weights = TRUE uses graph weight attribute (equal 1 for unweighted graph)
# can be specified vector of length E(G) or else weights = NULL to ignore weights attribute

NewmanSpectral<-function(G, weights = TRUE){
  suppressMessages(require(igraph))
  suppressMessages(require(rARPACK))
  
  N = length(V(G)); comm_num = 1
  if(is.null(E(G)$weight)){
    E(G)$weight = 1
  }else if(length(weights)==N){
    E(G)$weight = weights
  }
  
  A<-as.matrix(get.adjacency(G, attr = "weight"))
  m = sum(A)/2
  B<-A-apply(A, 1, sum)%*%t(apply(A, 1, sum))/2/m
  u<-eigs_sym(B, 1, which = "LA", options = list(maxitr = 1000000))$vectors[,1]
  s<-(u>=0)*2-1; s0 = ft_method(s,B,m)
  comm = c()
  
  if(abs(sum(s0))==length(s0)){
    print("No splits made")
  }
  else{
    go = T; splt = c(); nn = 2
    comm_temp = list((1:N)[s0>0], (1:N)[s0<0])
    
    while(go){
      
      go = F; k = 0
      comm_temp1 = list()
      
      for(i in 1:nn){
        res = SubSplit(B,m,comm_temp[[i]])
        splt = c(splt, res[[1]])
        
        if(res[[1]]){
          comm_temp1 = append(comm_temp1, list(res[[2]]))
          comm_temp1 = append(comm_temp1, list(res[[3]]))
          go = T; k = k+2
        }
        else{
          comm = append(comm, list(comm_temp[[i]]))
        }
      }
      comm_temp = comm_temp1; nn = k
    }
  }
  
  K = length(comm); KK = K*2-2
  
  if(K==0){
    members = list(1:N)
    membership = rep(1,N)
    merge = matrix(ncol = 2)
    mod = 0
  }
  else{
    v1 = 1:KK; v2 = c(v1[!splt], v1[splt][(K-2):1])
    merge = matrix(match(v1, v2)[KK:1], ncol = 2, byrow = T)
    members = comm
    membership = 1:N 
    for(c in 1:K) membership[unlist(comm[[c]])]=c
    AA = A-apply(A, 1, sum)%*%t(apply(A, 1, sum))/2/m
    mod = sum(sapply(1:K, function(z){
      sum(AA[membership==z, membership==z])/2/m}))
  }
  
  result = list(
    membership = membership,
    merges = merge,
    members = members,
    modularity = mod
  )
  return(result)
}
NewmanSpectral = compiler::cmpfun(NewmanSpectral)
