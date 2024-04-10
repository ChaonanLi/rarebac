## 功能冗余计算函数
# 李超男 20210110
# 多线程优化版本

### 距离矩阵计算
parallelDistCal <- function(species2Fun, distance = 'bray', threads = 100){
  dis <- parallelDist::parallelDist(as.matrix(species2Fun), 
                                    method = distance, 
                                    diag = FALSE, upper = FALSE, 
                                    threads = threads)
  return (dis)
}

### 定义功能冗余计算执行函数
funK <- function(aListIn) {
  sample <- as.character(aListIn$sample)
  v <- as.numeric(aListIn$values)
  tol <- as.numeric(aListIn$tol)
  D <- as.matrix(aListIn$distance)
  p <- v/sum(v)
  K <- apply(D, 1, function(x) sum(x*p))
  K[p<tol] <- NA
  res <- c(sample, as.vector(K))
  return(res)
}
funKbar <- function(aListIn){
  sample <- as.character(aListIn$sample)
  v <- as.numeric(aListIn$values)
  tol <- as.numeric(aListIn$tol)
  D <- as.matrix(aListIn$distance)
  p <- v/sum(v)
  Kbar <- sapply(1:nrow(D), function(i) sum(D[i,]*v/sum(v[-i])))
  Kbar[p<tol] <- NA
  res <- c(sample, as.vector(Kbar))
  return(res)
}
funQ <- function(aListIn){
  sample <- as.character(aListIn$sample)
  v <- as.numeric(aListIn$values)
  D <- as.matrix(aListIn$distance)
  p <- v/sum(v)
  Q <- t(p)%*%D%*%p
  res <- c(sample, as.vector(Q))
  return(res)
}
funSim <- function(aListIn){
  sample <- as.character(aListIn$sample)
  v <- as.numeric(aListIn$values)
  p <- v/sum(v)
  S <- 1-sum(p^2)
  res <- c(sample, S)
  return(res)
} 
funN <- function(aListIn){
  sample <- as.character(aListIn$sample)
  v <- as.numeric(aListIn$values)
  tol <- as.numeric(aListIn$tol)
  p <- v/sum(v)
  N <- length(p[p>tol])
  res <- c(sample, N)
  return(res)
}

### 定义多线程执行框架
parallellyRun <- function(listIn, threads, func){
  #cat("          -- Register parallel cores...", fill = T)
  cl <- makeCluster(threads, type = "PSOCK")
  #cat("          -- Computation start...", fill = T)
  system.time({
    res <- parLapply(cl, listIn, func)
  })
  final.res <- do.call(rbind, res)
  final.res <- as.data.frame(final.res)
  stopCluster(cl)
  return(final.res) 
}

## 功能冗余函数（多线程优化）
FunctionalRedParallelDist <- function(comm, dis, threads = 100, tol = 1e-8, abundance = TRUE){
  # comm: species in row and sample in col 
  # dis: species x species 
  
  # 导入包
  library(doParallel)
  library(parallelDist)
  
  # 按照距离矩阵对群落数据进行排序
  cat(paste("** [", Sys.time(), '] ** start to match species orders', sep = ''), fill = T)
  speciesKeep <- attributes(dis)$Labels
  idx <- sapply(speciesKeep, function(x){
    which(rownames(comm) == x)
  })
  comm <- as.data.frame(comm[idx,])
  if(length(unique(rownames(comm) == speciesKeep)) > 1 | !unique(rownames(comm) == speciesKeep)){
    stop("物种顺序不一致！请检查数据集！")
  }
  
  # 将distance转化为矩阵
  D <- as.dist(dis)
  
  # 是否考虑丰度
  if(!abundance) {
    comm[comm>0] <- 1
  }
  
  # 转置群落表
  commt <- as.data.frame(t(comm))
  
  # 准备多线程计算的数据集
  cat(paste("** [", Sys.time(), '] ** start to format datasets', sep = ''), fill = T)
  listIn <- list()
  for (i in seq(nrow(commt))){
    listIn[[i]] <- list(
      sample = rownames(commt)[i],
      values = as.numeric(commt[i,]),
      tol = tol,
      distance = D
    )
  }
  
  # 开始多线程计算
  ## 计算K
  cat(paste("** [", Sys.time(), '] ** start to calculate [K]', ' using ', threads, ' cores...', sep = ''), fill = T)
  K <- parallellyRun(listIn = listIn, threads = threads, func = funK)
  rownames(K) <- K[,1]; K <- K[,-1]
  colnames(K) <- colnames(commt)
  
  ## 计算Kbar
  cat(paste("** [", Sys.time(), '] ** start to calculate [Kbar]', ' using ', threads, ' cores...', sep = ''), fill = T)
  Kbar <- parallellyRun(listIn = listIn, threads = threads, func = funKbar)
  rownames(Kbar) <- Kbar[,1]; Kbar <- Kbar[,-1]
  colnames(Kbar) <- colnames(commt)
  
  ## 计算Q
  cat(paste("** [", Sys.time(), '] ** start to calculate [Q]', ' using ', threads, ' cores...', sep = ''), fill = T)
  Q <- parallellyRun(listIn = listIn, threads = threads, func = funQ)
  colnames(Q) <- c("SampleID", "Q")
  Q$Q <- as.numeric(as.vector(Q$Q))
  
  ## 计算Sim
  cat(paste("** [", Sys.time(), '] ** start to calculate [D]', ' using ', threads, ' cores...', sep = ''), fill = T)
  Sim <- parallellyRun(listIn = listIn, threads = threads, func = funSim)
  colnames(Sim) <- c("SampleID", "Sim")
  Sim$Sim <- as.numeric(as.vector(Sim$Sim))
  
  ## 计算U
  cat(paste("** [", Sys.time(), '] ** start to calculate [U]', ' using ', threads, ' cores...', sep = ''), fill = T)
  U <- data.frame(SampleID = Q$SampleID, U = Q$Q/Sim$Sim)
  
  ## 计算N
  cat(paste("** [", Sys.time(), '] ** start to calculate [N]', ' using ', threads, ' cores...', sep = ''), fill = T)
  N <- parallellyRun(listIn = listIn, threads = threads, func = funN)
  colnames(N) <- c("SampleID", "N")
  N$N <- as.numeric(as.vector(N$N))
  
  # 合并数据并返回结果
  cat(paste("** [", Sys.time(), '] ** start to combines results...', sep = ''), fill = T)
  red <- data.frame(N=N$N, D=Sim$Sim, Q=Q$Q, U=U$U)
  rownames(red) <- rownames(commt)
  res <- list()
  res$Kbar <- Kbar
  res$V <- K
  res$red <- red
  res 
}

## 功能冗余计算函数(原始函数)
uniqueness <- function(comm, dis, tol = 1e-8, abundance = TRUE) {
  if(!is.null(colnames(comm)) & !is.null(attributes(dis)$Labels)) {
    if(any(!colnames(comm)%in%attributes(dis)$Labels)) stop("One or several species in comm are not in dis; check species names in comm and in dis")
    else dis <- as.dist(as.matrix(dis)[colnames(comm), colnames(comm)])
  }
  else if(ncol(comm)!=attributes(dis)$Size) stop("the number of species in comm must be equal to that in dis")
  D <- as.matrix(dis)
  if(!abundance) {
    comm[comm>0] <- 1
  }
  commt <- as.data.frame(t(comm))
  
  funK <- function(v) {
    p <- v/sum(v)
    K <- apply(D, 1, function(x) sum(x*p))
    K[p<tol] <- NA
    return(K)
  }
  V <- cbind.data.frame(sapply(commt, funK))
  rownames(V) <- colnames(comm)
  colnames(V) <- rownames(comm)
  funKbar <- function(v){
    p <- v/sum(v)
    Kbar <- sapply(1:nrow(D), function(i) sum(D[i,]*v/sum(v[-i])))
    Kbar[p<tol] <- NA
    return(Kbar)
  }
  Kbar <- cbind.data.frame(sapply(commt, funKbar))
  rownames(Kbar) <- colnames(comm)
  colnames(Kbar) <- rownames(comm)
  funQ <- function(v){
    p <- v/sum(v)
    Q <- t(p)%*%D%*%p
    return(Q)
  }
  Q <- unlist(sapply(commt, funQ))
  
  funSim <- function(v){
    p <- v/sum(v)
    S <- 1-sum(p^2)
    return(S)
  } 
  Sim <- unlist(sapply(commt, funSim))
  
  funN <- function(v){
    p <- v/sum(v)
    N <- length(p[p>tol])
    return(N)
  }
  N <- unlist(sapply(commt, funN))
  U <- Q/Sim
  
  red <- cbind.data.frame(N=N, D=Sim, Q=Q, U=U)
  rownames(red) <- rownames(comm)
  
  res <- list()
  res$Kbar <- Kbar
  res$V <- V
  res$red <- red
  return(res)
  
}