#' GEBV-O Strategy
#'
#' Identify parental lines based on GEBV-O strategy and simulate their offsprings.
#'
#' @param phe.t matrix. An n*t matrix denotes the phenotypic values of the
#' training population with n individuals and t target traits.
#' @param geno.t matrix. An n*p matrix denotes the marker score matrix of the
#' training population. The markers must be coded as 1, 0, or -1 for alleles
#' AA, Aa, or aa. The missing value must have been already imputed.
#' @param marker matrix. An p*2 matrix whose first column indicates the chromosome
#' number to which a marker belongs; and second column indicates the position of
#' the marker in centi-Morgan (cM).
#' @param geno.c matrix. An nc*p matrix denotes the marker score matrix of the
#' candidate population with nc individuals and p markers. The markers must be
#' coded as 1, 0, or -1 for alleles AA, Aa, or aa. The missing value must have
#' been already imputed. If geno.c is set to be NULL, the candidate population
#' is exactly the training population.
#' @param npl vector. A vector denotes how many parental lines with the top GEBVs
#' will be chosen from each trait. If npl is set to be NULL, there will be 4
#' parental lines chosen from each trait.
#' @param weight vector. A vector with length t indicates the weights of target
#' traits in selection index. If weight is set to be NULL, the equal weight will
#' be assigned to all the target traits.
#' @param direction vector. A vector with length t indicates the selecting
#' directions for target traits. The elements of direction are 1, or -1
#' representing the rule that the larger the better; or the smaller the better.
#' If direction is set to be NULL, the selecting direction will be the same as
#' weight.
#' @param nprog integer. An integer indicates the number of progenies which
#' will be produced for each of the best individuals at every generation.
#' @param nsele integer. An integer indicates the number of the best individuals
#' which will be selected at each generation. If nsele is set to be NULL, the
#' number will be the same as the number of F1 individuals.
#' @param ngen integer. An integer indicates the number of generations in the
#' simulation process.
#' @param nrep integer. An integer indicates the number of repetitions in the
#' simulation process.
#' @param console logical. A logical variable, if console is set to be TRUE,
#' the simulation process will be shown in the R console.
#'
#' @return
#' \item{method}{The GEBV-O strategy.}
#' \item{weight}{The weights of target traits in selection index.}
#' \item{direction}{The selecting directions of target traits in selection index.}
#' \item{mu}{The mean vector of target traits.}
#' \item{sd}{The standard deviation vector of target traits.}
#' \item{GEBV.value}{The GEBVs of target traits in each generation and each
#' repetition.}
#' \item{parental.lines}{The IDs and D-score of parental lines selected in
#' each repetition.}
#' \item{suggested.subset}{The most frequently selected parental lines by this
#' strategy.}
#'
#' @note
#' The function output.best and output.gain can be used to summarize the result.
#'
#' @export
#'
#' @references
#'
#' Chung PY, Liao CT. 2020. Identification of superior parental lines for
#' biparental crossing via genomic prediction. PLoS ONE 15(12):e0243159.
#'
#' @seealso
#' \code{\link[IPLGP]{phe.sd}}
#' \code{\link[IPLGP]{GBLUP.fit}}
#' \code{\link[IPLGP]{GA.Dscore}}
#' \code{\link[IPLGP]{simu.gamete}}
#' \code{\link[IPLGP]{simu.GDO}}
#' \code{\link[IPLGP]{simu.GEBVGD}}
#' \code{\link[IPLGP]{output.best}}
#' \code{\link[IPLGP]{output.gain}}
#'
#' @examples
#' # generate simulated data
#' set.seed(2000)
#' phe.test <- data.frame(trait1 = rnorm(10,30,10), trait2 = rnorm(10,10,5))
#' geno.test <- matrix(sample(c(1, -1), 200, replace = TRUE), 10, 20)
#' marker.test <- cbind(rep(1:2, each=10), rep(seq(0, 90, 10), 2))
#' geno.candidate <- matrix(sample(c(1,-1), 300, replace = TRUE), 15, 20)
#'
#' # run and output
#' result <- simu.GEBVO(phe.test, geno.test, marker.test, geno.candidate,
#' nprog = 5, nsele = 10, ngen = 5, nrep = 5)
#' result$suggested.subset
simu.GEBVO <- function(phe.t, geno.t, marker, geno.c = NULL, npl = NULL, weight = NULL,
                     direction = NULL, nprog = 50, nsele = NULL, ngen = 10, nrep = 30,
                     console = TRUE){

  nt <- ncol(phe.t)
  ind.t <- nrow(phe.t)
  if(is.null(npl)){npl <- rep(4, nt)}
  if(is.null(weight)){weight <- rep(1, nt)/nt}
  if(is.null(direction)){direction <- weight/abs(weight)}
  if(console != TRUE & console != FALSE){console <- TRUE}

  n0 <- sum(npl)
  nt <- length(npl)
  nf1 <- choose(n0, 2)
  if(is.null(nsele)){nsele <- nf1}
  if(length(weight) > nt){weight <- weight[1:nt]}
  if(length(weight) < nt){weight <- c(weight,rep(0, nt-length(weight)))}
  if(length(direction) > nt){direction <- direction[1:nt]}
  if(length(direction) < nt){direction <- c(direction, rep(1, nt-length(direction)))}

  markertest <- c(nrow(marker) != ncol(geno.t), NA%in%marker[,2], marker[,1] != sort(marker[,1]))
  datatry <- try(marker[,2]%*%marker[,2], silent=TRUE)
  if(class(datatry)[1] == "try-error" | TRUE%in%markertest){
    stop("Marker data error, please cheak your marker data. Or the number of marker does not match the genetype data.", call. = FALSE)
  }

  datatry <- try(npl*weight*direction*nprog*nsele*ngen*nrep, silent = TRUE)
  if(class(datatry)[1] == "try-error" | NA%in%datatry){
    stop("Argument error, please cheak your argument.", call. = FALSE)
  }

  if(nprog*nf1 < nsele){
    stop("Argument error, 'nprog' too small or 'nsele' too large.", call. = FALSE)
  }

  phe1 <- phe.sd(phe.t)
  mu0 <- phe1[[2]]
  sd0 <- phe1[[3]]
  phe2 <- phe1[[1]]

  fit <- GBLUP.fit(phe = phe2, geno = geno.t)
  row.names(fit) <- 1:nrow(fit)
  K0 <- geno.t%*%t(geno.t)/ncol(geno.t)
  diag(K0) <- 1
  d0 <- det(K0)
  i0 <- 10^-10
  while(d0 == 0){
    i0 <- i0*10
    d0 <- det(K0+diag(i0, nrow(K0)))
  }
  K00 <- solve(K0+diag(i0, nrow(K0)))

  if(is.null(geno.c)){
    geno.c <- geno.t
    p.c < -fit
  } else {
    datatry <- try(geno.t%*%t(geno.c), silent=TRUE)
    if(class(datatry)[1] == "try-error" | length(geno.c[geno.c != 1 & geno.c != 0 & geno.c != -1]) > 0){
      stop("Candidate set genotype data error, please cheak your candidate set genotype data.", call. = FALSE)
    }

    geno.c2 <- rbind(geno.t,geno.c)
    Kptp <- t(geno.t%*%t(geno.c2))/ncol(geno.c2)
    p.c <- Kptp%*%K00%*%fit
    npc <- nrow(p.c)-ind.t
    p.c <- cbind(matrix(0, npc, ind.t), diag(npc))%*%p.c
  }

  if(is.null(row.names(geno.c))){
    row.names(p.c) <- 1:nrow(p.c)
    row.names(geno.c) <- row.names(p.c)
  } else {row.names(p.c) <- row.names(geno.c)}

  p0 <- c()
  direction0 <- direction > 0
  k <- 1
  while(length(p0) < n0){
    for(ps in 1:nt){
      p1 <- rownames(p.c)[order(p.c[,ps], decreasing = direction0[ps])[1:npl[ps]]]
      p0 <- union(p0, p1)
    }
    npl[k] <- npl[k]+1
    k <- k+1
    if(k > nt){k <- 1}
  }
  p0 <- p0[1:n0]

  phe.m0 <- matrix(0, length(p0), nrow(p.c))
  for(i in 1:length(p0)){
    phe.m0[i, p0[i]==rownames(p.c)] <- 1
  }
  phe.p0 <- phe.m0%*%p.c
  rownames(phe.p0) <- p0
  geno.p0 <- geno.c[rownames(p.c)%in%p0,]

  kp0 <- geno.p0%*%t(geno.p0)/ncol(geno.p0)
  diag(kp0) <- 1
  Dscore <- rep(det(kp0),nrep)

  GEBVO.p0 <- cbind(1:nrow(phe.p0), phe.p0, geno.p0)

  gvalue.result <- list()
  p.result <- matrix(0, nrep, n0)

  for(m in 1:nrep){
    nf1 <- choose(n0, 2)
    GEBVO.gvalue <- list()
    GEBVO.SNP <- list()

    p.result[m,] <- rownames(GEBVO.p0)

    GEBV0 <- phe.p0
    if(is.null(colnames(phe.t))){
      colnames(GEBV0)<-paste("t", 1:nt, sep = "")
    } else {colnames(GEBV0) <- colnames(phe.t)}
    GEBVO.gvalue[[1]] <- GEBV0*matrix(sd0, nrow(GEBVO.p0), nt, byrow = TRUE)+matrix(mu0, nrow(GEBVO.p0), nt, byrow = TRUE)
    GEBVO.p <- GEBVO.p0[,(nt+2):ncol(GEBVO.p0)]

    GEBVO.SNP[[1]] <- GEBVO.p
    GEBVO.F1 <- list()
    k1 <- 1
    for(i in 1:(n0-1)){
      for(j in (i+1):n0){
        GEBVO.F1[[k1]] <- as.matrix(GEBVO.p[c(i, j),])
        k1 <- k1+1
      }
    }

    GEBVO.F1.SNP <- c()
    for(i in 1:nf1){
      F1_0 <- (GEBVO.F1[[i]][1,]+GEBVO.F1[[i]][2,])/2
      GEBVO.F1.SNP <- cbind(GEBVO.F1.SNP, F1_0)
    }
    GEBVO.SNP[[2]] <- GEBVO.F1.SNP
    GEBVO.F1.SNP2 <- rbind(geno.t, t(GEBVO.F1.SNP))

    Kpt <- t(geno.t%*%t(GEBVO.F1.SNP2))/nrow(GEBVO.F1.SNP)

    fity.F1 <- Kpt%*%K00%*%fit
    fity.F1 <- cbind(matrix(0, nf1, ind.t), diag(nf1))%*%fity.F1
    row.names(fity.F1) <- NULL
    if(is.null(colnames(phe.t))){
      colnames(fity.F1) <- paste("t",1:nt,sep="")
    } else {colnames(fity.F1) <- colnames(phe.t)}
    GEBVO.gvalue[[2]] <- fity.F1*matrix(sd0, nrow(fity.F1), nt, byrow = TRUE)+matrix(mu0, nrow(fity.F1), nt, byrow = TRUE)

    if(console){
      cat("Method", "Repeat", "Generation", "\n")
      cat("GEBVO", m, paste("F", 1, sep = ""), "\n", sep = "\t")
    }

    p.snp <- GEBVO.F1

    for(g in 3:(ngen+1)){
      GEBVO.F2 <- list()
      k2 <- 1
      for(i in 1:nf1){
        for(j in 1:nprog){
          marker.ind <- cbind(marker, t(p.snp[[i]]))
          F2_1 <- simu.gamete(marker.ind)
          F2_2 <- simu.gamete(marker.ind)
          F2 <- cbind(F2_1, F2_2)
          GEBVO.F2[[k2]] <- t(F2)
          k2 <- k2+1
        }
      }
      GEBVO.F2.SNP <- matrix(0, length(GEBVO.F2[[1]][1,]), (nf1*nprog))
      for(i in 1:(nf1*nprog)){
        F2_0 <- (GEBVO.F2[[i]][1,]+GEBVO.F2[[i]][2,])/2
        GEBVO.F2.SNP[,i] <- F2_0
      }
      GEBVO.SNP[[g]] <- GEBVO.F2.SNP
      GEBVO.F2.SNP2 <- rbind(geno.t, t(GEBVO.F2.SNP))

      Kpt2 <- t(geno.t%*%t(GEBVO.F2.SNP2))/nrow(GEBVO.F2.SNP)

      fity.F2 <- Kpt2%*%K00%*%fit
      nf2 <- nrow(fity.F2)-ind.t
      fity.F2 <- cbind(matrix(0, nf2, ind.t), diag(nf2))%*%fity.F2
      row.names(fity.F2) <- NULL
      if(is.null(colnames(phe.t))){
        colnames(fity.F2) <- paste("t", 1:nt, sep = "")
      } else {colnames(fity.F2) <- colnames(phe.t)}
      GEBVO.gvalue[[g]] <- fity.F2*matrix(sd0, nrow(fity.F2), nt, byrow = TRUE)+matrix(mu0, nrow(fity.F2), nt, byrow = TRUE)

      if(console){cat("GEBVO", m, paste("F", g-1, sep = ""), "\n", sep = "\t")}

      sele <- order(fity.F2%*%weight, decreasing = TRUE)
      p.snp <- list()
      nf1 <- nsele
      for(i in 1:nf1){
        p.snp[[i]] <- GEBVO.F2[[sele[i]]]
      }
    }
    names(GEBVO.gvalue) <- c("p", paste("F", 1:ngen, sep = ""))
    gvalue.result[[m]] <- GEBVO.gvalue
  }

  bestsub <- table(p.result)
  bestsub <- data.frame(parental.lines = names(bestsub), chosen.ratio = as.numeric(bestsub)/nrep)
  parental.lines <- list(parental.lines = p.result, D.score = Dscore)

  return(list(method = "GEBV-O", weight = weight, direction = direction, mu = mu0, sd = sd0, GEBV.value = gvalue.result, parental.lines = parental.lines, suggested.subset = bestsub))
}

