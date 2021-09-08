#' GEBV-O Strategy
#'
#' Identify parental lines based on GEBV-O strategy and simulate their offsprings.
#'
#' @param fitted.t matrix. An n*t matrix denotes the fitted values of each traits
#' of the training population. The missing value must have been already imputed.
#' @param geno.t matrix. An n*p matrix denotes the marker score matrix of the
#' training population. The markers must be coded as 1, 0, or -1 for alleles
#' AA, Aa, or aa. The missing value must have been already imputed.
#' @param marker matrix. A p*2 matrix whose first column indicates the chromosome
#' number to which a marker belongs; and second column indicates the position of
#' the marker in centi-Morgan (cM).
#' @param geno.c matrix. An nc*p matrix denotes the marker score matrix of the
#' candidate population with nc individuals and p markers. The markers must be
#' coded as 1, 0, or -1 for alleles AA, Aa, or aa. The missing value must have
#' been already imputed. If geno.c is set to be NULL, the candidate population
#' is exactly the training population.
#' @param npl integer. An integer indicates how many parental lines with the
#' top GEBV index will be chosen from each trait. If npl is set to be NULL, there
#' will be be 4 times the number of traits.
#' @param weight vector. A vector with length t indicates the weights of target
#' traits in selection index. If weight is set to be NULL, the equal weight will
#' be assigned to all the target traits.
#' @param direction vector. A vector with length t indicates the selecting
#' directions for target traits. The elements of direction are 1, or 0
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
#' The fitted value data in the input data can be obtained by the function
#' GBLUP.fit and mmer, that can be seen in the Examples shown below.
#'
#' @export
#'
#' @references
#'
#' Chung PY, Liao CT. 2020. Identification of superior parental lines for
#' biparental crossing via genomic prediction. PLoS ONE 15(12):e0243159.
#'
#' @seealso
#' \code{\link[sommer]{mmer}}
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
#' t1 <- rnorm(10,30,10)
#' t2 <- rnorm(10,10,5)
#' t3 <- NULL
#' t4 <- NULL
#' t5 <- NULL
#' geno.test <- matrix(sample(c(1, -1), 200, replace = TRUE), 10, 20)
#' marker.test <- cbind(rep(1:2, each=10), rep(seq(0, 90, 10), 2))
#' fit <- GBLUP.fit(t1, t2, t3, t4, t5, geno = geno.test)
#'
#' geno.candidate <- matrix(sample(c(1,-1), 300, replace = TRUE), 15, 20)
#'
#' # run and output
#' result <- simu.GEBVO(fit, geno.t = geno.test, marker = marker.test,
#' geno.c = geno.candidate, nprog = 5, nsele = 10, ngen = 5, nrep = 5)
#' result$suggested.subset
#'
#'
#'
#' # other method: use mmer to obtain the fitted value
#' \dontrun{
#' set.seed(2000)
#' t1 <- rnorm(10,30,10)
#' t2 <- rnorm(10,10,5)
#' phe <- cbind(t1, t2)
#' nt <- ncol(phe)
#' geno.test <- matrix(sample(c(1, -1), 200, replace = TRUE), 10, 20)
#' marker.test <- cbind(rep(1:2, each=10), rep(seq(0, 90, 10), 2))
#' rownames(geno.test) <- 1:nrow(geno.test)
#' id <- rownames(geno.test)
#' K0 <- geno.test%*%t(geno.test)/ncol(geno.test)
#' diag(K0) <- 1
#'
#' dat <- data.frame(id, phe)
#' fit0 <- sommer::mmer(cbind(t1, t2)~1,
#'       random = ~sommer::vs(id, Gu = K0, Gtc = sommer::unsm(nt)),
#'       rcov = ~sommer::vs(units, Gtc = sommer::unsm(nt)),
#'       data = dat,
#'       tolparinv = 0.01)
#'
#' u0 <- fit0$U$`u:id`
#' fit <- matrix(unlist(u0), ncol = nt)
#' colnames(fit) <- names(u0)
#'
#' fit <- fit+matrix(fit0$fitted[1,], nrow(fit), nt, byrow = TRUE)
#' fit <- fit[order(as.numeric(names((u0[[1]])))),]
#' }
simu.GEBVO <- function(fitted.t, geno.t, marker, geno.c = NULL, npl = NULL, weight = NULL,
                       direction = NULL, nprog = 50, nsele = NULL, ngen = 10, nrep = 30, console = TRUE){

  phe.t <- fitted.t
  datatry <- try(phe.t*phe.t, silent=TRUE)
  if(class(datatry)[1] == "try-error" | NA %in% phe.t){
    stop("Phenotype data error, please cheak your phenotype data and impute the missing value.", call. = FALSE)
  }

  nt <- ncol(phe.t)
  ind.t <- nrow(phe.t)
  if(is.null(npl)){npl <- 4*nt}
  if(is.null(weight)){weight <- rep(1, nt)/nt}
  if(is.null(direction)){direction <- as.numeric(weight >= 0)}
  if(!console[1] %in% c(0,1) | length(console) > 1){console <- TRUE}

  markertest <- c(nrow(marker) != ncol(geno.t), NA%in%marker[,2], marker[,1] != sort(marker[,1]))
  datatry <- try(marker[,2]%*%marker[,2], silent = TRUE)
  if(class(datatry)[1] == "try-error" | TRUE %in% markertest){
    stop("Marker data error, please cheak your marker data. Or the number of marker does not match the genetype data.", call. = FALSE)
  }

  datatry <- try(npl*weight*direction*nprog*nsele*ngen*nrep, silent = TRUE)
  if(class(datatry)[1] == "try-error" | NA %in% datatry){
    stop("Argument error, please cheak your argument.", call. = FALSE)
  }
  direction <- as.numeric(direction > 0)

  nt <- ncol(phe.t)
  nf1 <- choose(npl, 2)
  if(is.null(nsele)){nsele <- nf1}
  if(length(weight) > nt){weight <- weight[1:nt]}
  if(length(weight) < nt){weight <- c(weight,rep(0, nt-length(weight)))}
  if(length(direction) > nt){direction <- direction[1:nt]}
  if(length(direction) < nt){direction <- c(direction, rep(1, nt-length(direction)))}

  if(nprog*nf1 < nsele){
    stop("Argument error, 'nprog' too small or 'nsele' too large.", call. = FALSE)
  }

  phe1 <- phe.sd(phe.t)
  mu0 <- phe1[[2]]
  sd0 <- phe1[[3]]
  fit <- phe1[[1]]

  row.names(fit) <- 1:nrow(fitted.t)
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
    p.c <- fit
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

  inde <- p.c%*%weight
  p0 <- rownames(inde)[order(inde, decreasing = TRUE)][1:npl]

  phe.m0 <- matrix(0, length(p0), nrow(p.c))
  for(i in 1:length(p0)){
    phe.m0[i, p0[i]==rownames(p.c)] <- 1
  }
  phe.p0 <- phe.m0%*%p.c
  rownames(phe.p0) <- p0
  geno.p0 <- geno.c[rownames(p.c) %in% p0,]

  kp0 <- geno.p0%*%t(geno.p0)/ncol(geno.p0)
  diag(kp0) <- 1
  Dscore <- rep(det(kp0),nrep)

  GEBVO.p0 <- cbind(1:nrow(phe.p0), phe.p0, geno.p0)

  gvalue.result <- list()
  p.result <- matrix(0, nrep, npl)

  for(m in 1:nrep){
    nf1 <- choose(npl, 2)
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
    for(i in 1:(npl-1)){
      for(j in (i+1):npl){
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

