#' Muti-trait GBLUP Model
#'
#' Built the muti-trait GBLUP model using the phenotypic and genotypic data of a
#' training population by 'mmer' from R package 'sommer'. Then, output the fitted
#' values of the training population.
#'
#' @param t1 vector. The phenotype of trait1. The missing value must be coded as NA.
#' The length of all triat must be the same.
#' @param t2 vector. The phenotype of trait2. The missing value must be coded as NA.
#' The length of all triat must be the same.
#' @param t3 vector. The phenotype of trait3. The missing value must be coded as NA.
#' The length of all triat must be the same.
#' @param t4 vector. The phenotype of trait4. The missing value must be coded as NA.
#' The length of all triat must be the same.
#' @param t5 vector. The phenotype of trait5. The missing value must be coded as NA.
#' The length of all triat must be the same.
#' @param geno matrix. An n*p matrix with n individuals and p markers of the
#' training population. The markers must be coded as 1, 0, or -1 for alleles AA,
#' Aa, or aa. The missing value must have been already imputed.
#' @param K matrix. An n*n matrix denotes the genomic relationship matrix of the
#' training population if geno is set to be NULL.
#'
#' @return
#' The fitted values of the training population.
#'
#' @export
#'
#' @seealso
#' \code{\link[sommer]{mmer}}
#'
#' @references
#'
#' Habier D, Fernando RL, Dekkers JCM. 2007. The impact of genetic relationship
#' information on genome-assisted breeding values. Genetics 177:2389-2397.
#'
#' VanRaden PM. 2008. Efficient methods to compute genomic predictions.
#' J Dairy Sci. 91:4414-4423.
#'
#' @examples
#' # generate simulated data
#' t1 = rnorm(50,30,10)
#' t2 = rnorm(50,10,5)
#' t3 = rnorm(50,20,20)
#' t4 = NULL
#' t5 = NULL
#'
#' # run with the marker score matrix
#' geno.test <- matrix(sample(c(1, -1), 5000, replace = TRUE), 50, 100)
#' result1 <- GBLUP.fit(t1, t2, t3, t4, t5, geno = geno.test)
#' result1
#'
#' # run with the genomic relationship matrix
#' K.test <- geno.test%*%t(geno.test)/ncol(geno.test)
#' result2 <- GBLUP.fit(t1, t2, t3, t4, t5, K = K.test)
#' result2
GBLUP.fit <- function(t1, t2, t3, t4, t5, geno = NULL, K = NULL){

  phetest <- c(length(t1), length(t2), length(t3), length(t4), length(t5))
  if(length(table(phetest[phetest != 0])) != 1){
    stop("Phenotype data error, please input number vectors with the same length.", call. = FALSE)
  }

  phe <- cbind(t1, t2, t3, t4, t5)
  datatry <- try(phe*phe, silent=TRUE)
  if(class(datatry)[1] == "try-error"){
    stop("Phenotype data error, please input number vectors with the same length.", call. = FALSE)
  }
  nt <- ncol(phe)

  if(is.null(geno) & is.null(K)){
    stop("One of the arguments 'geno' and 'K' must be assigned.", call. = FALSE)
  }

  if(is.null(geno)){
    K0 <- K
    rownames(K0) <- 1:nrow(K0)
    colnames(K0) <- 1:nrow(K0)
    id <- rownames(K0)
  } else {
    datatry <- try(geno%*%t(geno), silent = TRUE)
    if(class(datatry)[1] == "try-error" | length(geno[geno != 1 & geno != 0 & geno != -1]) > 0){
      stop("Genotype data error or have not been imputed.", call. = FALSE)
    }
    rownames(geno) <- 1:nrow(geno)
    id <- rownames(geno)
    K0 <- geno%*%t(geno)/ncol(geno)
    diag(K0) <- 1
  }

  datatry <- try(K0%*%as.matrix(phe), silent = TRUE)
  if(class(datatry)[1] == "try-error" | NA%in%K0){
    stop("Input data error, please cheak your input data.", call. = FALSE)
  }

  dat <- data.frame(id,phe)
  if(nt == 1){
    t1 <- c(phe)
    fit <- sommer::mmer(t1~1,
                        random = ~sommer::vs(id, Gu = K0, Gtc = sommer::unsm(nt)),
                        rcov = ~sommer::vs(units, Gtc = sommer::unsm(nt)),
                        data = dat)

    tol <- 10^-5
    while(length(fit) == 0){
      cat("Try bigger tolparinv. tolparinv:", tol, "\n")
      fit <- sommer::mmer(t1~1,
                          random = ~sommer::vs(id, Gu = K0, Gtc = sommer::unsm(nt)),
                          rcov = ~sommer::vs(units, Gtc = sommer::unsm(nt)),
                          data = dat,
                          tolparinv = tol)
      tol <- tol*10
    }
  } else {
    fit <- sommer::mmer(cbind(t1, t2, t3, t4, t5)~1,
                        random = ~sommer::vs(id, Gu = K0, Gtc = sommer::unsm(nt)),
                        rcov = ~sommer::vs(units, Gtc = sommer::unsm(nt)),
                        data = dat)

    tol <- 10^-5
    while(length(fit) == 0){
      cat("Try bigger tolparinv. tolparinv:", tol, "\n")
      fit <- sommer::mmer(cbind(t1, t2, t3, t4, t5)~1,
                          random = ~sommer::vs(id, Gu = K0, Gtc = sommer::unsm(nt)),
                          rcov = ~sommer::vs(units, Gtc = sommer::unsm(nt)),
                          data = dat,
                          tolparinv = tol)
      tol <- tol*10
    }
  }

  u0 <- fit$U$`u:id`
  fitted.value <- c()
  for(i in 1:length(u0)){
    fitted.value <- cbind(fitted.value,u0[[i]])
  }
  colnames(fitted.value) <- names(u0)

  fitted.value <- fitted.value+matrix(fit$fitted[1,], nrow(fitted.value), nt, byrow = TRUE)
  fitted.value <- fitted.value[order(as.numeric(rownames((fitted.value)))),]
  if(!is.null(row.names(phe))){row.names(fitted.value) <- row.names(phe)}
  fitted.value <- as.matrix(fitted.value)

  return(fitted.value)
}

