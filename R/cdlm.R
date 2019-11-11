#' Helper function for fitting linear model for compound decision problems
#'
#' Calculates the design matrix X and other quantities that are used by the \code{\link{cdlm}} function
#'
#' @param Y observed data (n x 1)
#' @param s2 unbiased estimate of variance of Y (n x 1)
#' @param C main effects in the linear regression (n x p_C)
#' @param CY terms in the linear regression that interact with Y (n x p_Y)
#' @param N neighbors of each index (n x q). This encodes structural information, where the ith row of N contains indices that are expected to be similar to the ith parameter of interest.
#' @param CN terms in the linear regression that interact with Y_{i_k} for i_k in the ith row of N. This should be a list of length q, where the kth item in the list is an n x p_k matrix.
#'
#' @return
#' \item{X}{design matrix}
#' \item{Yk}{matrix of neighboring Y_k indicated by the N matrix}
#' \item{Z}{vector Z used in the empirical risk function}
#'
#' @examples
#' n = 10
#' set.seed(1)
#' theta = sort(rnorm(n))
#' s2 = abs(theta) + runif(n)
#' Y = theta + rnorm(n, sd = sqrt(s2))
#' C = cbind(1, s2)
#' CY = cbind(1, 1 / s2)
#' N = cbind(c(n, 1:(n - 1)), c(2:n, 1))
#' CN = rep(list(matrix(1, n, 1)), ncol(N))
#' XYkZ(Y, s2, C, CY, CN, N)
#'
#' @export

XYkZ = function(Y, s2,
                C = NULL, CY, CN = NULL, N = NULL) {

    n = length(Y)
    
    ## make everything matrices
    if (is.null(C)) {
        C = matrix(NA, n, 0)
    } else {
        C = as.matrix(C)
    }

    CY = as.matrix(CY)
    
    if (is.null(CN)) {
        CN = matrix(NA, n, 0)
        Yk = matrix(NA, n, 0)
    } else {
        CN = lapply(CN, as.matrix)
        pk = sapply(CN, ncol) ## number of interactions per neighbor
        CN = Reduce(cbind, CN)
        N = as.matrix(N)
        Yk = rep(list(NA), length(pk))
        for (k in 1:length(pk)) {
            Yk[[k]] = replicate(pk[k], Y[N[, k]])
        }
        Yk = Reduce(cbind, Yk)
    }
    
    ## design matrix
    X = cbind(C, CY * Y, CN * Yk)

    ## Z
    Z = c(rep(0, ncol(C)), colSums(CY * s2), rep(0, ncol(CN)))

    return(list(X = X, Yk = Yk, Z = Z))
    
}

#' Linear model for compound decision problems
#'
#' Fits linear regression model that can incorporate covariates and structural information into a compound decision rule
#'
#' @param Y observed data (n x 1)
#' @param s2 unbiased estimate of variance of Y (n x 1)
#' @param C main effects in the linear regression (n x p_C)
#' @param CY terms in the linear regression that interact with Y (n x p_Y)
#' @param N neighbors of each index (n x q). This encodes structural information, where the ith row of N contains indices that are expected to be similar to the ith parameter of interest.
#' @param CN terms in the linear regression that interact with Y_{i_k} for i_k in the ith row of N. This should be a list of length q, where the kth item in the list is an n x p_k matrix.
#' @param k3 third central moment of Y (n x 1)
#' @param k4 fourth central moment of Y (n x 1)
#' @param penalty a vector c(p, M), constrains beta of regression model to have Lp norm at most M
#'
#' @return
#' \item{b}{estimated regression coefficients}
#' \item{S}{estimated covariance matrix of sqrt(n) (estimated beta - oracle beta)}
#' \item{est}{estimated parameters of interest}
#'
#' @examples
#' n = 10
#' set.seed(1)
#' theta = sort(rnorm(n))
#' s2 = abs(theta) + runif(n)
#' Y = theta + rnorm(n, sd = sqrt(s2))
#' C = cbind(1, s2)
#' CY = cbind(1, 1 / s2)
#' N = cbind(c(n, 1:(n - 1)), c(2:n, 1))
#' CN = rep(list(matrix(1, n, 1)), ncol(N))
#'
#' ## no penalty
#' fit = cdlm(Y, s2, C, CY, CN, N, k3 = rep(0, n), k4 = 3 * s2^2)
#' mean((theta - fit$est)^2)
#' ## z-scores
#' fit$b / sqrt(diag(fit$S) / n)
#'
#' ## add penalty
#' fit = cdlm(Y, s2, C, CY, CN, N, penalty = c(1, 0.1))
#' mean((theta - fit$est)^2)
#'
#' @import CVXR
#' @export

cdlm = function(Y, s2,
                C = NULL, CY, CN = NULL, N = NULL,
                k3 = NULL, k4 = NULL,
                penalty = NULL) {

    ## calculate X, Yk, and Z
    xykz = XYkZ(Y, s2, C, CY, CN, N)
    X = xykz$X
    Yk = xykz$Yk
    Z = xykz$Z
    
    ## calculate betahat
    if (is.null(penalty)) {
        ## unconstrained
        
        XX = crossprod(X, X)
        b = solve(XX, crossprod(X, Y) - Z)
        
        ## also do inference if k3 and k4 are available
        if(!is.null(k3) && !is.null(k4)) {
            
            n = length(Y)
            
            ## make everything matrices
            if (is.null(C)) {
                C = matrix(NA, n, 0)
            } else {
                C = as.matrix(C)
            }
            
            CY = as.matrix(CY)
            
            if (is.null(CN)) {
                CN = matrix(NA, n, 0)
                I = matrix(NA, 0, 0)
            } else {
                
                ## calculate I first, then cbind CN
                CN = lapply(CN, as.matrix)
                N = as.matrix(N)
                pk = sapply(CN, ncol)
                
                I1 = matrix(0, sum(pk), sum(pk))
                for (k1 in 1:length(pk)) {
                    for (k2 in 1:length(pk)) {
                        
                        ## get k1th neighbors of all i
                        l = N[, k1]
                        ## which have k2th neighbor equal to i
                        i = which(N[l, k2] == 1:n)
                        
                        if (length(i) > 0) {
                            ## pairs of i and l such that i_k1 = l and l_k2 = i
                            pairs = cbind(i, N[i, k1])
                            ## print(c(k1, k2))
                            ## print(pairs)
                            rows = ifelse(k1 == 1, 1, cumsum(pk)[k1 - 1] + 1) : cumsum(pk)[k1]
                            cols = ifelse(k2 == 1, 1, cumsum(pk)[k2 - 1] + 1) : cumsum(pk)[k2]
                            I1[rows, cols] = I1[rows, cols] +
                                crossprod(CN[[k1]][pairs[,1],] * s2[pairs[,1]],
                                          CN[[k2]][pairs[,2],] * s2[pairs[,2]])
                        }
                        
                    }
                }
                I1 = I1 / n
                
                CN = Reduce(cbind, CN)
                I = 1 / n * crossprod(CN * Yk, CN * Yk * s2) + I1
                
            }
            
            A = 1 / n * crossprod(C, C * s2)
            B = 1 / n * crossprod(C, CY * (s2 * Y + k3))
            C = 1 / n * crossprod(C, CN * s2 * Yk)
            E = 1 / n * crossprod(CY, CY * (s2 * Y^2 + 2 * Y * k3 + k4 - 2 * s2^2))
            F = 1 / n * crossprod(CY * (s2 * Y + k3), CN * Yk)
            
            V = rbind(Reduce(cbind, list(A, B, C)),
                      Reduce(cbind, list(t(B), E, F)),
                      Reduce(cbind, list(t(C), t(F), I)))

            ## covariance of sqrt(n) * (betahat - b)
            S = solve(XX / n) %*% V %*% solve(XX / n)
            
        } else {
            S = NA
        }
            
    } else  {
        ## constrained
        
        vars = Variable(ncol(X))
        objective = Minimize(2 * sum(Z * vars) + sum((Y - X %*% vars)^2))
        constraint = list(p_norm(vars, penalty[1]) <= penalty[2])
        problem = Problem(objective, constraint)
        ## result = solve(problem)

        ## directly call solver because problem is convex
        prob_data = get_problem_data(problem, solver = "ECOS")
        solver_output = ECOSolveR::ECOS_csolve(c = prob_data[["c"]],
                                               G = prob_data[["G"]],
                                               h = prob_data[["h"]],
                                               dims = prob_data[["dims"]],
                                               A = prob_data[["A"]],
                                               b = prob_data[["b"]])
        result = unpack_results(problem, "ECOS", solver_output)
        
        b = result$getValue(vars)
        S = NA
        
    }
    
    return(list(b = drop(b), S = S, est = drop(X %*% b)))
    
}


#' Cross validation for penalized linear model for compound decision problems
#'
#' Chooses tuning parameter M in the penalized version of \code{\link{cdlm}}, using cross validation.
#'
#' @param Y observed data (n x 1)
#' @param s2 unbiased estimate of variance of Y (n x 1)
#' @param C main effects in the linear regression (n x p_C)
#' @param CY terms in the linear regression that interact with Y (n x p_Y)
#' @param N neighbors of each index (n x q). This encodes structural information, where the ith row of N contains indices that are expected to be similar to the ith parameter of interest.
#' @param CN terms in the linear regression that interact with Y_{i_k} for i_k in the ith row of N. This should be a list of length q, where the kth item in the list is an n x p_k matrix.
#' @param Lp p norm
#' @param nfolds number of folds of CV
#' @param M.min minimum value of M to consider
#' @param M.max largest value of M to consider
#' @param nM number of M to fit, equally spaced between M.min and M.max inclusive
#' @param verbose if TRUE, prints a message every time a fold is fit
#'
#' @return
#' \item{Ms}{list of M values that were tried}
#' \item{M.opt}{optimal M value}
#' \item{b}{estimated regression coefficients}
#' \item{est}{estimated parameters of interest}
#'
#' @examples
#' n = 10
#' set.seed(1)
#' theta = sort(rnorm(n))
#' s2 = abs(theta) + runif(n)
#' Y = theta + rnorm(n, sd = sqrt(s2))
#' C = cbind(1, s2)
#' CY = cbind(1, 1 / s2)
#' N = cbind(c(n, 1:(n - 1)), c(2:n, 1))
#' CN = rep(list(matrix(1, n, 1)), ncol(N))
#'
#' fit = cv.cdlm(Y, s2, C, CY, CN, N, Lp = 1, nfolds = 5)
#' mean((theta - fit$est)^2)
#' 
#' @import CVXR
#' @export

cv.cdlm = function(Y, s2, C, CY, CN = NULL, N = NULL, Lp = 1, nfolds = 2,
                   M.min = 0, M.max = 5, nM = 11, verbose = FALSE) {

    n = length(Y)
    
    ## make everything matrices
    if (is.null(C)) {
        C = matrix(NA, n, 0)
    } else {
        C = as.matrix(C)
    }

    CY = as.matrix(CY)
    
    ## issues come up when doing CV with neighbors due to indexing
    ## just tack neighbors on to C matrix for purposes of CV
    if (is.null(CN)) {
        C.cv = C
    } else {
        CN.cv = lapply(CN, as.matrix)
        pk = sapply(CN.cv, ncol) ## number of interactions per neighbor
        CN.cv = Reduce(cbind, CN.cv)
        N.cv = as.matrix(N)
        Yk = rep(list(NA), length(pk))
        for (k in 1:length(pk)) {
            Yk[[k]] = replicate(pk[k], Y[N.cv[, k]])
        }
        Yk = Reduce(cbind, Yk)
        C.cv = cbind(C, CN.cv * Yk)
    }
    
    ## possible tuning parameters
    Ms = seq(M.min, M.max, length = nM)
    
    ## cross validation
    
    ind = sample(1:n)
    inc = floor(n / nfolds)
    
    errs = matrix(NA, nfolds, length(Ms))
    for (i in 1:nfolds) {
        if (verbose) cat("fold", i, "\n")
        
        test = ind[((i - 1) * inc + 1):(i * inc)]
        if (i == nfolds) {
            test = ind[((i - 1) * inc + 1):n]
        }
        
        ## fit on training
        bs = sapply(Ms, function(M) {
            tryCatch(cdlm(Y[-test], s2[-test], C.cv[-test,], CY[-test,], penalty = c(Lp, M))$b,
                     error = function(e) rep(NA, ncol(C.cv) + ncol(CY)))
        })
        
        ## testing error
        errs[i,] = colMeans((cbind(C.cv[test,], (CY * Y)[test,]) %*% bs - Y[test])^2)
    }
    
    ## refit using all M and return optimal
    if (sum(!is.na(errs)) == 0) {
        M.opt = NA
        b = NA
        est = NA
    } else {
        M.opt = Ms[which.min(colMeans(errs))]
        fit = cdlm(Y, s2, C, CY, CN, N, penalty = c(Lp, M.opt))
        b = fit$b
        est = fit$est
    }
    return(list(Ms = Ms, M.opt = M.opt, b = b, est = est))
}
