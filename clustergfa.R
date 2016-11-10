countC <<- 0
ALPHA <<- 100
GFA_cluster <- function(Y,R,K,S,opts) {
    cat("GFA_cluster starting...\n")
    #
    # Store dimensionalities of data sets 
    #
    M <- length(Y)              # The number of views; M=2 corresponds to BIBFA
    D <- unlist(lapply(Y,ncol)) # Collect the number of features in vector D
    Ds <- sum(D)                # Total number of features
    N <- nrow(Y[[1]])

    alpha_0 <- opts$prior.alpha_0   # Easier access for hyperprior values
    beta_0 <- opts$prior.beta_0
    alpha_0t <- opts$prior.alpha_0t
    beta_0t <- opts$prior.beta_0t
    alpha_0g <- opts$alpha_0g
    beta_0g  <- opts$beta_0g

    ## Initialize the model parameters:

    # Latent variables
    #Z <- matrix(rnorm(1*K,0,1),1,K) # The mean 
    Z <- zz <- covZ <- ZZ <- chol_covZinv <- vector("list",length=S)               # The covariance

    Zn <- vector("list",length=S) #matrix(rnorm(N*K),N,K) * 0 #TEST
    ZZn <- vector("list",length=S)
    covZn <- vector("list",length=S)


    tau <- a_tau <- b_tau <- vector("list",length=S)  # The mean noise precisions
    logtau <- logalpha <- vector("list",length=S)
    alpha <- a_ard <- b_ard <- vector("list",length=S) #ARD
    for (s in 1:S) {
        # ARD parameters
        alpha[[s]] <- matrix(1,M,K)          # The mean of the ARD precisions
        b_ard[[s]] <- matrix(1,M,K)          # The parameters of the Gamma distribution
        a_ard[[s]] <- alpha_0 + D/2               # for ARD precisions
        # tau noise precisions
        tau[[s]] <- rep(opts$init.tau,M)
        logtau[[s]] <- log(tau[[s]])
        a_tau[[s]] <- alpha_0t               # The parameters of the Gamma distribution
        b_tau[[s]] <- rep(0,M)               #     for the noise precisions
        # ZZ^T
        Z[[s]]    <- matrix(rnorm(N*K),N,K)
        Zn[[s]]    <- matrix(rnorm(N*K),N,K)  #TEST
        covZ[[s]] <- diag(1,K)
        covZn[[s]]<- diag(1,K) #TEST
        ZZ[[s]] <- crossprod(Z[[s]],Z[[s]]/S)+N/S*covZ[[s]]
        ZZn[[s]]<- crossprod(Zn[[s]],Zn[[s]]/S)+N/S*covZn[[s]]
    }
    a_ardn <- alpha_0 + D/2
    b_ardn <- matrix(1,M,K)
    alphan <- matrix(1,M,K)
    logalphan <- digamma(matrix(a_ardn,M,K,byrow=T)) - log(b_ardn)

    # The projections
    W <- vector("list",length=S)    # The means
    covW <- vector("list",length=S) # The covariances
    WW <- vector("list",length=S)   # The second moments
    for (s in 1:S) for(m in 1:M) {
      W[[s]][[m]] <- matrix(rnorm(D[m]*K),D[m],K)
      covW[[s]][[m]] <- diag(1,K)
      WW[[s]][[m]] <- crossprod(W[[s]][[m]]) + covW[[s]][[m]]*D[m]
    }
    Wn <- covWn <- WWn <- list()
    for (m in 1:M) {
        Wn[[m]] <- matrix(rnorm(D[m]*K),D[m],K)  #TEST
        covWn[[m]] <- diag(1,K)
        WWn[[m]] <- crossprod(Wn[[m]])+diag(1,K)
    }

    # initializa q(pi|alpha_pi) = Dir(pi|alpha_pi)
    alpha_pi <- rep(1,S) # TODO: put prior in opts
    PI <- alpha_pi / sum(alpha_pi)
    log_pi <- digamma(alpha_pi)-digamma(sum(alpha_pi))

    # init gamma q(gamma_k|a_gamma,b_gamma) ~ Beta(a_g,b_g)
    gamma <- rep(.5,S)
    a_gamma <- b_gamma <- rep(1,S)
    log_gamma <- log_1gamma <- digamma(a_gamma) - digamma(a_gamma+b_gamma)

    # init C
    C <- matrix(1/S,N,S)

    ## Main VB loop:
    cost <- c()
    for (iter in 1:opts$iter.max) {
        if (iter>1) lb1 <- calc_lb(Y,R,C,Z,covZ,ZZ,Zn,covZn,ZZn,logtau,logalpha,logalphan,tau,a_tau,b_tau,a_ard,b_ard,a_ardn,b_ardn,W,WW,covW,Wn,WWn,covWn,alpha,alphan,alpha_0,beta_0,alpha_0t,beta_0t,alpha_0g,beta_0g,gamma,log_gamma,log_1gamma,a_gamma,b_gamma,PI,log_pi,alpha_pi)
        # update cluster assignments C
        res <- up_C(Y,Z,covZ,Zn,covZn,W,WW,Wn,WWn,logtau,tau,
                    R,log_gamma,log_1gamma,log_pi) 
        tmp_c <- C
        C <- res$C
        if (iter>1) lb2 <- calc_lb(Y,R,C,Z,covZ,ZZ,Zn,covZn,ZZn,logtau,logalpha,logalphan,tau,a_tau,b_tau,a_ard,b_ard,a_ardn,b_ardn,W,WW,covW,Wn,WWn,covWn,alpha,alphan,alpha_0,beta_0,alpha_0t,beta_0t,alpha_0g,beta_0g,gamma,log_gamma,log_1gamma,a_gamma,b_gamma,PI,log_pi,alpha_pi)
        if (iter>1 && lb2 < lb1) {print("ERROR in updating C");print(lb2-lb1);countC <<- countC+1}
        # update cluster specific latents Z
        if (iter>1) lb1 <- calc_lb(Y,R,C,Z,covZ,ZZ,Zn,covZn,ZZn,logtau,logalpha,logalphan,tau,a_tau,b_tau,a_ard,b_ard,a_ardn,b_ardn,W,WW,covW,Wn,WWn,covWn,alpha,alphan,alpha_0,beta_0,alpha_0t,beta_0t,alpha_0g,beta_0g,gamma,log_gamma,log_1gamma,a_gamma,b_gamma,PI,log_pi,alpha_pi)
        for (s in 1:S) {
            res       <- up_Z(Y,W[[s]],WW[[s]],tau[[s]],Zn[[s]],Wn,C[,s])
            Z[[s]]    <- res$Z
            ZZ[[s]]   <- res$ZZ
            covZ[[s]] <- res$covZ
        }
        if (iter>1) lb2 <- calc_lb(Y,R,C,Z,covZ,ZZ,Zn,covZn,ZZn,logtau,logalpha,logalphan,tau,a_tau,b_tau,a_ard,b_ard,a_ardn,b_ardn,W,WW,covW,Wn,WWn,covWn,alpha,alphan,alpha_0,beta_0,alpha_0t,beta_0t,alpha_0g,beta_0g,gamma,log_gamma,log_1gamma,a_gamma,b_gamma,PI,log_pi,alpha_pi)
        if (iter>1 && lb2 < lb1) {print("ERROR in updating Z");print(lb2-lb1);}
        # update noise latents Zn 
        if (iter>1) lb1 <- calc_lb(Y,R,C,Z,covZ,ZZ,Zn,covZn,ZZn,logtau,logalpha,logalphan,tau,a_tau,b_tau,a_ard,b_ard,a_ardn,b_ardn,W,WW,covW,Wn,WWn,covWn,alpha,alphan,alpha_0,beta_0,alpha_0t,beta_0t,alpha_0g,beta_0g,gamma,log_gamma,log_1gamma,a_gamma,b_gamma,PI,log_pi,alpha_pi)
        for (s in 1:S) {
            res   <- up_Z(Y,Wn,WWn,tau[[s]],Z[[s]],W[[s]],C[,s])
            Zn[[s]] <- res$Z
            ZZn[[s]]   <- res$ZZ
            covZn[[s]] <- res$covZ
        }
        if (iter>1) lb2 <- calc_lb(Y,R,C,Z,covZ,ZZ,Zn,covZn,ZZn,logtau,logalpha,logalphan,tau,a_tau,b_tau,a_ard,b_ard,a_ardn,b_ardn,W,WW,covW,Wn,WWn,covWn,alpha,alphan,alpha_0,beta_0,alpha_0t,beta_0t,alpha_0g,beta_0g,gamma,log_gamma,log_1gamma,a_gamma,b_gamma,PI,log_pi,alpha_pi)
        if (iter>1 && lb2 < lb1) {print("ERROR in updating Zn");print(lb2-lb1);}
        # update estimated labels R
        res   <- up_R(C,gamma)
        Rest  <- res$R
        # update cluster "labeling" gamma
        if (iter>1) lb1 <- calc_lb(Y,R,C,Z,covZ,ZZ,Zn,covZn,ZZn,logtau,logalpha,logalphan,tau,a_tau,b_tau,a_ard,b_ard,a_ardn,b_ardn,W,WW,covW,Wn,WWn,covWn,alpha,alphan,alpha_0,beta_0,alpha_0t,beta_0t,alpha_0g,beta_0g,gamma,log_gamma,log_1gamma,a_gamma,b_gamma,PI,log_pi,alpha_pi)
        res       <- up_gamma(alpha_0g,beta_0g,R,C)
        a_gamma   <- res$a_gamma
        b_gamma   <- res$b_gamma
        gamma     <- res$gamma # E[gamma]
        log_gamma <- res$log_gamma # E[log(gamma)]
        log_1gamma<- res$log_1gamma # E[log(1-gamma)]
        if (iter>1) lb2 <- calc_lb(Y,R,C,Z,covZ,ZZ,Zn,covZn,ZZn,logtau,logalpha,logalphan,tau,a_tau,b_tau,a_ard,b_ard,a_ardn,b_ardn,W,WW,covW,Wn,WWn,covWn,alpha,alphan,alpha_0,beta_0,alpha_0t,beta_0t,alpha_0g,beta_0g,gamma,log_gamma,log_1gamma,a_gamma,b_gamma,PI,log_pi,alpha_pi)
        if (iter>1 && lb2 < lb1) {print(paste("ERROR in updating gamma:",lb2,lb1,lb2-lb1));}
        # update cluster prior probabilities PI
        if (iter>1) lb1 <- calc_lb(Y,R,C,Z,covZ,ZZ,Zn,covZn,ZZn,logtau,logalpha,logalphan,tau,a_tau,b_tau,a_ard,b_ard,a_ardn,b_ardn,W,WW,covW,Wn,WWn,covWn,alpha,alphan,alpha_0,beta_0,alpha_0t,beta_0t,alpha_0g,beta_0g,gamma,log_gamma,log_1gamma,a_gamma,b_gamma,PI,log_pi,alpha_pi)
        res   <- up_pi(C,alpha_pi,log_pi)
        PI    <- res$PI
        alpha_pi <- res$alpha_pi
        log_pi<- res$log_pi
        if (iter>1) lb2 <- calc_lb(Y,R,C,Z,covZ,ZZ,Zn,covZn,ZZn,logtau,logalpha,logalphan,tau,a_tau,b_tau,a_ard,b_ard,a_ardn,b_ardn,W,WW,covW,Wn,WWn,covWn,alpha,alphan,alpha_0,beta_0,alpha_0t,beta_0t,alpha_0g,beta_0g,gamma,log_gamma,log_1gamma,a_gamma,b_gamma,PI,log_pi,alpha_pi)
        if (iter>1 && lb2 < lb1) {print("ERROR in updating pi");print(lb2-lb1)}
        # update factor loadings for "noise model"
        if (iter>1) lb1 <- calc_lb(Y,R,C,Z,covZ,ZZ,Zn,covZn,ZZn,logtau,logalpha,logalphan,tau,a_tau,b_tau,a_ard,b_ard,a_ardn,b_ardn,W,WW,covW,Wn,WWn,covWn,alpha,alphan,alpha_0,beta_0,alpha_0t,beta_0t,alpha_0g,beta_0g,gamma,log_gamma,log_1gamma,a_gamma,b_gamma,PI,log_pi,alpha_pi)
        res   <- up_Wn(Y,Zn,ZZn,Z,W,C,alphan,tau)
        Wn    <- res$Wn
        WWn   <- res$WWn
        covWn <- res$covWn
        if (iter>1) lb2 <- calc_lb(Y,R,C,Z,covZ,ZZ,Zn,covZn,ZZn,logtau,logalpha,logalphan,tau,a_tau,b_tau,a_ard,b_ard,a_ardn,b_ardn,W,WW,covW,Wn,WWn,covWn,alpha,alphan,alpha_0,beta_0,alpha_0t,beta_0t,alpha_0g,beta_0g,gamma,log_gamma,log_1gamma,a_gamma,b_gamma,PI,log_pi,alpha_pi)
        if (iter>1 && lb2 < lb1) {print("ERROR in updating Wn");print(lb2-lb1)}
        # update ARD alpha for "noise model" loadings
        res   <- up_alpha(Wn,WWn,alpha_0,beta_0)
        alphan    <- res$alpha
        a_ardn    <- res$a_ard
        b_ardn    <- res$b_ard
        logalphan <- res$logalpha
        # update factor loadings for each cluster
        for (s in 1:S) {
            res      <- up_W(Y,Z[[s]],ZZ[[s]],Zn[[s]],Wn,
                             alpha[[s]],tau[[s]],C[,s])
            W[[s]]   <- res$W
            WW[[s]]  <- res$WW
            covW[[s]]<- res$covW
        }
        # update ARD alpha for each cluster
        for (s in 1:S) {
            res            <- up_alpha(W[[s]],WW[[s]],
                                       alpha_0,beta_0)
            alpha[[s]]     <- res$alpha
            a_ard[[s]]     <- res$a_ard
            b_ard[[s]]     <- res$b_ard
            logalpha[[s]]  <- res$logalpha
        }
        # update tau for each cluster
        for (s in 1:S) {
            res        <- up_tau(alpha_0t,beta_0t,Y,C[,s],Zn[[s]],ZZn[[s]],
                            Wn,WWn,Z[[s]],ZZ[[s]],W[[s]],WW[[s]])
            a_tau[[s]] <- res$a_tau
            b_tau[[s]] <- res$b_tau
            tau[[s]]   <- res$tau
            logtau[[s]]<- res$logtau
        }
        cost[iter] <- calc_lb(Y,R,C,Z,covZ,ZZ,Zn,covZn,ZZn,logtau,logalpha,logalphan,tau,a_tau,b_tau,a_ard,b_ard,a_ardn,b_ardn,W,WW,covW,Wn,WWn,covWn,alpha,alphan,alpha_0,beta_0,alpha_0t,beta_0t,alpha_0g,beta_0g,gamma,log_gamma,log_1gamma,a_gamma,b_gamma,PI,log_pi,alpha_pi)
        print(paste("Iteration",iter,"lb:",cost[iter]))
        #browser()
        if(iter>1 && abs((cost[iter]-cost[iter-1])/cost[iter]) < opts$iter.crit) {
            print("Converged!")
            break
        }
        if(iter>1 && cost[iter]<cost[iter-1])
            print("ERROR IN VB ITERATION LB!!!")
    }
    # return model:
    list(cost=cost,gamma=gamma,W=W,covW=covW,ZZ=ZZ,WW=WW,Z=Z,
         covZ=covZ,tau=tau,alpha=alpha,D=D,K=K,S=S,R=R,C=C,
         Rest=Rest,alpha_pi=alpha_pi,alphan=alphan,
         Wn=Wn,WWn=WWn,PI=PI,log_pi=log_pi,log_gamma=log_gamma,
         log_1gamma=log_1gamma,a_gamma=a_gamma,b_gamma=b_gamma,
         Zn=Zn,logtau=logtau,covZn=covZn,ZZn=ZZn)
}


up_R <- function(C,gamma) {
    R <- C %*% gamma
    list(R=R)
}

up_Z <- function(Y,W,WW,tau,Zn,Wn,C) {
    ## Update latent Z
    M <- length(Y)

    # Efficient and robust way of computing
    # solve(diag(1,K) + tau * WW^t)
    covZ <- diag(1,ncol(W[[1]]))
    for(m in 1:M) {
      covZ <- covZ + tau[m]*WW[[m]]
    }
    covZ <- chol2inv(chol(covZ))

    Z <- (Y[[1]]-tcrossprod(Zn,Wn[[1]]))%*%W[[1]]*tau[1]
    for(m in 2:M)
      Z <- Z + (Y[[m]]-tcrossprod(Zn,Wn[[m]]))%*%W[[m]]*tau[m]
    Z  <- Z %*% covZ
    ZZ <- crossprod(C*Z, Z) + sum(C) * covZ

    list(Z=Z,covZ=covZ,ZZ=ZZ) #return
}


up_Zn <- function(Y,Wn,WWn,tau,W,Z,C) {
    ## Update latent Zn

    K  <- ncol(Wn[[1]])
    M <- length(Y)
    S  <- ncol(C)
    Zn <- matrix(0,N,K)
    ZZn<- vector("list",length=S)
    covZn <- vector("list",length=S)

    for (s in 1:S) {
        covZn[[s]] <- diag(1,K)
        for (m in 1:M) 
            covZn[[s]] <- covZn[[s]] + tau[[s]][m]*WWn[[m]]
        covZn[[s]] <- chol2inv(chol(covZn[[s]]))
        tmp <- 0
        for (m in 1:M) {
            tmp <- tmp + C[,s]*((Y[[m]]-tcrossprod(Z[[s]],W[[s]][[m]]))%*%Wn[[m]]*tau[[s]][m])
        }
        Zn <- Zn + tmp%*%covZn[[s]]
    }
    for (s in 1:S) {
        ZZn[[s]] <- crossprod(C[,s]*Zn,Zn) + sum(C[,s])*covZn[[s]]
    }

    list(Zn=Zn,covZn=covZn,ZZn=ZZn) #return
}

up_C <- function(Y,Z,covZ,Zn,covZn,W,WW,Wn,WWn,logtau,tau,
                 R,log_gamma,log_1gamma,log_pi) {
    ## update cluster assignments C
    
    D <- unlist(lapply(Y,ncol)) # data dimensionality
    M <- length(Y)
    N <- nrow(Y[[1]])  # num of samples
    S <- length(log_gamma) # num of clusters
    K <- ncol(W[[1]][[1]]) # num of GFA components
    C <- matrix(NA,N,S)
    M <- length(Y)

    for (s in 1:S) {
        #logdet_covZ <- tmp <- determinant(covZ[[s]],logarithm=TRUE)$modulus
        logdet_covZ <- 2*sum(log(diag(chol(covZ[[s]]))))
        #print(paste('logdet:',logdet_covZ-tmp))
        logdet_covZn<- 2*sum(log(diag(chol(covZn[[s]]))))
        # prior prob. + R likelihood term
        C[,s] <- log_pi[s] + ALPHA * (R * log_gamma[s] + (1-R) * log_1gamma[s])
        # log-likelihood
        C[,s] <- C[,s] + sum(D*logtau[[s]])/2

        for (tt in 1:N) {
            for (m in 1:M) {
                E  <- Y[[m]][tt,,drop=F] - tcrossprod(Zn[[s]][tt,],Wn[[m]]) 
                # E2 = < |x-Wn*zn|^2 >
                E2 <- sum((Y[[m]][tt,])^2) + sum(WWn[[m]]*(crossprod(Zn[[s]][tt,,drop=F])+covZn[[s]])) - 
                      2*Y[[m]][tt,,drop=F] %*% Wn[[m]] %*% t(Zn[[s]][tt,,drop=F]) 
                C[tt,s] <- C[tt,s] - tau[[s]][m]/2 * 
                    (E2+sum(WW[[s]][[m]]*(crossprod(Z[[s]][tt,,drop=F])+covZ[[s]]))
                     -2*E%*%W[[s]][[m]]%*%t(Z[[s]][tt,,drop=F]))
            }
            C[tt,s] <- C[tt,s] - .5*(-logdet_covZ + sum(diag(covZ[[s]])) + sum(Z[[s]][tt,]^2) - K)
            C[tt,s] <- C[tt,s] - .5*(-logdet_covZn + sum(diag(covZn[[s]])) + sum(Zn[[s]][tt,]^2) - K) 
        }
    }
    # normalize
    C <- t(apply(C,1,function(cc) exp(cc-max(cc))/sum(sort(exp(cc-max(cc))))))
    list(C=C) #return
}

up_gamma <- function(alpha_0g,beta_0g,R,C) {
    ## update gamma, the "labeling" of clusters
    a_gamma <- alpha_0g + colSums(R*C)
    b_gamma <- beta_0g + colSums((1-R)*C)
    gamma <- a_gamma / (a_gamma + b_gamma)
    log_gamma <- digamma(a_gamma) - digamma(a_gamma+b_gamma)
    log_1gamma<- digamma(b_gamma) - digamma(a_gamma+b_gamma)
    list(a_gamma=a_gamma,b_gamma=b_gamma,gamma=gamma,
         log_gamma=log_gamma,log_1gamma=log_1gamma) #return
}

up_pi <- function(C,alpha_pi,log_pi) {
    # E[ln p(pi)] - E[ln q(pi)]
    lb <- 0
    logB <- function(a) sum(sort(lgamma(a))) - lgamma(sum(sort(a))) # logarithm of the beta function
    lb.pi <- 0 # log(Unif()) # -logB(rep(1,S)) # all coefficients in rest (alpha_k-1) = 0
    lb.qi <- -(logB(alpha_pi) + (sum(alpha_pi)-S)*digamma(sum(alpha_pi)) - sum((alpha_pi-1)*digamma(alpha_pi)))
    lb <- lb + lb.pi - lb.qi

    # E[ln p(c)] - E[ln q(c)]
    lb.pc <- sum(sort(C%*%log_pi)) # = <log prod_k pi_k^{P(c_t = k}> = sum_k P(c_t = k) <log pi_k>
    #lb.qc <- sum(rowSums(C*log(C))) # would be a NaN too easily
    lb.qc <- sum(sort(apply(C,1,function(ct) sum(sort(ct*log(ct)),na.rm=T))))
    lb1 <- lb + lb.pc - lb.qc

    old_alpha <- alpha_pi
    old_log   <- log_pi

    ## update PI, cluster probabilities

    alpha_pi <- 1 + colSums(C)
    PI <- alpha_pi/sum(alpha_pi)
    log_pi <- digamma(alpha_pi) - digamma(sum(alpha_pi))
    
    # E[ln p(pi)] - E[ln q(pi)]
    lb <- 0
    logB <- function(a) sum(lgamma(a)) - lgamma(sum(a)) # logarithm of the beta function
    lb.pi <- 0 # log(Unif()) # -logB(rep(1,S)) # all coefficients in rest (alpha_k-1) = 0
    lb.qi <- -(logB(alpha_pi) + (sum(alpha_pi)-S)*digamma(sum(alpha_pi)) - sum((alpha_pi-1)*digamma(alpha_pi)))
    lb <- lb + lb.pi - lb.qi

    # E[ln p(c)] - E[ln q(c)]
    lb.pc <- sum(C%*%log_pi) # = <log prod_k pi_k^{P(c_t = k}> = sum_k P(c_t = k) <log pi_k>
    lb.pc2 <- sum(sort(C%*%log_pi))
    if (lb.pc != lb.pc2) {print("!!!!!!!")}
    #lb.qc <- sum(rowSums(C*log(C))) # would be a NaN too easily
    lb.qc <- sum(apply(C,1,function(ct) sum(ct*log(ct),na.rm=T)))
    lb2 <- lb + lb.pc - lb.qc
    #print(paste(lb2,lb1))
    if (lb2 < lb1) { 
        print(paste('up_pi:',lb2,lb1,lb2-lb1))
        #browser()
    }

    list(alpha_pi=alpha_pi, PI=PI, log_pi=log_pi) #return
}

up_W <- function(Y,Z,ZZ,Zn,Wn,alpha,tau,C) {
    ## update W, the loading matrix for a single cluster

    M <- length(tau) # num of views
    K <- ncol(Z)     # num of GFA components
    D <- unlist(lapply(Y,ncol))
    W <- WW <- covW <- vector("list",length=M)
    D <- unlist(lapply(Wn,nrow))

    for (m in 1:M) {
        # Efficient and robust way of computing
        # solve(sqrt(diag(alpha)) + tau * ZZ^T)
        tmp <- 1/sqrt(alpha[m,])
        covW[[m]] <- 1/tau[m] * outer(tmp,tmp) *
            chol2inv(chol(outer(tmp,tmp)*ZZ + diag(1/tau[m],K)))

        # remember to subtract the noise model!
        W[[m]] <- crossprod(Y[[m]]-tcrossprod(Zn,Wn[[m]]),C*Z) %*% covW[[m]] * tau[m]		

        WW[[m]] <- crossprod(W[[m]]) + covW[[m]]*D[m]
    }

    list(W=W,covW=covW,WW=WW) #return
}

up_Wn <- function(Y,Zn,ZZn,Z,W,C,alphan,tau) {
    ## update noise model loading matrix

    S <- length(W) # num of clusters
    K <- ncol(alphan)
    M <- nrow(alphan) # num of views
    D <- unlist(lapply(W[[1]],nrow))

    Wn <- lapply(1:M, function(m) matrix(0,D[m],K))
    covWn <- vector("list",length=M)
    WWn <- vector("list",length=M)
    for (m in 1:M) {
        covWn[[m]] <- diag(alphan[m,])
        for (s in 1:S) {
            # ZZn should contain latents weighted by c_t
            covWn[[m]] <- covWn[[m]] + tau[[s]][m]*ZZn[[s]]
        }
        covWn[[m]] <- chol2inv(chol(covWn[[m]]))

        for (s in 1:S) {
            Wn[[m]] <- Wn[[m]] + crossprod(Y[[m]]-tcrossprod(Z[[s]],W[[s]][[m]]),C[,s]*Zn[[s]]) * tau[[s]][m]
        }
        Wn[[m]] <- Wn[[m]] %*% covWn[[m]]
        WWn[[m]] <- crossprod(Wn[[m]]) + D[m]*covWn[[m]]
    }
    list(Wn=Wn,WWn=WWn,covWn=covWn) #return
}

up_alpha <- function(W,WW,alpha_0,beta_0) {
    ## update ARD parameters
    M <- length(WW)    #num of views
    K <- ncol(W[[1]])
    D <- unlist(lapply(W,nrow)) #data dimensions
    a_ard <- alpha_0 + D/2
    b_ard <- matrix(NA,M,K)
    alpha <- matrix(NA,M,K)
    logalpha <- matrix(NA,M,K)
    for (m in 1:M) {
        b_ard[m,]    <- beta_0 + diag(WW[[m]])/2
        alpha[m,]    <- a_ard[m] / b_ard[m,]
        logalpha[m,] <- digamma(a_ard[m]) - log(b_ard[m,])
    }
    list(a_ard=a_ard,b_ard=b_ard,alpha=alpha,logalpha=logalpha)
}

up_tau <- function(alpha_0t,beta_0t,Y,C,Zn,ZZn,Wn,WWn,Z,ZZ,W,WW) {
    ## Update tau, the noise precisions
    M <- length(Y)
    D <- unlist(lapply(Wn,nrow)) #data dimensions

    a_tau <- alpha_0t + sum(C)*D/2
    b_tau <- rep(NA,length(Y))
    for(m in 1:M) {
        E  <- Y[[m]]-tcrossprod(Zn,Wn[[m]])
        E2 <- sum(C*Y[[m]]^2) + sum(WWn[[m]]*ZZn) - 2*sum(C*Y[[m]]*tcrossprod(Zn,Wn[[m]]))
        E22 <- sum(C*E^2)
        if (E22 > E2) print("ERROR: E[X]^2 < E[X^2]")
        b_tau[m] <- beta_0t + ( E2 + sum(WW[[m]]*ZZ) - 
               2*sum(C*E*tcrossprod(Z,W[[m]])) )/2
    }
    tau <- a_tau/b_tau
    logtau <- digamma(a_tau) - log(b_tau)
    list(a_tau=a_tau,b_tau=b_tau,tau=tau,logtau=logtau) #return
}


classify_test <- function(Y, model) {
    # Need to estimate the latent variables given the global model posterior and new data

    tau  <- model$tau
    logtau<- model$logtau
    W    <- model$W
    WW   <- model$WW
    PI   <- model$PI
    log_pi<- model$log_pi
    Wn   <- model$Wn
    WWn  <- model$WWn
    N    <- nrow(Y[[1]])
    S    <- length(model$W)
    K    <- ncol(W[[1]][[1]])
    M    <- length(Y)
    D    <- unlist(lapply(Y,function(y) ncol(y)))
    gamma<- model$gamma
    log_gamma <- model$log_gamma
    log_1gamma<- model$log_1gamma
    
    # need to init these
    C    <- matrix(1/S,N,S)
    R    <- rep(0.5,N)
    Z    <- lapply(1:S,function(s) matrix(rnorm(N*K), N, K)) 
    Zn   <- lapply(1:S,function(s) matrix(rnorm(N*K), N, K)) 
    covZ <- lapply(1:S,function(s) diag(1,K))
    covZn<- lapply(1:S,function(s) diag(1,K)) 
    ZZ <- ZZn <- list()

    # find the posterior of the latent variables
    # for new data Y
    for (iter in 1:20) {
        # update cluster assignments C
        res <- up_C(Y,Z,covZ,Zn,covZn,W,WW,Wn,WWn,logtau,tau,
                    R,log_gamma,log_1gamma,log_pi)
        C   <- res$C
        # update cluster specific latents Z
        for (s in 1:S) {
            res       <- up_Z(Y,W[[s]],WW[[s]],tau[[s]],Zn[[s]],Wn,C[,s])
            Z[[s]]    <- res$Z
            ZZ[[s]]   <- res$ZZ
            covZ[[s]] <- res$covZ
        }
        # update noise latents Zn
        for (s in 1:S) { 
            res   <- up_Z(Y,Wn,WWn,tau[[s]],Z[[s]],W[[s]],C[,s])
            Zn[[s]]    <- res$Z
            ZZn[[s]]   <- res$ZZ
            covZn[[s]] <- res$covZ
        }
    }
    R <- rowSums(C %*% diag(gamma)) # estimate relevance
    list(C=C,R=R,Z=Z,Zn=Zn)
}

calc_lb <- function(Y,R,C,Z,covZ,ZZ,Zn,covZn,ZZn,logtau,logalpha,logalphan,
                    tau,a_tau,b_tau,a_ard,b_ard,a_ardn,b_ardn,
                    W,WW,covW,Wn,WWn,covWn,alpha,alphan,
                    alpha_0,beta_0,alpha_0t,beta_0t,alpha_0g,
                    beta_0g,gamma,log_gamma,log_1gamma,a_gamma,
                    b_gamma,PI,log_pi,alpha_pi) {
    # data size
    S <- length(W)
    M <- length(Y)
    N <- nrow(Y[[1]])
    K <- ncol(Wn[[1]])
    D <- unlist(lapply(Y,ncol))
    Ds <- sum(D)                # Total number of features
    const <- - N*Ds/2*log(2*pi) # Constant factors for the lower bound

    # Calculate the lower bound.
    # Consists of calculating the likelihood term and 
    # KL-divergences between the factorization and the priors
    #
    lb <- 0
    
    # update statistic
    ZZ <- ZZn <- vector("list",length=S)
    for (s in 1:S) {
        ZZ[[s]]<- crossprod(C[,s]*Z[[s]],Z[[s]]) + sum(C[,s])*covZ[[s]]
        ZZn[[s]] <- crossprod(C[,s]*Zn[[s]],Zn[[s]]) + sum(C[,s])*covZn[[s]]
    }

    # The precision terms
    for (s in 1:S) {
        temp <- rep(0,M)
        for(m in 1:M) {
          E  <- Y[[m]]-tcrossprod(Zn[[s]],Wn[[m]])
          E2 <- sum(C[,s]*Y[[m]]^2) + sum(WWn[[m]]*ZZn[[s]]) - 2*sum(C[,s]*Y[[m]]*tcrossprod(Zn[[s]],Wn[[m]]))
          temp[m] <- ( E2 + sum(WW[[s]][[m]]*ZZ[[s]]) - 
                                      2*sum(C[,s]*E*tcrossprod(Z[[s]],W[[s]][[m]])) )/2
        }
        # likelihood from this cluster
        lb.p <- const+sum(sum(C[,s])*D/2*logtau[[s]])-sum(temp*tau[[s]])
        lb <- lb + lb.p
    }

    # E[ln p(r)]
    for (s in 1:S) {
        lb.r <- ALPHA*sum(sort(R*C[,s]*log_gamma[s] + (1-R)*C[,s]*log_1gamma[s])) 
        lb <- lb + lb.r
    }

    # E[ ln p(Z)] - E[ ln q(Z) ]
    for (s in 1:S) {
        lb.px <- - sum(diag(ZZ[[s]]))/2 - sum(diag(ZZn[[s]]))/2
        lb.qx <- - sum(C[,s])*(sum(log(svd(covZ[[s]],nu=0,nv=0)$d))/2 - K/2) -
                   sum(C[,s])*(sum(log(svd(covZn[[s]],nu=0,nv=0)$d))/2 - K/2) 
        lb <- lb + lb.px - lb.qx
    }

    # E[ ln p(W)] - E[ ln q(W)]
    for (s in 1:S) {
        lb.pw <- D[1]/2*sum(logalpha[[s]][1,]) -  sum(diag(WW[[s]][[1]])*alpha[[s]][1,])/2
        for(m in 2:M) {
          lb.pw <- lb.pw +
            D[m]/2*sum(logalpha[[s]][m,]) -  sum(diag(WW[[s]][[m]])*alpha[[s]][m,])/2
        }
        lb.qw <- -D[1]*sum(log(svd(covW[[s]][[1]],nu=0,nv=0)$d))/2 -D[1]*K/2
        for(m in 2:M) {
          lb.qw <- lb.qw - D[m]*sum(log(svd(covW[[s]][[m]],nu=0,nv=0)$d))/2 -D[m]*K/2
        }
        lb <- lb + lb.pw - lb.qw
    }
    lb.pw <- D[1]/2*sum(logalphan[1,]) -  sum(diag(WWn[[1]])*alphan[1,])/2
    for(m in 2:M) {
      lb.pw <- lb.pw +
        D[m]/2*sum(logalphan[m,]) -  sum(diag(WWn[[m]])*alphan[m,])/2
    }
    lb.qw <- -D[1]*sum(log(svd(covWn[[1]],nu=0,nv=0)$d))/2 -D[1]*K/2
    for(m in 2:M) {
      lb.qw <- lb.qw - D[m]*sum(log(svd(covWn[[m]],nu=0,nv=0)$d))/2 -D[m]*K/2
    }
    lb <- lb + lb.pw - lb.qw

    # E[ln p(alpha)] - E[ln q(alpha)]
    for (s in 1:S) {
        lb.pa <- M*K*( -lgamma(alpha_0) + alpha_0*log(beta_0) ) + (alpha_0-1)*sum(logalpha[[s]]) - beta_0*sum(alpha[[s]])
        lb.qa <- -K*sum(lgamma(a_ard[[s]])) + sum(a_ard[[s]]*rowSums( log(b_ard[[s]]) )) + sum((a_ard[[s]]-1)*rowSums(logalpha[[s]])) - sum(b_ard[[s]]*alpha[[s]])
        lb <- lb + lb.pa - lb.qa
    }

    lb.pa <- M*K*( -lgamma(alpha_0) + alpha_0*log(beta_0) ) + (alpha_0-1)*sum(logalphan) - beta_0*sum(alphan)
    lb.qa <- -K*sum(lgamma(a_ardn)) + sum(a_ardn*rowSums( log(b_ardn) )) + sum((a_ardn-1)*rowSums(logalphan)) - sum(b_ardn*alphan)
    lb <- lb + lb.pa - lb.qa

    # E[ln p(tau)] - E[ln q(tau)]
    for (s in 1:S) {
        lb.pt <- -M*lgamma(alpha_0t) + M*alpha_0t*log(beta_0t) + sum((alpha_0t-1)*logtau[[s]]) - sum(beta_0t*tau[[s]])
        lb.qt <- -sum(lgamma(a_tau[[s]])) + sum(a_tau[[s]]*log(b_tau[[s]])) + sum((a_tau[[s]]-1)*logtau[[s]]) - sum(b_tau[[s]]*tau[[s]])
        lb <- lb + lb.pt - lb.qt
    }

    # E[ln p(pi)] - E[ln q(pi)]
    logB <- function(a) sum(lgamma(a)) - lgamma(sum(a)) # logarithm of the beta function
    lb.pi <- 0 # log(Unif()) # -logB(rep(1,S)) # all coefficients in rest (alpha_k-1) = 0
    lb.qi <- -(logB(alpha_pi) + (sum(alpha_pi)-S)*digamma(sum(alpha_pi)) - sum((alpha_pi-1)*digamma(alpha_pi)))
    lb <- lb + lb.pi - lb.qi

    # E[ln p(c)] - E[ln q(c)]
    lb.pc <- sum(sort(C%*%log_pi)) # = <log prod_k pi_k^{P(c_t = k}> = sum_k P(c_t = k) <log pi_k>
    #lb.qc <- sum(rowSums(C*log(C))) # would be a NaN too easily
    lb.qc <- sum(sort(apply(C,1,function(ct) sum(ct[order(ct)]*log(ct[order(ct)]),na.rm=T))))
    lb <- lb + lb.pc - lb.qc

    # E[ln p(gamma)] - E[ln q(gamma)]
    for (s in 1:S) {
        lb.pgamma <- (alpha_0g-1)*log_gamma[s] + (beta_0g-1)*log_1gamma[s] - logB(c(alpha_0g,beta_0g))
        lb.qgamma <- -(logB(c(a_gamma[s],b_gamma[s])) - (a_gamma[s]-1)*digamma(a_gamma[s]) -
                     (b_gamma[s]-1)*digamma(b_gamma[s]) + (a_gamma[s]+b_gamma[s]-2)*digamma(a_gamma[s]+b_gamma[s]))
        lb <- lb + lb.pgamma - lb.qgamma
    }
    lb
}

getDefaultOpts <- function(){
  #
  # A function for generating a default set of parameters.
  #
  # To run the algorithm with other values:
  #   opts <- getDefaultOpts()
  #   opts$opt.method <- "BFGS"
  #   model <- GFA(Y,K,opts)

  #
  # Initial value for the noise precisions. Should be large enough
  # so that the real structure is modeled with components
  # instead of the noise parameters (see Luttinen&Ilin, 2010)
  #  Values: Positive numbers, but generally should use values well
  #          above 1
  #
  init.tau <- 10^3

  #
  # Parameters for controlling when the algorithm stops.
  # It stops when the relative difference in the lower bound
  # falls below iter.crit or iter.max iterations have been performed.
  #
  iter.crit <- 10^-6
  iter.max <- 10^5
  
  #
  # Hyperparameters
  # - alpha_0, beta_0 for the ARD precisions
  # - alpha_0t, beta_0t for the residual noise predicions
  # - alpha_0g, beta_0g for the output distribution
  #
  prior.alpha_0 <- prior.beta_0 <- 1e-14
  prior.alpha_0t <- prior.beta_0t <- 1e-14
  alpha_0g <- beta_0g <- .5

  #
  # Verbosity level
  #  0: Nothing
  #  1: Final cost function value for each run of GFAexperiment()
  #  2: Cost function values for each iteration
  #
  verbose <- 2
 
  return(list(init.tau=init.tau, iter.crit=iter.crit,
              iter.max=iter.max,alpha_0g=alpha_0g,beta_0g=beta_0g,
              prior.alpha_0=prior.alpha_0,prior.beta_0=prior.beta_0,
              prior.alpha_0t=prior.alpha_0t,prior.beta_0t=prior.beta_0t,
              verbose=verbose))
}
