subhmm <- function(logR, logOR, purity=0.7, ploidy=1.5, clonal.prop=0.5, 
                    logR.var=0.5, logOR.var=0.5, df=3,
                    genoStates=c("0", "A", "AA", "AB","AAA","AAB","AAAA","AABB",
                                 "AAAB","AAAAA","AAABB","AAAAB"), prob0=NULL, 
                    mainclone.trans=NULL, subclone.trans=NULL,
                    subclone.prob=c(0.5, 0.5), subclone.probMat=NULL,
                    maxiter=1000, logOR2.min=1e-6, logOR.var.min=1e-5,
                    df.min=1e-5, df.max=100,
                    loglike.eps=0.01, parm.eps=1e-7, nLociRegion=50, print=1) {
  if (length(logR) != length(logOR)) stop("ERROR: logR and logOR must have the same length")
  if ((purity >= 1) || (purity <= 0)) stop("ERROR: purity must be betwen 0 and 1")
  if (ploidy <= 0) stop("ERROR: ploidy must be positive")
  if ((clonal.prop <= 0) || (clonal.prop >= 1)) stop("ERROR: clonal.prop must be positive")
  if (logR.var <= 0) stop("ERROR: logR.var must be positive")
  if (logOR.var <= 0) stop("ERROR: logOR.var must be positive")
  if (df <= 0) stop("ERROR: df must be positive")

  main_zk                     <- genoStates
  tmp                         <- callGenoStates(main_zk)
  main_zk                     <- tmp$zk
  alleles                     <- tmp$unq_alleles
  gs.null                     <- tmp$gs.null
  names(main_zk)              <- NULL
  mainJ                       <- length(main_zk)
  if (mainJ < 2) stop("ERROR: too few genotype states")
  tauJ                        <- mainJ - 1
  allzk_main                  <- c(main_zk, rep(main_zk,each=tauJ))
  sub_zk                      <- c(sapply(1:mainJ, function(x) main_zk[-x]))
  mz                          <- str_count(allzk_main, alleles[1]) # A allele frequency corresponding to z
  mz_sub                      <- str_count(sub_zk, alleles[1])     # A allele of subclone genotype
  pz                          <- str_count(allzk_main, alleles[2]) # B allele frequency corresponding to z
  pz_sub                      <- str_count(sub_zk, alleles[2])     # B allele of subclone genotype
  ctzs                        <- nchar(allzk_main) # main genotype copy number
  ctzs[which(allzk_main==gs.null)] <- 0
  ctms                         <- nchar(sub_zk)  # subclone genotype copy number
  ctms[which(sub_zk==gs.null)] <- 0
  J                            <- length(allzk_main)
  lr                           <- logR
  logor2                       <- logOR^2
  pr0                          <- purity
  lpr0                         <- log((1-pr0)/pr0)
  sigma20_lr                   <- logR.var
  lsigma20_lr                  <- log(sigma20_lr) 
  lpsi0                        <- log(ploidy, base=2)
  lsigma20_lor                 <- log(logOR.var) 
  r0                           <- prob0
  lsz0                         <- log((1-clonal.prop)/clonal.prop)
  P5_LB                        <- log(logOR.var.min) 
  v0                           <- df
  V0_LB                        <- df.min 
  if (V0_LB <= 0) stop("ERROR: with df.min")
  if (v0 <= V0_LB) stop("ERROR: with df or df.min")
  V0_UB                       <- df.max 
  if (V0_UB < V0_LB) stop("ERROR: with df.max")
  if (v0 >= V0_UB) stop("ERROR: with df or df.max")
  if (!is.finite(gamma(df.min/2))) warning("df.min may be too small")
  if (!is.finite(gamma(df.max/2))) warning("df.max may be too large")
  if (!is.finite(lpr0)) stop("ERROR with purity")
  if (!is.finite(lpsi0)) stop("ERROR with ploidy")
  if (!is.finite(lsigma20_lr)) stop("ERROR with logR.var")
  if (!is.finite(lsigma20_lor)) stop("ERROR with logOR.var")
  if (!is.finite(lsz0)) stop("ERROR with clonal.prop")
  if (!is.finite(P5_LB)) stop("ERROR: logOR.var.min is too small")

  if (is.null(r0)) r0 <- rep(1/mainJ, mainJ)
  if (length(r0) != mainJ) stop("ERROR with prob0")
  if (maxiter < 0) stop("ERROR: with maxiter")

  # Prevent log(0)
  LOG0ARG <- 1e-300
  LOG0    <- log(LOG0ARG)
  DEBUG   <- print > 9
  
  main_geneidx <- c(1:mainJ,rep(1:mainJ, each=tauJ))
  u0           <- subclone.prob
  hatp0_M      <- mainclone.trans
  if (is.null(hatp0_M)) hatp0_M <- matrix(c(rep(c(1-tauJ/5000,rep(1/5000, J/mainJ)),tauJ),1-tauJ/5000), J/mainJ, J/mainJ)
  if (nrow(hatp0_M) != mainJ)  stop("ERROR: with mainclone.trans")
  if (ncol(hatp0_M) != mainJ)  stop("ERROR: with mainclone.trans")

  hatp0_m      <- as.numeric(subclone.trans)
  if (!length(hatp0_m)) hatp0_m <- c(0.999, 0.001, 0.001, 0.999)
  if (length(hatp0_m) != 4) stop("ERROR: with subclone.trans")

  hatp0_mgtype <- subclone.probMat
  if (is.null(hatp0_mgtype)) {
    hatp0_mgtype <- matrix(1/tauJ, mainJ, tauJ)

    #hatp0_mgtype <- matrix((1-0.99)/(tauJ-1),mainJ, tauJ)
    #hatp0_mgtype[1,1] <- 0.99
    #if (mainJ > 1) hatp0_mgtype[2,3] <- 0.99
    #if (mainJ > 2) hatp0_mgtype[3,2] <- 0.99
    #if (mainJ > 3) hatp0_mgtype[4,3] <- 0.99
    #if (mainJ > 4) hatp0_mgtype[5:mainJ,] <- 1/tauJ
  }

  hatp0i_mgtype <- function(i, hh) {
    vec <- hatp0_M[i,]*hatp0_m[hh]
    ret <- matrix(rep(vec, each=tauJ), nrow=mainJ, ncol=tauJ, byrow=TRUE)*hatp0_mgtype
    ret <- as.vector(t(ret))
    ret
  }

  # hatp0
  hatp0 <- matrix(0,J,J)
  #hatp0 <- matrix(c(rep(c(1-(J-1)/5000,rep(1/5000,J)),(J-1)),1-(J-1)/5000),J,J) 
  for(i in 1:mainJ) hatp0[i,] <- c(hatp0_M[i,]*hatp0_m[1], hatp0i_mgtype(i,2))
  hatp0low <- NULL
  for(i in 1:mainJ){
    temp0    <- matrix(rep(c(hatp0_M[i,]*hatp0_m[3],hatp0i_mgtype(i,4)),each=tauJ),tauJ,J)
    hatp0low <- rbind(hatp0low, temp0)
  }
  hatp0[(mainJ+1):J,1:J] <- hatp0low

  # Remove missing. logOR can have missing values
  tmp             <- is.finite(lr)
  tmp[is.na(tmp)] <- FALSE
  if (!all(tmp)) {
    lr     <- lr[tmp]
    logor2 <- logor2[tmp] 
  }  

  # Replace small values of logor2 with minLogOR2
  tmp <- logor2 < logOR2.min
  tmp[is.na(tmp)] <- FALSE
  if (any(tmp)) logor2[tmp] <- logOR2.min
  n <- length(lr)

  notMissing  <- is.finite(logor2)
  noMissing   <- all(notMissing)
  if (!noMissing) logor2[!notMissing] <- -999999999.0

  

  ###################
  # Local functions #
  ###################
  logL_ <- function(parm){

    #if (DEBUG) print("BEGIN: logL_")
    p6 <- parm[6]
    p1 <- parm[1]
    p2 <- parm[2]
    p3 <- parm[3]
    p4 <- parm[4]
    p5 <- parm[5]

    # Reparameterize v0
    p62 <- p6*p6
    p6  <- V0_LB + (V0_UB - V0_LB)*p62/(1 + p62) 

    mu0_lrs    <- mu0_lr(J, p1, p2, p3)
    mu0_lors   <- mu0_lor(J, p2, p3)
    sigma20_lr <- exp(p4)
    ep5        <- exp(-p5)
    ncp        <- mu0_lors^2*ep5

    # Problem when ep5 is large (Inf), for now set dmat to 0.
    # Perhaps set bounds for parm[5]
    if (is.finite(ep5)) {    
      #dmat       <- dchisq(rep(logor2*ep5, each=J), df=1, ncp=rep(ncp, times=n), log=T)
      vec1             <- rep(NA, n)
      vec1[notMissing] <- logor2[notMissing]
      vec1             <- rep(vec1*ep5, each=J)
      vec2             <- rep(ncp, times=n)
      dmat             <- dchisq(vec1, df=1, ncp=vec2, log=T)
      dmat             <- matrix(dmat, nrow=n, ncol=J, byrow=TRUE)
    } else {
      stop("ERROR: set bounds on parm[5]")
      dmat             <- matrix(0, nrow=n, ncol=J, byrow=TRUE)
    }
    mat  <- (-0.5*log(2*pi*sigma20_lr) + lgamma((p6+1)/2) - lgamma(p6/2) + 0.5*p6*log(p6/2) -
             0.5*(p6+1)*log(0.5*(p6+(lr-matrix(mu0_lrs, nrow=n, ncol=J, byrow=TRUE))^2/sigma20_lr)) )*pzks1
  
    if (noMissing) {
      mat <- mat - pzks1*(p5 - dmat)
    } else {
      mat[notMissing, ] <- mat[notMissing, , drop=FALSE] - pzks1[notMissing, , drop=FALSE]*(p5 - dmat[notMissing, , drop=FALSE])
    }
    logf <- sum(mat)

    #if (DEBUG) print("END: logL_")

    return(-logf)   

  } # END: logL_

  

  akh_all.ca_all <- function() {
    if (DEBUG) print("BEGIN: akh_all.ca_all")

    c1              <- 1/sum(a1d) 
    a1h             <- c1*a1d
    #akh_all.new     <- matrix(NA, n,J)
    #akh_all.new[1,] <- a1h
    #ca_all.new      <- numeric(n)
    #ca_all.new[1]   <- c1
    a1h             <- matrix(a1h, nrow=1)

    pt5v0                  <- 0.5*v0
    pt5v0POWpt5v0          <- pt5v0^pt5v0
    v0Plus1Over2           <- (v0+1)/2
    gamma_v0Over2          <- gamma(v0/2)

    tmp1  <- pt5v0POWpt5v0*gamma((v0+1)/2)*(2*pi*sigma20_lr)^-0.5
    tmp2  <- 2*sigma20_lr
    tmp3  <- exp(-lsigma20_lor)
    tmp4  <- mu0_lors^2*tmp3
    ret_code    <- as.integer(-1)
    ret_akh_all <- as.numeric(rep(-9999, n*J))
    ret_ca_all  <- as.numeric(rep(-9999, n))
    tmp <- .C("akh_ca", as.integer(n), as.integer(J), as.integer(mainJ), as.integer(tauJ), 
                   as.numeric(pt5v0), as.numeric(pt5v0POWpt5v0), as.numeric(v0Plus1Over2),
                   as.numeric(gamma_v0Over2), as.numeric(tmp1), as.numeric(tmp2), 
                   as.numeric(tmp3), as.numeric(tmp4), as.numeric(logor2),
                   as.numeric(lr), as.integer(notMissing), as.numeric(mu0_lrs), 
                   as.numeric(a1h), as.numeric(hatp0), as.numeric(c1), 
                   ret_code=ret_code, ret_akh_all=ret_akh_all, ret_ca_all=ret_ca_all,
                   PACKAGE="subHMM");
    if (tmp$ret_code) stop("ERROR in akh_ca")
    ret_ca_all  <- tmp$ret_ca_all
    ret_akh_all <- matrix(tmp$ret_akh_all, byrow=TRUE, nrow=n, ncol=J)
    if (DEBUG) print("END: akh_all.ca_all")

    return(list(akh_all=ret_akh_all, ca_all=ret_ca_all))
  
  } # END: akh_all.ca_all

  bkh_all.cb_all <- function() {
    if (DEBUG) print("BEGIN: bkh_all.cb_all")
    tmp0 <- 0.5*v0
    tmp1 <- (tmp0)^(tmp0)*gamma((v0+1)/2)*(2*pi*sigma20_lr)^-0.5
    tmp2 <- gamma(v0/2)
    tmp3 <- 2*sigma20_lr
    tmp4 <- (v0+1)/2
    tmp5 <- exp(-lsigma20_lor)
    tmp6 <- mu0_lors^2*tmp5
  
    ret_code    <- as.integer(-1)
    bkh_all.new <- as.numeric(rep(-9999, n*J))
    cb_all.new  <- as.numeric(rep(-9999, n))

    tmp <- .C("bkh_cb", as.integer(n), as.integer(J), as.numeric(tmp0), as.numeric(tmp1), 
              as.numeric(tmp2), as.numeric(tmp3), as.numeric(tmp4), as.numeric(tmp5),
              as.numeric(tmp6), as.numeric(logor2), as.numeric(lr), as.integer(notMissing),
              as.numeric(mu0_lrs), as.numeric(hatp0),  
              ret_code=ret_code, ret_bkh_all=bkh_all.new, ret_cb_all=cb_all.new,
              PACKAGE="subHMM")
    if (tmp$ret_code) stop("ERROR: bad return code in bkh_cb")
    if (DEBUG) print("END: bkh_all.cb_all")
    list(cb_all=tmp$ret_cb_all, bkh_all=matrix(tmp$ret_bkh_all, nrow=n, ncol=J, byrow=TRUE))

  } # END: bkh_all.cb_all



  locFunc_hatp0 <- function() {
    if (DEBUG) print("BEGIN: locFunc_hatp0")  
    t1        <- 2*sigma20_lr 
    t2        <- exp(-lsigma20_lor)
    t3        <- sum(akh_all[n,]) 
    tveca     <- reverseCumSum(lca_all)
    tvecb     <- reverseCumSum(lcb_all)
    tvec      <- exp(tveca[3:n]-tvecb[2:(n-1)])*ca_all[2:(n-1)]/t3
    ncpVecJ   <- mu0_lors^2*t2
    tmp0      <- 0.5*v0
    tmp1      <- tmp0^tmp0*gamma((v0+1)/2)*(2*pi*sigma20_lr)^-0.5
    tmp2      <- gamma(v0/2)
    tmp3      <- (v0+1)/2
    
    if (notMissing[n]) {
      vec <- ca_all[n]*tmp1/(tmp2*(tmp0+((lr[n]-mu0_lrs)^2/t1))^tmp3)*t2*dchisq(logor2[n]*t2,df=1, ncp=ncpVecJ)/sum(akh_all[n,])
    } else {
      vec <- ca_all[n]*tmp1/(tmp2*(tmp0+((lr[n]-mu0_lrs)^2/t1))^tmp3)/sum(akh_all[n,])
    }   

    ret_code         <- -1
    ret_hatp0_M      <- rep(0, mainJ*mainJ)
    ret_hatp0_m      <- rep(0, 4)
    ret_hatp0_mgtype <- rep(0, mainJ*tauJ)
    ret_r0s          <- rep(0, mainJ*tauJ + mainJ)
    rm(tveca, tvecb)
    gc()
    if (DEBUG) print("BEGIN: C code")
    ret <- .C("locFunc_hatp0", as.integer(n), as.integer(J), as.integer(mainJ), 
                 as.integer(tauJ), as.numeric(t1), as.numeric(t2), 
                 as.numeric(tvec), as.numeric(ncpVecJ),
                 as.numeric(tmp0), as.numeric(tmp1), as.numeric(tmp2), as.numeric(tmp3), 
                 as.numeric(logor2), as.numeric(lr), as.integer(notMissing),
                 as.numeric(akh_all), as.numeric(bkh_all), as.numeric(mu0_lrs), 
                 as.numeric(hatp0), as.numeric(vec), 
                 as.numeric(hatp0_mgtype), as.numeric(pzks1),
                 ret_code=as.integer(ret_code), ret_hatp0_M=as.numeric(ret_hatp0_M), 
                 ret_hatp0_m=as.numeric(ret_hatp0_m), ret_r0s=as.numeric(ret_r0s), 
                 ret_hatp0_mgtype=as.numeric(ret_hatp0_mgtype),
                 PACKAGE="subHMM")
    if (DEBUG) print("END: C code") 
    if (ret$ret_code) stop("ERROR in locFunc_hatp0")
    ret_hatp0_M      <- matrix(ret$ret_hatp0_M, byrow=FALSE, nrow=mainJ, ncol=mainJ)
    ret_hatp0_mgtype <- matrix(ret$ret_hatp0_mgtype, byrow=TRUE, nrow=mainJ, ncol=tauJ)

    ret <- list(hatp0_M=ret_hatp0_M, hatp0_m=ret$ret_hatp0_m, r0s=ret$ret_r0s, 
                  hatp0_mgtype=ret_hatp0_mgtype)
    if (DEBUG) print("END: locFunc_hatp0")
    return(ret)

  } # END: locFunc_hatp0

  get_a1d <- function() {
    if (DEBUG) print("BEGIN: get_a1d")
    pt5v0                  <- 0.5*v0
    pt5v0POWpt5v0          <- pt5v0^pt5v0
    v0Plus1Over2           <- (v0+1)/2
    gamma_v0Over2          <- gamma(v0/2)
    gamma_v0Plus1Over2     <- gamma(v0Plus1Over2)
    twoPi                  <- 2*pi
    twoPiSimgma20lrPOWmpt5 <- (twoPi*sigma20_lr)^-0.5

    tmp1    <- pt5v0POWpt5v0*gamma_v0Plus1Over2*twoPiSimgma20lrPOWmpt5
    tmp2    <- gamma_v0Over2*(pt5v0+((lr[1]-mu0_lrs)^2/(2*sigma20_lr)))^v0Plus1Over2
    a1d     <- (tmp1/tmp2)*r0s

    if (notMissing[1]) {
      tmp3     <- dchisq(logor2[1]*exp(-lsigma20_lor),df=1, ncp=(mu0_lors)^2*exp(-lsigma20_lor))
      tmp4     <- exp(-lsigma20_lor)*tmp3
      a1d      <- a1d*tmp4
    }
    if (DEBUG) print("END: get_a1d") 
    a1d

  } # END: get_a1d

  get_cb_bkh <- function() {
    if (DEBUG) print("BEGIN: get_cb_bkh")
    # Scaling HMM with forward-backward alg. in E-step
    tmp0 <- 0.5*v0
    tmp1 <- (tmp0)^(tmp0)*gamma((v0+1)/2)*(2*pi*sigma20_lr)^-0.5
    tmp2 <- gamma(v0/2)
    tmp3 <- 2*sigma20_lr
    tmp4 <- (v0+1)/2
    tmp5 <- exp(-lsigma20_lor)
    tmp6 <- mu0_lors^2*tmp5

    ret_code    <- as.integer(-1)
    bkh_all.new <- as.numeric(rep(-9999, n*J))
    cb_all.new  <- as.numeric(rep(-9999, n))

    tmp <- .C("bkh_cb", as.integer(n), as.integer(J), as.numeric(tmp0), as.numeric(tmp1), 
              as.numeric(tmp2), as.numeric(tmp3), as.numeric(tmp4), as.numeric(tmp5),
              as.numeric(tmp6), as.numeric(logor2), as.numeric(lr), as.integer(notMissing),
              as.numeric(mu0_lrs), as.numeric(hatp0),  
              ret_code=ret_code, ret_bkh_all=bkh_all.new, ret_cb_all=cb_all.new,
              PACKAGE="subHMM")
    if (tmp$ret_code) stop("ERROR in bkh_cb")
    if (DEBUG) print("END: get_cb_bkh")
    list(cb_all=tmp$ret_cb_all, bkh_all=matrix(tmp$ret_bkh_all, nrow=n, ncol=J, byrow=TRUE))
    
  } # END: get_cb_bkh 

  logL_reg <- function(parm){
    #if (DEBUG) print("BEGIN: logL_reg")
    mu0_lrs    <- mu0_lr(J,lpsi0, lpr0, parm)
    mu0_lors   <- mu0_lor(J, lpr0, parm)
    sigma20_lr <- exp(lsigma20_lr)
  


    logf <- 0
    tmp1 <- -0.5*log(2*pi*sigma20_lr)+lgamma((v0+1)/2)-lgamma(v0/2)+0.5*v0*log(v0/2)
    mat  <- matrix(lr, nrow=n, ncol=length(mu0_lrs), byrow=FALSE) - matrix(mu0_lrs, nrow=n, ncol=length(mu0_lrs), byrow=TRUE)
    mat  <- 0.5*(v0+1)*log(0.5*(v0+mat^2/sigma20_lr))
    mat  <- (tmp1 - mat)*pzks1
    if (!noMissing) {
      tmp2 <- exp(-lsigma20_lor)
      m    <- length(mu0_lors)
      vec1 <- logor2*tmp2
      vec1[!notMissing] <- NA
      mat2 <- dchisq(matrix(vec1, nrow=n, ncol=m, byrow=FALSE), df=1, 
                     ncp=matrix(mu0_lors^2*tmp2, nrow=n, ncol=m, byrow=TRUE), log=T)
      mat[notMissing, ] <- mat[notMissing, , drop=FALSE] - pzks1[notMissing,, drop=FALSE]*lsigma20_lor + 
                           pzks1[notMissing, , drop=FALSE]*mat2[notMissing,, drop=FALSE]
    }
    logf <- sum(mat[sub_reg, ])
    #if (DEBUG) print("END: logL_reg")
    return(-logf)   
  } # END: logL_reg

  
  # These functions are needed to evaluate logL_ function
  mu0_lr <- function(nj, lpsi0, lpr0, lsz0) {
    #if (DEBUG) print("BEGIN: mu0_lr")
    ret      <- rep(NA, nj)
    pr0      <- 1/(1+exp(lpr0))
    sz0      <- 1/(1+exp(lsz0))
    jv1      <- 1:mainJ
    jv2      <- (mainJ+1):nj
    ret[jv1] <- log((2+pr0*(ctzs[jv1]-2)), base=2)-lpsi0
    tmp1     <- ctms[jv2 - mainJ]
    ret[jv2] <- log(2+(tmp1-2)*pr0+(ctzs[jv2]-tmp1)*sz0*pr0, base=2)-lpsi0
    #if (DEBUG) print("END: mu0_lr")
    ret
  
  } # END: mu0_lr

  mu0_lor <- function(nj, lpr0, lsz0) {
    #if (DEBUG) print("BEGIN: mu0_lor")
    ret      <- rep(NA, nj)
    pr0      <- 1/(1+exp(lpr0))
    sz0      <- 1/(1+exp(lsz0))
    jv1      <- 1:mainJ
    jv2      <- (mainJ+1):nj
    arg1     <- mz[jv1]*pr0+(1-pr0)
    arg2     <- pz[jv1]*pr0+(1-pr0)
    tmp      <- arg1 < LOG0ARG
    if (any(tmp)) arg1[tmp] <- LOG0ARG
    tmp      <- arg2 < LOG0ARG
    if (any(tmp)) arg2[tmp] <- LOG0ARG
    ret[jv1] <- log(arg1) - log(arg2)

    tmpj     <- jv2 - mainJ
    tmp1     <- mz_sub[tmpj]
    tmp2     <- pz_sub[tmpj]
    arg1     <- 1+(mz[jv2]-tmp1)*sz0*pr0+(tmp1-1)*pr0   
    arg2     <- 1+(pz[jv2]-tmp2)*sz0*pr0+(tmp2-1)*pr0
    tmp      <- arg1 < LOG0ARG
    if (any(tmp)) arg1[tmp] <- LOG0ARG
    tmp      <- arg2 < LOG0ARG
    if (any(tmp)) arg2[tmp] <- LOG0ARG
    ret[jv2] <- log(arg1) - log(arg2)
    #if (DEBUG) print("END: mu0_lor") 
    ret
  
  } # END: mu0_lor

  get_sub_prob <- function() {
    if (DEBUG) print("BEGIN: get_sub_prob")
    nn    <- length(ok_subzoneh_spl)
    if (!nn) {
      if (DEBUG) print("END: get_sub_prob")
      return(NULL)
    }
    ret   <- matrix(0, nn, mainJ+1) 
    colnames(ret) <- c("SubRegion", main_zk)
    wlist <- list()
    nzk   <- length(main_zk)
    for (i in 1:nzk) wlist[[i]] <- which(sub_zk %in% main_zk[i])+mainJ
    jvec  <- (mainJ+1):J
    for(ss in 1:length(ok_subzoneh_spl)){
      sub_reg <- ok_subzoneh_spl[[ss]]
      mat     <- pzks1[sub_reg, , drop=FALSE]
      denom   <- rowSums(mat[, jvec, drop=FALSE])
      vec     <- rep(NA, nzk)
      for (i in 1:nzk) {
        cols   <- wlist[[i]]
        vec[i] <- mean(rowSums(mat[, cols, drop=FALSE])/denom)
      }
      ret[ss,] <- c(sz_ests[ss,2], vec)
    }
    if (DEBUG) print("END: get_sub_prob")
    ret

  } # END: get_sub_prob

  hes_logL1 <- function(parm){
    #if (DEBUG) print("BEGIN: hes_logL1")
    mu0_lrs     <- mu0_lr(J, parm[1], parm[2], lsz0)
    mu0_lors    <- mu0_lor(J, parm[2], lsz0)
    ret_code    <- as.integer(-1)
    ret_loglike <- as.numeric(0)

  

    tmp <- .C("C_hes_logL1", as.integer(n), as.integer(J), as.numeric(parm), 
              as.numeric(mu0_lrs), as.numeric(mu0_lors), as.numeric(pi), 
              as.numeric(logor2), as.numeric(lr), as.integer(notMissing), 
              as.numeric(r0s), as.numeric(hatp0),  
              ret_code=ret_code, ret_loglike=ret_loglike, PACKAGE="subHMM")
    if (tmp$ret_code) stop("ERROR in C_hes_logL1")
    #if (DEBUG) print("END: hes_logL1")

    return(tmp$ret_loglike)

  } # END: hes_logL1

  myoptimR <- function(parm) {
    if (DEBUG) print("BEGIN: optim")

    # From reparameterizing v0 in likelihood
    p6      <- parm[6]
    parm[6] <- sqrt((p6 - V0_LB)/(V0_UB - p6)) 

    # From reparameterizing logOR.var in likelihood
    parm[5] <- sqrt(parm[5] - P5_LB) 

    ret_code <- as.integer(-1)
    nparm    <- as.integer(length(parm))
    reltol   <- sqrt(.Machine$double.eps)
    tmp <- .C("myoptimC", as.numeric(parm), nparm, as.numeric(P5_LB),
              as.numeric(V0_LB), as.numeric(V0_UB), 
              as.integer(n), as.integer(J), as.integer(mainJ), as.numeric(lr), as.numeric(logor2),
              as.integer(notMissing), as.numeric(pzks1), as.numeric(ctzs), as.numeric(ctms), 
              as.numeric(mz), as.numeric(pz), as.numeric(mz_sub), as.numeric(pz_sub), 
              as.numeric(reltol), ret_code=ret_code, PACKAGE="subHMM")

    if (tmp$ret_code) stop("ERROR in myoptimC")
    ret    <- tmp[[1]]
    if (!all(is.finite(ret))) stop("ERROR with optim: non-finite estimated parameters")
    p6     <- ret[6]
    p62    <- p6*p6
    ret[6] <- V0_LB + (V0_UB - V0_LB)*(p62/(1 + p62)) # From rep
    p5     <- ret[5]
    ret[5] <- P5_LB + p5*p5 
    if (DEBUG) print("END: optim")

    return(ret)

  } # END: myoptimR


  # Compute r0s
  tmp       <- matrix(r0*u0[2], nrow=mainJ, ncol=tauJ, byrow=FALSE)
  tmp       <- tmp*hatp0_mgtype
  r0s       <- c(r0*u0[1], as.vector(t(tmp)))
  mu0_lrs   <- mu0_lr(J,lpsi0, lpr0, lsz0)
  mu0_lors  <- mu0_lor(J, lpr0, lsz0)
  a1d       <- get_a1d()
  tmp       <- akh_all.ca_all()
  ca_all    <- tmp$ca_all
  akh_all   <- tmp$akh_all

  logL_Y    <- rep(NA, maxiter+1)
  logL0     <- sum(-log(ca_all[1:n])) + log(sum(akh_all[n,]))
  logL_Y[1] <- logL0
  if (print) {
    cat("\nBegin stage 1\n")
    cat(paste("Initial loglikelihood = ", formatValue(logL0), "\n", sep=""))
    if (print > 1) {
      PARM0 <- formatValue(c(ploidy,purity,clonal.prop,logR.var,logOR.var, df))
      names(PARM0) <- getParmNames(0, all=1)
      cat("Initial estimates:\n")
      print(PARM0)
    }
  }

  converged  <- FALSE
  for (mm in 1:maxiter) {
    gc()    

    # Scaling HMM with forward-backward alg. in E-step
    tmp     <- get_cb_bkh()
    cb_all  <- tmp$cb_all
    bkh_all <- tmp$bkh_all

    pzks             <- akh_all*bkh_all/sum(akh_all[n,])
    lca_all          <- log(ca_all)
    lcb_all          <- log(cb_all)
    pzks1            <- matrix(data=NA, nrow=n, ncol=J)
    pzks1[n ,]       <- pzks[n,]/cb_all[n]
    tmp1             <- reverseCumSum(lca_all)
    tmp2             <- reverseCumSum(lcb_all)
    pzks1[1:(n-1), ] <- exp(tmp1[2:n] - tmp2[1:(n-1)])*pzks[1:(n-1), ]
  
    tmp              <- locFunc_hatp0() # uses lca_all
    rm(lca_all, lcb_all, pzks)
    gc()
    hatp0_M          <- tmp$hatp0_M
    hatp0_m          <- tmp$hatp0_m
    r0s              <- tmp$r0s
    hatp0_mgtype     <- tmp$hatp0_mgtype

    hatp0            <- matrix(0,J,J)
    tmp1             <- hatp0_M*hatp0_m[1]
    for(i in 1:mainJ) hatp0[i,] <- c(tmp1[i,], hatp0i_mgtype(i,2))

     hatp0low    <- matrix(NA, nrow=mainJ*tauJ, ncol=mainJ*mainJ)
     tmp1        <- hatp0_M*hatp0_m[3]
     tmp3        <- 0
     for(i in 1:mainJ){
       temp0 <- matrix(rep(c(tmp1[i,],hatp0i_mgtype(i,4)),each=tauJ),tauJ,J)
       tmp2  <- tmp3 + 1
       tmp3  <- tmp2 + nrow(temp0) - 1
       hatp0low[tmp2:tmp3, ] <- temp0
     }
     hatp0[(mainJ+1):J,1:J] <- hatp0low


     PARM0 <- c(lpsi0,lpr0,lsz0,lsigma20_lr,lsigma20_lor, v0)
     rm(tmp, temp0, tmp1, tmp2, tmp3)
     gc()
  
     PARM1        <- myoptimR(PARM0)
     lpsi0        <- PARM1[1]
     lpr0         <- PARM1[2]
     lsz0         <- PARM1[3]
     lsigma20_lr  <- PARM1[4]
     lsigma20_lor <- PARM1[5]
     v0           <- PARM1[6]
     sigma20_lr   <- exp(lsigma20_lr)
     pr0          <- 1/(1+exp(lpr0))
     mu0_lrs      <- mu0_lr(J,lpsi0, lpr0, lsz0)
     mu0_lors     <- mu0_lor(J, lpr0, lsz0)
  
     if (!notMissing[1]) {
       a1d.new <- (0.5*v0)^(0.5*v0)*gamma((v0+1)/2)*(2*pi*sigma20_lr)^-0.5/(gamma(v0/2)*(0.5*v0+((lr[1]-mu0_lrs)^2/(2*sigma20_lr)))^((v0+1)/2))*r0s 
     } else {
       a1d.new <- (0.5*v0)^(0.5*v0)*gamma((v0+1)/2)*(2*pi*sigma20_lr)^-0.5/(gamma(v0/2)*(0.5*v0+((lr[1]-mu0_lrs)^2/(2*sigma20_lr)))^((v0+1)/2))*exp(-lsigma20_lor)*dchisq(logor2[1]*exp(-lsigma20_lor),df=1, ncp=(mu0_lors)^2*exp(-lsigma20_lor))*r0s                  
     }

     c1           <- 1/sum(a1d)
     a1h          <- c1*a1d
     akh_all      <- matrix(NA, n,J)
     akh_all[1,]  <- a1h
     ca_all       <- numeric(n)
     ca_all[1]    <- c1
     tmp          <- akh_all.ca_all()  
     ca_all       <- tmp$ca_all
     akh_all      <- tmp$akh_all

     logL1        <- sum(-log(ca_all)) + log(sum(akh_all[n, ]))
     logL_Y[mm+1] <- logL1
     converged    <- checkStop(logL0, logL1, PARM0, PARM1, loglike.eps, parm.eps, mm, print) 
     if (converged) break
     logL0        <- logL1

  } # END: for (mm in 1:M)

  if (!converged) {
    warning("Maximum number of iterations reached without convergence")
  } else if (print) {
    cat("Algorithm converged in ", mm, " iterations\n", sep="")
    logL_Y <- logL_Y[1:(mm+1)] 
    cat("Begin stage 2\n")
  }
  
  # Obtain estimates from the first step - global estimates
  pr0         <- 1/(1+exp(PARM1[2]))
  psi0        <- 2^PARM1[1]
  sz0         <- 1/(1+exp(PARM1[3]))
  sigma20_lr  <- exp(PARM1[4])
  sigma20_lor <- exp(PARM1[5])
  v0          <- PARM1[6] 

  idx_hgtype        <- max.col(pzks1, ties.method="first")
  subc_hidx         <- ifelse(idx_hgtype >= (mainJ+1), 1, 0)
  main_hgtype       <- numeric(n)
  tmp1              <- subc_hidx == 0
  if (any(tmp1)) main_hgtype[tmp1] <- max.col(pzks1[tmp1, 1:mainJ], ties.method="first")
  tmp1              <- max.col(pzks1[,(mainJ+1):J], ties.method="first")
  tmp2              <- main_geneidx[-c(1:mainJ)]
  tmp3              <- subc_hidx == 1
  if (any(tmp3)) main_hgtype[tmp3] <- tmp2[tmp1[tmp3]]

  ###########
  # Stage 2 #
  ###########
  # identify subclone regions based on the model
  subzoneh_spl    <- NULL
  ok_subzh_idc    <- NULL
  ok_subzoneh_spl <- NULL
  sz_ests         <- NULL
  subzoneh        <- which(subc_hidx==1)

  if (length(subzoneh)) {
    subzoneh_spl    <- split(subzoneh, cumsum(c(1,diff(subzoneh)!=1 | diff(main_hgtype[subzoneh])!=0)))
    ok_subzh_idc    <- sapply(1:length(subzoneh_spl), function(x) ifelse(length(subzoneh_spl[[x]])>=nLociRegion,1,0))
    ok_subzoneh_spl <- subzoneh_spl[which(ok_subzh_idc==1)]
    names(ok_subzoneh_spl) <- NULL
  }

  # region-specific Sz0
  lsz_est <- function(sub_reg) {
    if (DEBUG) print("BEGIN: lsz_est")
    m_est_s <- optim(lsz0, logL_reg, method="BFGS", control=list(trace=0))
    lsz0    <- m_est_s$par
    if (DEBUG) print("END: lsz_est")
    return(lsz0)
  } # END: lsz_est

  if (length(ok_subzoneh_spl)) {
    sz_ests <- matrix(0,nrow=length(ok_subzoneh_spl), ncol=2)
    colnames(sz_ests) <- c("Proportion", "SubRegion")
    for(ss in 1:length(ok_subzoneh_spl)){
      sub_reg      <- ok_subzoneh_spl[[ss]]
      lsz1         <- lsz_est(sub_reg)
      sz_ests[ss,] <- c(1/(1+exp(lsz1)),ss)
    }
  }

  hat_logr  <- mu0_lrs[main_hgtype]
  hat_logor <- mu0_lors[main_hgtype]
  

  ##################
  # BEGIN sub_prob #
  ##################
  sub_prob       <- get_sub_prob()
  sub_hgtype_dec <- NULL
  sub_hgtype_fin <- NULL 

  # final copy number profile considering subclone 
  if (!is.null(sub_prob)) {
    if (DEBUG) print("BEGIN: final copy number")
    sub_hgtype_dec <- apply(matrix(sub_prob[,2:(mainJ+1)],ncol=mainJ),1, function(x) which.max(x))
    sub_hgtype_fin <- rep(NA, n)
    if (length(ok_subzoneh_spl)) {
      for(ss in 1:length(ok_subzoneh_spl)){
        sub_reg <- ok_subzoneh_spl[[ss]]
        if (length(sub_reg)) {
          sub_hgtype_fin[sub_reg] <- sub_hgtype_dec[ss] 
          temp0   <- (tauJ*(unique(main_hgtype[sub_reg])-1)+1):(tauJ*(unique(main_hgtype[sub_reg])-1)+tauJ)
          state_i <- temp0[which(sub_zk[temp0]==main_zk[sub_hgtype_dec[ss]])]+mainJ
          if (length(state_i)) {
            hat_logr[sub_reg]  <- mu0_lrs[state_i]
            hat_logor[sub_reg] <- mu0_lors[state_i]
          }
        }
      }
    }
    if (DEBUG) print("END: final copy number")
  }

  prm1          <- getParmNames(1)
  prm0          <- getParmNames(0)
  theta1        <- c(lpsi0,lpr0,lsigma20_lr,lsigma20_lor,v0)
  names(theta1) <- prm1
  se_theta      <- NULL
  asymse_varlr  <- NULL
  if (DEBUG) print("BEGIN: hessian")
  hes_theta1    <- hessian(hes_logL1, theta1)
  if (DEBUG) print("END: hessian")
  var_theta1    <- try(solve(-hes_theta1), silent=TRUE)
  if ("try-error" %in% class(var_theta1)) {
    warning("Hessian matrix could not be inverted, covariances will not be computed")
    var_theta1  <- NULL
    var_theta   <- NULL
    asymvar1_lr <- NULL
  } else {
    theta1_dev   <- diag(c(2^lpsi0*log(2), -(1+exp(lpr0))^-2*exp(lpr0),
                         exp(lsigma20_lr),exp(lsigma20_lor),1))
    var_theta    <- t(theta1_dev)%*%var_theta1%*%theta1_dev
    se_theta     <- sqrt(diag(var_theta)) # asymtotic standard errors of the global estimators given data and lsz0
    var1_s2lr    <- diag(var_theta)[3]
    var1_s2v     <- diag(var_theta)[5]

    if (v0 < 2) {
      warning("df was estimated to be less than two, asymptotic variance of logR will not be computed")
      asymvar1_lr <- NULL 
    } else {
      asymvar1_lr  <- 1/(v0-2)^2*((v0^2*var1_s2lr)-4*v0/(v0-2)*exp(lsigma20_lr)*var_theta[3,5]+4/((v0-2))^2*(exp(lsigma20_lr))^2*var1_s2v)
      asymse_varlr <- sqrt(asymvar1_lr)
    }
    rownames(var_theta) <- prm0
    colnames(var_theta) <- prm0
    names(se_theta)     <- prm0
  }
 
  theta.orig             <- c(psi0, pr0, sigma20_lr, sigma20_lor, v0)
  names(theta.orig)      <- getParmNames(0)
  hatp0_m                <- matrix(hatp0_m, nrow=2)
  rownames(hatp0_m)      <- c("P(no subclone)", "P(subclone)")
  colnames(hatp0_m)      <- rownames(hatp0_m)
  rownames(hatp0_M)      <- genoStates
  colnames(hatp0_M)      <- genoStates
  rownames(hatp0_mgtype) <- genoStates
  colnames(hatp0_mgtype) <- genoStates[-1]
  tmp <- paste("P(M=", genoStates, "|no S)", sep="")
  for (i in 1:mainJ) tmp <- c(tmp, paste("P(S=", genoStates[-1], "|M=", genoStates[i], ")", sep=""))
  colnames(pzks1)        <- tmp

  list(converged=converged, parms=theta.orig, 
       cov.parms=var_theta, se.parms=se_theta, asym.se.varlogR=asymse_varlr,
       logR.est=hat_logr, logOR.est=hat_logor,
       prob.stage1=pzks1, mainclone.genotype.index=main_hgtype, 
       subclone.prob=sub_prob, subclone.ind=subc_hidx,  
       subclone.genotype=sub_hgtype_fin, subclone.regions=ok_subzoneh_spl,
       clonal.prop.region=sz_ests, clonal.prop.est=sz0, 
       loglike.vec=logL_Y[1:mm],
       mainclone.trans.est=hatp0_M, subclone.trans.est=hatp0_m, 
       subclone.probMat.est=hatp0_mgtype
       ) 

} # END: subhmm

getParms.origScale <- function(PARM1) {

  ret    <- PARM1
  ret[2] <- 1/(1+exp(PARM1[2]))
  ret[1] <- 2^PARM1[1]
  ret[3] <- 1/(1+exp(PARM1[3]))
  ret[4] <- exp(PARM1[4])
  ret[5] <- exp(PARM1[5])

  ret
 
} # END: getParms.origScale

getParmNames <- function(which, all=0) {

  if (all) {
    if (!which) {
      # Orignal scale
      ret <- c("ploidy", "purity", "clonal.prop", "logR.var", "logOR.var", "df")
    } else {
      ret <- c("logBase2Ploidy", "logOneMinusPurityOverPurity", 
              "logOneMinusClonalpOverClonalp",
              "logLogR.var", "logLogOR.var", "df")
    }
  } else {
    if (!which) {
      # Orignal scale
      ret <- c("ploidy", "purity", "logR.var", "logOR.var", "df")
    } else {
      ret <- c("logBase2Ploidy", "logOneMinusPurityOverPurity", 
              "logLogR.var", "logLogOR.var", "df")
    }
  }

  ret

} # END: getParmNames

checkStop <- function(logL0, logL1, parm0, parm1, loglikeTol, parmTol, iter, print) {

  v1  <- abs(logL0 - logL1)
  v2  <- max(abs(parm0 - parm1))
  if (print) {
    str <- paste("Iter ", iter, " loglike = ", formatValue(logL1, digits=4),
                                " loglike diff = ", formatValue(v1),
                                " max parm diff = ", formatValue(v2), "\n", sep="")
    cat(str)
    if (print > 1) {
      cat("Estimates:\n")
      parm1 <- getParms.origScale(parm1) 
      tmp   <- formatValue(parm1)
      names(tmp) <- getParmNames(0, all=1)
      print(tmp)
    }
  }
  ret <- (v1 < loglikeTol) || (v2 < parmTol)

  ret

} # END: checkStop

formatValue <- function(val, digits=4) {
  as.numeric(formatC(val, format="f", digits=digits))
} # END: formatValue

# Function to remove leading/trailing white space
removeWhiteSpace <- function(str, leading=1, trailing=1) {

  if ((leading) && (trailing)) {
    ret <- gsub("^\\s+|\\s+$", "", str, perl=TRUE)
  } else if (leading) {
    ret <- gsub("^\\s+", "", str, perl=TRUE)
  } else if (trailing) {
    ret <- gsub("\\s+$", "", str, perl=TRUE)
  } else {
    ret <- str
  }

  ret

} # END: removeWhiteSpace

# Function to compute a cummaulative sum in "reverse order"
reverseCumSum <- function(vec) {

  n   <- length(vec)
  ret <- vec[n:1]
  ret <- cumsum(ret)
  ret <- ret[n:1]

  ret
}


# Function to check the genotype states
checkGenoStates <- function(gs, unq) {

  temp <- nchar(gs) > 0
  gs   <- gs[temp]
  gs2  <- gsub(unq[1], "@", gs, fixed=TRUE)
  gs2  <- gsub(unq[2], unq[1], gs2, fixed=TRUE)
  gs2  <- gsub("@", unq[2], gs2, fixed=TRUE)
    
  temp <- gs2 %in% gs
  if (any(temp)) {
    err <- gs2[temp]
    str <- paste(err, collapse=",", sep="")
    str <- paste("ERROR: with genoStates ", str, sep="")
    stop(str)
  }

  NULL

} # END: checkGenoStates

# Function to compute copy number, etc
callGenoStates <- function(zk) {

  # Initially set "0" to "" and make it be first
  zk    <- unique(removeWhiteSpace(zk))
  zk    <- gsub("0", "", zk, fixed=TRUE)
  tmp   <- nchar(zk) > 0
  zk    <- zk[tmp]
  zk    <- c("", zk)
  J     <- length(zk)
  ctz0  <- nchar(zk)
  tlist <- strsplit(zk, "", fixed=TRUE)
  vec   <- unlist(tlist)
  unq   <- unique(vec)
  temp  <- nchar(unq) > 0
  unq   <- unq[temp]
  nunq  <- length(unq)
  if ((!nunq) || (nunq > 2)) stop("ERROR with genoStates")
  if (nunq == 1) unq <- c(unq, unq)
  mz    <- rep(NA, J)
  pz    <- rep(0, J)

  # Determine the major/minor allele. We want the major to be first
  n1 <- sum(vec %in% unq[1])
  n2 <- sum(vec %in% unq[2])
  if (n2 > n1) unq <- unq[2:1]

  # Compute mz, pz and also normalize the genotype states
  for (j in 1:J) {
    vec   <- tlist[[j]]
    mz[j] <- sum(vec %in% unq[1]) 
    if (nunq > 1) pz[j] <- sum(vec %in% unq[2])
    zk[j] <- paste(sort(vec), collapse="", sep="") 
  } 
  if (length(unique(unq)) > 1) checkGenoStates(zk, unq) 

  # Reset "" to "0"
  zk[1] <- "0"

  list(zk=zk, ctz0=ctz0, mz=mz, pz=pz, unq_alleles=unq, gs.null="0")

} # END: callGenoStates




