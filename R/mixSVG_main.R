mixSVG_main <- function (y, X, s_trans, pat_idx, pat_name, perm_sample, libsize, vtest_zero_prop){
    vtest = (mean(y == 0) < vtest_zero_prop)
    model_init = glm(y ~ X - 1 + offset(log(libsize)), family = poisson)
    model0 = fit_glmm(y, X, model_init, libsize)
    par = model0$par[1, ]
    beta = par[1:ncol(X)]
    tau = par['tau']
    w = model0$w
    vw = model0$vw
    mu =  model0$mu
    XVivX_iv = solve(t(X/vw)%*%X)
    res = (w - X %*% beta)/vw
    res2 = res^2
    res2_perm = matrix(res2[perm_sample], nrow = nrow(perm_sample))


if (vtest) {
  
  res2_perm = Vw_perm = numeric()
 
  for(i_perm in 1:ncol(perm_sample)){
    
    eps_perm =  rnorm(length(y),0,sqrt(tau)) # tau*res_perm[,i_perm]
    eta_perm = beta + eps_perm + log(libsize)  # as.vector(X %*% beta) 
    mu_perm = exp(eta_perm)
    
    y_perm = rpois(length(y), mu_perm)
    model0 = fit_glmm(y_perm, X, model_init, libsize)
    
    beta_perm = model0$par[1, ][1:ncol(X)]
    
    w_perm = model0$w
    vw_perm = model0$vw
    res2_perm = cbind(res2_perm,((w_perm - X %*% beta_perm)/vw_perm)^2)
    Vw_perm = cbind(Vw_perm,vw_perm)
  }
  
  
}


    
    test_func = function(i_pat) {
        s1 = s_trans[, (2 * i_pat - 1)]
        s2 = s_trans[, (2 * i_pat)]
        s = cbind(s1, s2)
        s1_sq = s1^2
        s2_sq = s2^2
        Tb = c(sum(res * s1), sum(res * s2))
        ZivX = t(s/vw) %*% X
        Vb = t(s/vw) %*% s - ZivX %*% XVivX_iv %*% t(ZivX)
        Tb = t(Tb) %*% solve(Vb) %*% t(t(Tb))
        pval_b = pchisq(Tb, 2, lower.tail = F)
      
        if (vtest) {
            s_sq = s1_sq + s2_sq
            Tv = sum(res2 * s_sq)
            #Tv_perm = colSums(res2_perm * s_sq)
            #ETv = mean(Tv_perm)
            #DTv = var(Tv_perm)


n = length(y)
J = rep(1,n)

moment = function(s1,s2,vw){
  JVinvJ=sum(1/vw)
  JVinv.X1=sum(s1/vw)
  A1=(s1^2+s2^2)/vw
  JVinvA1.J=sum(A1/vw)
  
  XVinX =sum(1/vw)
  XVin2X =sum(1/vw^2)
  XVin3X =sum(1/vw^3)
  XVin2XK =sum((s1^2+s2^2)/vw^2)
  XVin3XK =sum((s1^2+s2^2)/vw^3)
  
  trPK = sum(A1) - sum((s1^2+s2^2)/vw^2)/XVinX
  trPP = sum(1/vw^2) - 2*XVin3X/XVinX + (XVin2X/XVinX)^2
  trPKP = XVin2XK -2*XVin3XK/XVinX + XVin2X*XVin2XK/XVinX^2
  trPKPK  = sum(A1^2)-2*sum(A1^2/vw)/JVinvJ + (JVinvA1.J/JVinvJ)^2

  
  ETv = trPK
  DTv = 2*trPKPK  - 2*trPKP^2/trPP
  
  return(c(ETv,DTv))
}


mm = moment(s1,s2,vw)
ETv = mm[1]
DTv = mm[2]    


mm_perm = apply(Vw_perm, 2, moment, s1=s1, s2=s2)
ETv_perm =  mean(mm_perm[1,])
DTv_perm =  mean(mm_perm[2,])


   Tv_perm = colSums(res2_perm * s_sq)         
  ETv0_perm = mean(Tv_perm)
  DTv0_perm = var(Tv_perm)
            
ETv = ETv - mean(ETv_perm) + ETv0_perm
DTv = DTv - mean(DTv_perm) + DTv0_perm

            
            k = DTv/(2 * ETv)
            df = 2 * ETv^2/(DTv)
            pval_v = c(pchisq(Tv/k, df, lower.tail = FALSE), pchisq(Tv/k, df, lower.tail = TRUE))
            pval_v = 2*min(pval_v)
            
            pval = c(pval_b, pval_v)
            pval[which(pval == 0)] <- 5.55e-17
            pval[which((1 - pval) < 0.001)] <- 0.99
            T_omn = mean(tan(pi * (0.5 - pval)))
            pval = 1 - pcauchy(T_omn)
        }
        else {
            pval_v = 1
            pval = pval_b
            Tv =  k = df = 0
        }
        return(c(pval, pval_b, pval_v,  Tv, k, df))
    }
    pval_pat = t(apply(t(pat_idx), 2, FUN = test_func))
    colnames(pval_pat) = c("pval_omn", "pval_b", "pval_v",  'Tv', 'k', 'df')
    pval = pval_pat[, "pval_omn"]
    pval[which(pval == 0)] <- 5.55e-17
    pval[which((1 - pval) < 0.001)] <- 0.99
    T_final = mean(tan(pi * (0.5 - pval)))
    pval = 1 - pcauchy(T_final)
    out = list(model0 = par, pval = pval, pval_pat = pval_pat, 
        pattern = pat_name, res = res, vw = vw, mu = mu, w = w, Tv = Tv)
    return(out)
}
