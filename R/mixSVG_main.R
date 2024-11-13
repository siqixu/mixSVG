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

eps = rep(0,n)
beta = -6
tau = 1
offset = log(libsize)
# iteration 
for(iter in 1:100) {
  eta = X %*% beta  + eps + offset
  eta = as.vector(eta)
  mu = exp(eta)
  #mu[is.infinite(mu)] = median(mu)
  w = (1/mu)*(y-mu) + eta - offset
  vw = 1/mu + tau
  #w[abs(w)>1000] = median(w)
  #vw[vw>1000] = median(vw)
  
  res = (w - X %*% beta)/vw
  eps = tau*res    
  
}
XVivX_iv = solve(t(X/vw)%*%X)
 res2 = res^2
res2_perm = matrix(res2[perm_sample], nrow = nrow(perm_sample))


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
            Tv_perm = colSums(res2_perm * s_sq)
            ETv = mean(Tv_perm)
            DTv = var(Tv_perm)

  n = length(y)
  J = rep(1,n)
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
  DTv = 2*trPKPK - 2*trPKP^2/trPP
            
   
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
        pattern = pat_name, res = res, vw = vw, mu = mu, w = w  )
    return(out)
}
