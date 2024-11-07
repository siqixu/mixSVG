
mixSVG_main = function(y, X, s_trans, pat_idx, pat_name, perm_sample, libsize, vtest_zero_prop) {

  vtest = (mean(y==0) < vtest_zero_prop)

  # estimation under the null
  model_init = glm(y ~  X - 1 + offset(log(libsize)), family = poisson)
  model0 = fit_glmm(y, X, model_init, libsize)
  par = model0$par[1,]
  beta = par[1:ncol(X)]
  w = model0$w
  vw = model0$vw
  XVivX_iv = solve(t(X/vw)%*%X)
  res = (w - X %*% beta)/vw
  res2 = res^2
  res2_perm = matrix(res2[perm_sample], nrow = nrow(perm_sample))

  tau = par['tau']
  #mu = model0$mu
  res_perm = matrix(res[perm_sample], nrow = nrow(perm_sample))
  eps_perm =  tau*res_perm
  eta_perm = X %*% beta  + eps_perm + log(libsize)
  eta_perm = as.vector(eta_perm)
  mu_perm = exp(eta_perm)
  vw_perm = 1/mu_perm + tau
  y_perm = rpoi(n, libsize*mu_perm)
  w_perm = as.vector((1/mu_perm)*(y_perm-mu_perm) + X %*% beta) + eps_perm
  beta_perm = XVivX_iv %*% t(X/vw_perm) %*% w_perm
  res_perm = (w_perm - X %*% beta_perm)/vw_perm
  res2_perm = res_perm^2
  

  test_func = function(i_pat){
    # i_pat = 1
    s1 = s_trans[, (2*i_pat-1)]
    s2 = s_trans[, (2*i_pat)]
    s = cbind(s1,s2)
    s1_sq = s1^2
    s2_sq = s2^2

    # test for fix effect
    Tb = c(sum(res*s1), sum(res*s2))
    ZivX = t(s/vw)%*%X
    Vb = t(s/vw)%*%s - ZivX%*%XVivX_iv%*%t(ZivX)
    Tb =  t(Tb) %*% solve(Vb) %*% t(t(Tb))
    pval_b = pchisq(Tb, 2, lower.tail = F)

    if(vtest){
      # test for random effect
      s_sq = s1_sq + s2_sq
      Tv = sum(res2*s_sq) 
      Tv_perm = colSums(res2_perm*s_sq) 
   

      ETv = mean(Tv_perm)
      DTv = var(Tv_perm)
      k = DTv/(2 * ETv)
      df = 2*ETv^2/(DTv)
      pval_v = c(pchisq(Tv/k, df, lower.tail = FALSE), pchisq(Tv/k, df, lower.tail = TRUE))
      pval_v = 2*min(pval_v)

      # the omnibus test of mixed effects
      pval = c(pval_b, pval_v)
      pval[which(pval == 0)] <- 5.55e-17
      pval[which((1 - pval) < 0.001)] <- 0.99
      T_omn = mean(tan(pi*(0.5-pval)))
      pval = 1 - pcauchy(T_omn)

    }else{
      pval_v = 1
      pval = pval_b
    }

    return(c(pval, pval_b, pval_v))
  }

  # test for each spatial expression pattern
  pval_pat = t(apply(t(pat_idx), 2, FUN = test_func))
  colnames(pval_pat) = c('pval_omn', 'pval_b', 'pval_v')

  # combine the p-values of all patterns
  pval = pval_pat[, 'pval_omn']
  pval[which(pval == 0)] <- 5.55e-17
  pval[which((1 - pval) < 0.001)] <- 0.99
  T_final = mean(tan(pi*(0.5-pval)))
  pval = 1 - pcauchy(T_final)

  out = list(model0 = par, pval = pval,  pval_pat = pval_pat, pattern = pat_name, res=res)
  return(out)
}

