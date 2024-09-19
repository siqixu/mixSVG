
mixSVG = function(count,
                  coord,
                  X = NULL,
                  libsize_inc = TRUE,
                  libsize = NULL,
                  vtest_zero_prop = 0.995,
                  ncore = 10,
                  n_perm = 1000,
                  sig = 0.05,
                 c_gau, c_cos, l_gau, l_cos){

  n = ncol(count)
  ngene = nrow(count)
  
  if(n!=nrow(coord)){
    cat('The matrics count and coord do not have the same numbers of spots')
    break
  }
  if(!is.null(X)){
    if(n!=nrow(X)){
      cat('The matrics count and X do not have the same numbers of spots')
      break
    }
    if(ncol(X)!=2){
      cat('The matrix X should contains two columns for two-dimensioanl spatial coordinates.')
      break
    }
  }
  cat('\n-------------------------------------------------------')
  cat("\nDetecting Spatially Variable (SV) Genes by mixSVG")
  cat('\n-------------------------------------------------------')
  cat('\nNumber of genes:', ngene)
  cat('\nNumber of spots:', n)
  
  if(!is.null(X) & is.null(colnames(X))){colnames(X) = paste0("x", 1:ncol(X))}

  if(libsize_inc & is.null(libsize)){
    libsize = colSums(count)
  }
  if(!libsize_inc){
    libsize = rep(1, n)
  }

  # transformation of spatial coordinates
  s_trans = coord

for(transfunc in c('gaussian', 'cosine')){
  for(q in c(0.2, 0.8)){
    
    if(transfunc=='gaussian'){cx = cy = c_gau
    }else{
      cx_cos_0.2 = c_cos[1,]
      cx_cos_0.8 = c_cos[2,]
      cy_cos_0.2 = c_cos[3,]
      cy_cos_0.8 = c_cos[4,]
      
      if(q==0.2){cx = cx_cos_0.2;cy = cy_cos_0.2}
      if(q==0.8){cx = cx_cos_0.8;cy = cy_cos_0.8}
    }
    
  for(ci in 1:length(cx)){
      
      s_trans = cbind(s_trans, 
                      transcoord_func(coord[,1], transfunc = transfunc, q = q, c = cx[ci]),
                      transcoord_func(coord[,2], transfunc = transfunc, q = q, c = cy[ci]))
 
    }
  }
}


  pat_idx = 1:(ncol(s_trans)/2)

  # generate permutation samples
  perm_sample = apply(t(1:n_perm), 2, FUN = function(i){
    set.seed(i)
    sample(1:n, size = n, replace = F)
  })


  X = cbind("intercept" = rep(1,n), X)
  registerDoParallel(ncore)
  results <- foreach(gi = 1:ngene) %dopar% {
    # gi = 1;
    set.seed(gi)
    y = as.matrix(count[gi,])
    results = mixSVG_main(y, X, s_trans, pat_idx, perm_sample, libsize, vtest_zero_prop)
  }
  names(results) = rownames(count)

  pval = unlist(lapply(results, function(x){x$pval}))
  pval_adj = p.adjust(pval, method="BH")

  pval_all = cbind.data.frame(pval, pval_adj)
  pval_sig = pval_all[pval_adj < sig,]

  mixSVG_output = list(results=results, pval_all=pval_all, pval_sig=pval_sig)

  cat('\nNumber of detected SV genes:', nrow(mixSVG_output$pval_sig),
      '\n(significance level for adjusted P-values:', sig, ')' )
  cat('\n-------------------------------------------------------\n')
  return(mixSVG_output)

}








