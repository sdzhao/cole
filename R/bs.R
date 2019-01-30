lambda.tuning = function(x,y,beta)

bs.camel = function(x,y, lambda_seq = NULL, nlambda = 20, verbose = F, x_val = NULL, y_val = NULL, folder = 10, intercept = F){
  if(!is.null(x_val) && !is.null(y_val)){
    if(verbose){
      print("tuning lambda with validation data sets")
    }
    
    if(is.null(lambda_seq)){
      if(verbose){
        print(paste("lambda sequnce not provided, tuning with nlambda = ", nlambda, sep = ""))
      }
      fit_camel = camel.cmr(x, y, nlambda = nlambda, verbose = verbose)
    }else{
      if(verbose){
        print(paste("lambda sequnce provided with length ", length(lambda_seq), sep = ""))
      }
      fit_camel = camel.cmr(x, y, lambda = lambda_seq, verbose = verbose)
      nlambda = length(lambda_seq)
    }
    
    err_temp = numeric(nlambda)
    
    
    
    
    
    
    for(i in 1:nlambda){
      beta_temp = fit_camel$beta[[i]]
      esta0 = m_y - as.vector(t(as.matrix(colMeans(X_train))) %*% beta_temp)
      esta0_mat = matrix(rep(esta0, n_val), nrow = n_val, byrow = T)
      err_temp[i] = norm(Y_val - esta0_mat - X_val%*%beta_temp,'F')
    }
    orders = sort(err_temp, index.return = T)
    
    esta0 = m_y - as.vector(t(as.matrix(colMeans(X_train))) %*% fit_camel$beta[[orders$ix[1]]])
    esta0_mat = matrix(rep(esta0, n_te), nrow = n_te, byrow = T)
    PE_mat[j,1] = norm(Y_test-esta0_mat-X_test%*%fit_camel$beta[[orders$ix[1]]],'F')^2/n_te
    time_mat[j,1] = (proc.time() - t)["elapsed"]
    
    beta_camel = matrix(NA, nlambda*q, p+1)
    sigma_camel = numeric(nlambda*q)
    for(i in 1:nlambda){
      beta_camel[((i-1)*q+1):(i*q),2:(p+1)] = t(fit_camel$beta[[i]])
      beta_camel[((i-1)*q+1):(i*q),1] = m_y - as.vector(t(as.matrix(colMeans(X_train))) %*% fit_camel$beta[[i]])
      # beta_camel[((i-1)*q+1):(i*q),1] = 0
      sigma_camel[((i-1)*q+1):(i*q)] = colSums((Y_train - X_train%*%fit_camel$beta[[i]])^2)/n_tr
    }
    
  }
}