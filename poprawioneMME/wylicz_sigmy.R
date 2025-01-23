library(MASS)



y = IFL


mme = function(y, X, Z, A, sigma_a, sigma_e) {
  alpha = sigma_e / sigma_a
  invA = ginv(A)
  C = rbind(cbind(t(X)%*%X, t(X)%*%Z),
            cbind(t(Z)%*%X, t(Z)%*%Z+invA*c(alpha)))
  rhs = rbind(t(X)%*%y, t(Z)%*%y)
  invC = ginv(C)
  estimators = invC%*%rhs
  list(C = C, est = estimators)
}


var(y)
sigma_a = 121.01    #starting value for random effect
sigma_e = 10.51 #starting value for error variance


EM = function(y, X, Z, A, sigma_a, sigma_e) {
  n = nrow(X)
  p = ncol(X)
  q = nrow(A)

  
  t = 1 #iteration number 1
  tmp = 0.1 #test for convergance
  
  while (tmp > 0.05) {
    print(paste("Iteration:", t, "sigma_a:", sigma_a, "sigma_e:", sigma_e, "tmp:", tmp))
    mme_new = mme(y, X, Z, A, sigma_a, sigma_e)
    C_new = ginv(mme_new$C)
    Ck = C_new[(p+1):(p+q), (p+1):(p+q)]
    mme2 = mme_new$est
    
    a = as.matrix(mme2[(p+1):(p+q)])
    sigma_a_new = (t(a)%*%ginv(A)%*%a + sum(diag(ginv(A)%*%Ck))*c(sigma_e))/q
    
    res = as.matrix(y-X%*%as.matrix(mme2[1:p]) - Z%*%as.matrix(mme2[(p+1):(p+q)]))
    X.tmp1 = cbind(X,Z) %*% C_new
    X.tmp2 = t(cbind(X,Z))
    sigma_e_new = (t(res)%*%res + sum(diag(X.tmp1%*%X.tmp2))*c(sigma_e))/n
    
    tmp = max(abs(sigma_a - sigma_a_new), abs(sigma_e - sigma_e_new))
    sigma_a = sigma_a_new
    sigma_e = sigma_e_new
    
    t = t + 1
  }
  list(t = t, sigma_a = sigma_a, sigma_e = sigma_e, res = res)
}

wyniki_IFL = EM(y, X, Z, A, sigma_a, sigma_e)
wyniki_IFL

#simga_a = 6277.935, sigma_e = 55.59866
res_IFL <- wyniki_IFL$res

g11 = wyniki_IFL$sigma_a
z11 = wyniki_IFL$sigma_e

#-----------------------------------



y = ICF


mme = function(y, X, Z, A, sigma_a, sigma_e) {
  alpha = sigma_e / sigma_a
  invA = ginv(A)
  C = rbind(cbind(t(X)%*%X, t(X)%*%Z),
            cbind(t(Z)%*%X, t(Z)%*%Z+invA*c(alpha)))
  rhs = rbind(t(X)%*%y, t(Z)%*%y)
  invC = ginv(C)
  estimators = invC%*%rhs
  list(C = C, est = estimators)
}




var(y)
sigma_a = 121.01    #starting value for random effect
sigma_e = 10.51 #starting value for error variance


EM = function(y, X, Z, A, sigma_a, sigma_e) {
  n = nrow(X)
  p = ncol(X)
  q = nrow(A)
  
  t = 1 #iteration number 1
  tmp = 0.1 #test for convergance
  
  while (tmp > 0.05) {
    print(paste("Iteration:", t, "sigma_a:", sigma_a, "sigma_e:", sigma_e, "tmp:", tmp))
    mme_new = mme(y, X, Z, A, sigma_a, sigma_e)
    C_new = ginv(mme_new$C, tol = 1e-20)
    Ck = C_new[(p+1):(p+q), (p+1):(p+q)]
    mme2 = mme_new$est
    
    a = as.matrix(mme2[(p+1):(p+q)])
    sigma_a_new = (t(a)%*%ginv(A)%*%a + sum(diag(ginv(A)%*%Ck))*c(sigma_e))/q
    
    res = as.matrix(y-X%*%as.matrix(mme2[1:p]) - Z%*%as.matrix(mme2[(p+1):(p+q)]))
    X.tmp1 = cbind(X,Z) %*% C_new
    X.tmp2 = t(cbind(X,Z))
    sigma_e_new = (t(res)%*%res + sum(diag(X.tmp1%*%X.tmp2))*c(sigma_e))/n
    
    tmp = max(abs(sigma_a - sigma_a_new), abs(sigma_e - sigma_e_new))
    sigma_a = sigma_a_new
    sigma_e = sigma_e_new
    
    t = t + 1
  }
  list(t = t, sigma_a = sigma_a, sigma_e = sigma_e, res = res)
}



wyniki_ICF = EM(y, X, Z, A, sigma_a, sigma_e)

wyniki_ICF

#simga_a = 6017.857, sigma_e = 55.4629405

res_IFL <- wyniki_IFL$res
res_ICF <- wyniki_ICF$res


g11 = wyniki_IFL$sigma_a
z11 = wyniki_IFL$sigma_e

g22  = wyniki_ICF$sigma_a
g12 = cov(IFL, ICF)

z22 = wyniki_ICF$sigma_e
z12 = cov(res_IFL,res_ICF)



#G = rbind(c(6434.149, g12), c(g12, 37.19625))
G = rbind(c(g11, g12), c(g12, g22))
G

#R = rbind(c(53.38752, z12),c(z12, 4.768613))
R = rbind(c(z11, z12),c(z12, z22))
R




