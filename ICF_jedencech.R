library(MASS)



y_ICF = ICF$ICF


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
mme(y_ICF, X, Z, A, 235.55, 412.24)

var(y_ICF)
sigma_a = 121.01    #starting value for random effect
sigma_e = 10.51 #starting value for error variance

# standaryzacja Min-Max, inaczej nie dało się ustalić nowych sigm funkcją EM. To mi chat podpowiedział wiec do weryfikacji
y_ICF_norm = (as.numeric(y_ICF) - min(as.numeric(y_ICF))) / (max(as.numeric(y_ICF)) - min(as.numeric(y_ICF)))

EM = function(y, X, Z, A, sigma_a, sigma_e) {
  n = nrow(X)
  p = ncol(X)
  q = nrow(A)
  
  t = 1 #iteration number 1
  tmp = 0.1 #test for convergance
  
  while (tmp > 0.0001) {
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
  list(t = t, sigma_a = sigma_a, sigma_e = sigma_e)
}

(wyniki = EM(y_ICF_norm, X, Z, A, sigma_a, sigma_e))
wyniki$sigma_a+wyniki$sigma_e
var(y_ICF_norm)

#odziedziczalność

h2_ICF = (0.02629626 / (0.02629626 + 0.008426684)) * 100
h2_ICF

est_ICF = mme(y_ICF, X, Z, A, 0.02629626, 0.008426684)$est[19:499]
(m = mean(est_ICF))
(s = sd(est_ICF))

#standaryzacja
est_ICF = 10*(est_ICF - m) / s + 100
mean(est_ICF)
sd(est_ICF)
est_ICF


C = as.matrix(mme(y_ICF, X, Z, A, 0.02629626, 0.008426684)$C)
invC_ICF = ginv(C)

invC22_ICF = invC_ICF[19:499, 19:499]
(r2_ICF = diag(1 - invC_ICF*(0.008426684/0.02629626)))
r2_ICF[r2 < 0] = 0
(r_ICF = sqrt(r2_ICF))



