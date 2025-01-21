library(MASS)



y_ICF = ICF


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
mme_icf = mme(y_ICF, X, Z, A, g22, z22)
mme_icf


# ICF
random_ICF = mme_icf$est[20:length(mme_icf$est)]
random_observed_ICF = random_ICF[321:499]



#standaryzacja
ICF_std = 10*(ICF- mean(ICF)) / sd(ICF) + 100
random_observed_ICF_std = 10*(random_observed_ICF - mean(random_observed_ICF)) / sd(random_observed_ICF) + 100


C = as.matrix(mme(y_ICF, X, Z, A, g22, z22)$C)
invC_ICF = ginv(C)

invC22_ICF = invC_ICF[19:499, 19:499]
(r2_ICF = diag(1 - invC22_ICF*(as.numeric(z22)/as.numeric(g22))))
r2_ICF[r2_ICF < 0] = 0
(r_ICF = sqrt(r2_ICF))

plot(r_ICF, ylim = c(0, 1))



plot(ICF_std, random_observed_ICF_std, main = 'Obserwowane ICF vs estymowana wartosc hodowlana')

