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
mme_ifl = mme(y, X, Z, A, g11, z11)
mme_ifl





# IFL

random_observed_IFL = mme_ifl$est[340:518]

#standaryzacja do porownania
IFL_std = 10*(IFL- mean(IFL)) / sd(IFL) + 100
random_observed_IFL_std = 10*(random_observed_IFL - mean(random_observed_IFL)) / sd(random_observed_IFL) + 100



C = as.matrix(mme(y, X, Z, A,g11 , z11)$C)
invC = ginv(C)

invC22 = invC[19:518, 19:518]

(r2 = diag(1 - invC22*(as.numeric(z11)/as.numeric(g11))))
r2[r2 < 0] = 0
(r = sqrt(r2))

r_IFL = r

plot(r_IFL, ylim = c(0, 1), main = 'Dokładność oszacowania dla IFL')
abline(v = 320, col = "red", lwd = 2, lty = 2) 

#plot(IFL_std, random_observed_IFL_std,main = 'Obserwowane IFL vs estymowana wartosc hodowlana')

# Dopasowanie regresji liniowej dla IFL
model_IFL <- lm(IFL_std ~ random_observed_IFL_std)  # Regresja liniowa
cor_IFL <- cor(random_observed_IFL_std, IFL_std)    # Korelacja

# Wykres z linią trendu dla IFL
plot(random_observed_IFL_std, IFL_std, 
     main = "Obserwowane IFL vs estymowana wartość hodowlana",
     ylim = c(90, 160), 
     xlab = "Estymowana wartość IFL (standaryzowana)", 
     ylab = "Obserwowana wartość IFL (standaryzowana)")
abline(model_IFL, col = "red", lwd = 2)  # Dodanie linii trendu w kolorze czerwonym

