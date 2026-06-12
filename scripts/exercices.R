# Statistiques Bayésiennes
# 2026-06-11

# Packages ----

library(bayess)
library(rjags)
library(BayesFactor)
library(lars)
# library(ggstatsplot)

# Exercice 1 - normaldata ----

data("normaldata")
head(normaldata)

shift <- normaldata$x2

hist(shift, nclass = 10, col = "steelblue", probability = TRUE, main = "")

# Vraisemblance des données de `shift` en fonction de la valeur `m` de l'espérance
m <- seq(-0.6, 0.6, 0.01)
L <- NULL

for (i in 1:length(m)) { 
  L[i] <- prod(dnorm(shift, mean = m[i], sd = sd(shift)))
}

plot(m, L, type = "l", ylab = "vraisemblance",
     main = "vraisemblance en fonction de la moyenne m")

plot(m, L, type = "l", ylab = "density", main = "vraisemblance et posterior")
n <- length(shift)
d <- dnorm(m, n * mean(shift)/(n+1), sqrt(var(shift)/(n+1)))
lines(m, d/max(d) * max(L), col = 2)

cat('Maximum de vraisemblance: ', mean(shift), '\n')
cat('Maximum a posteriori: ', n*mean(shift)/(n+1), '\n')

# On suppose que la variance du prior est la même que celle des données.

# Impact de la presence d'une pisciculture ----

ForkLengthData <- read.csv2("data/ForkLengthData.csv")

ForkLengthData$Location <- as.factor(ForkLengthData$Location)
levels(ForkLengthData$Location) <- c("farm", "downstream",
                                     "other", "upstream")

boxplot(ForkLengthData$ForkLengthinmillimeters ~ ForkLengthData$Location)

# Approche fréquentiste
Long_aval <- ForkLengthData[ForkLengthData$Location == "downstream", 2]
Long_amont <- ForkLengthData[ForkLengthData$Location == "upstream", 2]

par(mfrow = c(1, 2))
hist(Long_amont); hist(Long_aval)

t.test(Long_aval, Long_amont, alternative = "greater")
wilcox.test(Long_aval, Long_amont, alternative = "greater")

tmp <- t.test(Long_aval)
cat("Estimation de l'espérance en aval : ", tmp$estimate, "\n")
cat("Intervalle de confiance : [", tmp$conf.int[1], ',', tmp$conf.int[2], "]\n")

tmp <- t.test(Long_amont)
cat("Estimation de l'espérance en amont : ", tmp$estimate, "\n")
cat("Intervalle de confiance : [", tmp$conf.int[1], ',', tmp$conf.int[2], "]\n")

# Approche bayésienne

extra=ForkLengthData[ForkLengthData$Location=="other",]
m=NULL
for (site in unique(extra$StationNumber)){
  m=c(m,mean(extra[extra$StationNumber==site,2]))
}
hist(m,probability = T)
mu=mean(m)
tau2=var(m)
cat('Moyenne des moyennes des tailles sur les autres sites :',mu,'\n')

## Moyenne des moyennes des tailles sur les autres sites : 97.0097

cat('Variance des moyennes des tailles sur les autres sites :',tau2,'\n')

## Variance des moyennes des tailles sur les autres sites : 70.53082
x=seq(min(m)-5,max(m)-5,0.01)
lines(x,dnorm(x,mu,sqrt(tau2)),col=2)

# Espérance et variance de la loi gaussienne a posteriori :
nd=length(Long_aval)
mu_dpost=(mean(Long_aval)*tau2 + mu*var(Long_aval)/nd)/(tau2 + var(Long_aval)/nd)
sigma_dpost=(tau2 * var(Long_aval)/nd)/(tau2 + var(Long_aval)/nd)
cat('Espérance a posteriori, aval :', mu_dpost,'\n')
## Espérance a posteriori, aval : 116.2898
cat('Variance a posteriori, aval :', sigma_dpost,'\n')
## Variance a posteriori, aval : 5.900625

# Même chose pour l'aval :
nu=length(Long_amont)
mu_upost=(mean(Long_amont)*tau2 + mu*var(Long_amont)/nu)/(tau2 + var(Long_amont)/nu)
sigma_upost=(tau2 * var(Long_amont)/nu)/(tau2 + var(Long_amont)/nu)
cat('Espérance a posteriori, amont :', mu_upost,'\n')
## Espérance a posteriori, amont : 104.4493
cat('Variance a posteriori, amont :', sigma_upost,'\n')
## Variance a posteriori, amont : 6.20247

# Représenter l'estimation fréquentiste et bayésienne
x=seq(90,130,0.01)
plot(x,dnorm(x,mu_upost,sqrt(sigma_upost)),type='l',col=2)
abline(v=mean(Long_amont),col=2)
tmp=t.test(Long_amont);
abline(v=tmp$conf.int[1],col=2,lty=2)
abline(v=tmp$conf.int[2],col=2,lty=2)
lines(x,dnorm(x,mu_dpost,sqrt(sigma_dpost)),type='l')
abline(v=mean(Long_aval))
tmp=t.test(Long_aval);
abline(v=tmp$conf.int[1],lty=2)
abline(v=tmp$conf.int[2],lty=2)

# Loi a posteriori de la différence
x=seq(-20,10,0.01)
plot(x,dnorm(x,mu_upost-mu_dpost,sqrt(sigma_upost+sigma_dpost)),type='l',col=2)

# Proba que la taille des poissons en aval
# soit plus grande d'au moins 5mm qu'en
# amont :
pnorm(-5,mu_upost-mu_dpost,sqrt(sigma_upost+sigma_dpost))

# Exercice poids de naissance ----

data <- read.table("data/poidsnaissance.txt", header = TRUE, sep = ",",
                    row.names = 1)

data$OBS <- NULL

sexe <- data$SEXE + 1
poids <- data$POIDNAIS

boxplot(poids ~ sexe)

t.test(poids ~ sexe, alternative = "greater")
wilcox.test(poids ~ sexe, alternative = "greater")

ttestBF(poids[sexe == 1], y = poids[sexe == 2], nullInterval = c(0, Inf))
log10(4.525689)

ttestBF(poids[sexe == 1], y = poids[sexe == 2], 
        nullInterval = c(0, Inf), rscale = "ultrawide")

fumeuse <- data$CIGJOUR > 0
boxplot(poids ~ fumeuse)
t.test(poids ~ fumeuse, alternative = "greater")
ttestBF(poids[fumeuse], poids[!fumeuse], nullInterval = c(-Inf, 0))

age <- data$AGEGEST
boxplot(age ~ fumeuse)
t.test(age ~ fumeuse, alternative = "greater")

# Régression linéaire bayésienne ----

data("caterpillar")
?caterpillar
y <- log(caterpillar$y)
x <- as.matrix(caterpillar[, 1:8])

# Données de départ
hist(exp(caterpillar$y))
hist(caterpillar$y)
hist(log(caterpillar$y))

x <- scale(x)

# Régression classique
summary(lm(y~x))  # 3 variables significatives

# Régression bayésienne
res1 <- BayesReg(y, x, betatilde = c(0, 0, 0, 0, 0, 0, 0, 0),
                 g = length(y))

# régression linéaire données poids de naissance
data <- read.table("data/poidsnaissance.txt", 
                   header = T, sep = ",", row.names = 1)
head(data)

# régression linéaire classique

m1 <- lm(POIDNAIS ~ ., data = data)
summary(m1)

# Supprimer les variables non significatives (-> 15n de variables)
m2 <- step(m1)
summary(m2)

# régression bayésienne

y <- data$POIDNAIS
x <- as.matrix(data[, -3])
x <- scale(x)
summary(lm(y~x))

m3 <- BayesReg(y, x, betatilde = c(0, 0, 0, 0, 0, 0, 0, 0),
                 g = length(y))

# Ridge ----

x1 <- rnorm(100)
x2 <- rnorm(100, x1, 0.01)
plot(x1, x2)
cor(x1, x2)

y <- 1 + 3 * x1 + 2 * x2 + rnorm(100, 0, 0.1)

lm(y ~ x1 + x2)

y <- 1 + 3 * x1 + 2 * x2 + rnorm(100, 0, 0.1)

lm(y ~ x1 + x2)

# grande variance des paramètres estimés, classique quand
# variables corrélées.
y <- 1 + 3 * x1 + 2 * x2 + rnorm(100, 0, 0.1)
summary(lm(y ~ x1 + x2))

# risque d'avoir sign/nonsign

x <- matrix(rnorm(50), nrow = 5, ncol = 10)
x
cor(x)
plot(x[, c(1, 8)])

x <- matrix(rnorm(5000), nrow = 500, ncol = 10)
cor(x)

# quand beaucoup de variables et peu d'observations, on aura des
# variables faussement corrélées.

# Regression penalisee - poids de naissance ----

data <- read.table("data/poidsnaissance.txt", header = TRUE, sep = ",",
                   row.names = 1)

model_lasso <- lars(
  x = as.matrix(data[, c(2:3, 5:8)]),
  y = data$POIDNAIS,
  type = "lasso",
  trace = FALSE,
  normalize = TRUE)

plot(model_lasso, xvar = "df", plottype = "coeff")

cv <- cv.lars(
  x = as.matrix(data[, c(2:3, 5:8)]),
  y = data$POIDNAIS, 
  K = 48, # leave-one-out
  type = "lasso",
  trace = FALSE,
  plot.it = TRUE, 
  se = TRUE,
  mode = "step",
  normalize = TRUE)

abline(h = cv$cv[5] + cv$cv.error[5])

# + faible: etape 5
cv$cv.error[5]
cv$cv[5] + cv$cv.error[5]
cv$cv  # etape 3: inf au seuil

print(model_lasso$lambda[3])
print(model_lasso$beta[3,])

# RJAGS ----

poids_naissace <- read.table(
  file = "data/poidsnaissance.txt", header = TRUE, 
  sep = ",", row.names = 1)

sexe <- poids_naissace$SEXE + 1
poids <- poids_naissace$POIDNAIS

data <- list(poids = poids, sexe = sexe, N = length(poids))

inits <- list(
  list(moyennes = c(2600, 4000), sigma = 500),
  list(moyennes = c(4500, 2700), sigma = 700),
  list(moyennes = c(4000, 4000), sigma = 300)
)

m1 <- jags.model(file = "data/modelepoidsnaissance.txt",
                 data = data,
                 inits = inits,
                 n.chains = 3)

update(m1, 3000)

mcmc1 <- coda.samples(m1, variable.names = c("moyennes", "sigma"),
                      n.iter = 2000)

mean(mcmc1[[1]][, "moyennes[1]"])

plot(mcmc1)

gelman.diag(mcmc1)

gelman.plot(mcmc1)

autocorr.plot(mcmc1)

summary(mcmc1)

dic.samples(m1, n.iter = 1000)

# On cherche à introduire une variable en plus
# AGEGEST

sexe <- poids_naissace$SEXE + 1
poids <- poids_naissace$POIDNAIS
age <- poids_naissace$AGEGEST

data <- list(
  poids = poids, 
  sexe = sexe,
  age = age,
  N = length(poids))

inits <- list(
  list(moyennes = c(2600, 4000), sigma = 500, delta = 150),
  list(moyennes = c(4500, 2700), sigma = 700, delta = 100),
  list(moyennes = c(4000, 4000), sigma = 300, delta = 200)
)

m2 <- jags.model(file = "data/modelepoidsnaissance_2.txt",
                 data = data,
                 inits = inits,
                 n.chains = 3)

update(m2, 3000)

mcmc2 <- coda.samples(m2, variable.names = c("moyennes", "sigma", "delta"),
                      n.iter = 2000)

dic.samples(m1, n.iter = 1000)
dic.samples(m2, n.iter = 1000)

plot(mcmc2)

gelman.diag(mcmc2)

gelman.plot(mcmc2)

autocorr.plot(mcmc2)

summary(mcmc2)
