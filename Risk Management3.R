library(quantmod)
#library(lubridate)
library(PerformanceAnalytics)

#### Log return ####
setwd("O:/FIN 627 - Investment Management/02 Project")
tickers <- read.csv("sym3.csv")
sym <- as.vector(c(tickers[, 2], "SPY"))
list <- lapply(sym, function(x)get(getSymbols(x, from = "2014-01-01", to = "2022-02-04")))#
adPrice <- na.omit(Reduce(merge, lapply(list, function(x)x[, 6])))

R <- na.omit(log(adPrice/lag(adPrice)))
rf <- 0.0075

#### CAPM related performance ####

CAPM <- table.CAPM(R,R$SPY.Adjusted)
CAPM

alpha1 <- CAPM["Alpha",]*252
alpha2 <- CAPM["Annualized Alpha",]
plot(as.numeric(alpha1)~as.numeric(alpha2))
abline(a=0,b=1)


#### Step 1: Efficient Frontier ####
mu <- apply(R, 2, mean)*252
mu <- mu[-51]
R <- R[, -51]
si <- cov(R)*252

#mu1 <- c(tickers[, 3], 0.085)
#mu <- mu1[-51]
##### Global Minimum Portfolio ####
l <-matrix(rep(1, 50), ncol = 1)
I <- diag(50)
w_0 <- (solve(si)%*%l)%*%(1/(t(l)%*%solve(si)%*%l))


##### Efficient Frontier ####
B <- solve(si)%*%(I - l%*%t(w_0))
w_1 <- B%*%mu


w_f <- function(A) {
  w_A <- w_0 + (1/A)*w_1
  mu_A <- t(w_A)%*%mu
  sig_A <- sqrt(t(w_A)%*%si%*%w_A)
  results <- c(w_A,mu_A,sig_A)  
  return(results)
}


A_seq <- seq(2, 100, by = 0.1)
ds <- data.frame(A = A_seq, t(sapply(A_seq,w_f)))

colnames(ds)[2:51] <- sym[-51]
colnames(ds)[52:53] <- c("mu_p", "sigma_p")

plot(mu_p ~ sigma_p, data = ds, col = "red", type = "l", main = c("CAL")
     , xlim = c(0, 0.8), ylim = c(0, 1.5))
points(mu_p ~ sigma_p, data = ds[ds$A == 100,], pch = 20, cex = 1)


#### Step 2: Capital Allocation Line ####
##### Optimal risky Portfolio ####
ds$SR <- (ds$mu_p - rf)/ds$sigma_p

abline(a = rf, b = max(ds$SR), col = "blue", lwd = 1)
points(mu_p ~ sigma_p, data = ds[which.max(ds$SR),], pch = 20, cex = 1)
legend("topright", legend = c("EF", "CAL"), fill = c("red", "blue"))

w_p <- t(ds[which.max(ds$SR), 2:51])
colnames(w_p) <- c("weight")


##### Summary in a table ####
w <- cbind(w_p, mu, sig = apply(R, 2, sd)*sqrt(252))
w <- rbind(w, Portfolio = c(1, ds$mu_p[which.max(ds$SR)], ds$sigma_p[which.max(ds$SR)]))


#others in FRM course
#SR_port <- t(ds[1,2:51])
GMV_port <- t(ds[100,2:51])

#w_conv <- function(a) {
#  w_a <- a*GMV_port + (1-a)*SR_port
#  mu_a <- t(w_a)%*%mu
#  sig_a <- sqrt(t(w_a)%*%si%*%w_a)
#  results <- c(w_a,mu_a,sig_a)
#  return(results)
#}

#a_seq <- seq(0,1,length = 1000)

#ds2 <- data.frame(a = a_seq, t(sapply(a_seq,w_conv)) )
#lines(X51~X52,data = ds2,col = 2,lty = 2)  


#### Step 3: Utility ####

r.o <- function(y)rf + y*(ds$mu_p[which.max(ds$SR)] - rf)
sd.o <- function(y)y*ds$sigma_p[which.max(ds$SR)]

# Utility ####
u <- function(y)r.o(y) - 0.5*A*sd.o(y)^2

###### Maximum Utility ####
u.optimal <- function(A){
  A <- A
  u <- function(y)r.o(y) - 0.5*A*sd.o(y)^2
  optimize(function(y)-u(y), c(0, 2), tol = 0.001)[2]
}

u_seq <- sapply(A_seq, u.optimal)
u_seq <- -unlist(u_seq)
u_seq <- cbind(A_seq, u_seq)
plot(u_seq, type = "l")

###### Optimal complete portfolio ####
y.optimal <- function(A){
  A <- A
  u <- function(y)r.o(y) - 0.5*A*sd.o(y)^2
  optimize(function(y)-u(y), c(0, 2), tol = 0.001)[1]
}

y_seq <- sapply(A_seq, y.optimal)
y_seq <- unlist(y_seq)
y_seq <- cbind(A_seq, y_seq)
plot(y_seq, type = "l")

##### Indifference line ####
A <- 10
umax <- -as.numeric(optimize(function(y)-u(y), c(0, 2), tol = 0.001)[2])

s_seq <- seq(0, 0.8, length = 101)
m_seq <- umax + 0.5*A*s_seq^2
lines(s_seq, m_seq)

y <- as.numeric(optimize(function(y)-u(y), c(0, 2), tol = 0.001)[1])
points(sd.o(y), r.o(y), col = "black", pch = 20, cex = 1)
legend("topright", legend = c("EF", "CAL", "Utility"), fill = c("red", "blue", "black"))




#### Geometric Brownian Motion ####
data <- cbind(mu_hat = w[, 2] + (w[, 3]^2)/2, si_hat = w[, 3])
rownames(data) <- rownames(w)
data <- split(t(data), f = rep(1:51, each = 2))

S0 <- getQuote(sym[-51])$Last
S0 <- 1000000



#w_p <- matrix(rep(1/50, 50), ncol = 1)
#R_p <- R%*%w_p
#si_hat_p <- sd(R_p)*sqrt(252)
#mu_hat_p <- mean(R_p)*252 + (si_hat_p^2)/2

#expected_p <- function(t){
#  St_hat <- S0*exp(mu_hat_p*t)
#  si_t <-S0*sqrt((exp(si_hat_p^2*t)-1)*exp(2*mu_hat_p*t))
#  ex <- c(St_hat, si_t)
#  return(ex)
#}


#### VaR ####
c <- 0.05
t <- 1/252

VaR.theory <- function(x){
  S0*sqrt((exp(x[2]^2*t)-1)*exp(2*x[1]*t))*qt(c, Inf)
}

w <- cbind(w, VaR.daily = sapply(data, VaR.theory))


VaR.compare <- function(t, c = 0.05){
  VaR.theory <- expected_p(t)[2]*qt(c, Inf)
  
  d <- rnorm(100000, mean = expected_p(t)[1], sd = expected_p(t)[2])
  VaR.empirical <- expected_p(t)[1] - quantile(d, c)
  
  VaR.c <- data.frame(VaR.theory, VaR.empirical)
  return(VaR.c)
}

VaR.compare(1/252)
VaR.compare(10/52)





w <- cbind(w, t(CAPM[c(1:2, 5, 9:12), ]))
write.csv(w, file = "MV Risky Portfolio.csv")



#### Current contribution ####
w_p1 <- c(tickers[, 4])
colnames(w_p1) <- c("weight")

mu_A1 <- t(w_p1)%*%mu
sig_A1 <- sqrt(t(w_p1)%*%si%*%w_p1)

mu_A1
sig_A1
sqrt((exp(sig_A1^2)-1)*exp(2*mu_A1 + (sig_A1^2)/2))*qt(c, Inf)

##### optimal portfolio with same sd ####
w_f2 <- function(A) {
  w_A <- w_0 + (1/A)*w_1
  mu_A <- t(w_A)%*%mu
  sig_A <- sqrt(t(w_A)%*%si%*%w_A)
  return(sig_A - sig_A1)
}

A2 <- uniroot(w_f2, c(1, 100))$root
w_A2 <- w_0 + (1/A2)*w_1
mu_A2 <- t(w_A2)%*%mu

w_A2
mu_A2
sqrt((exp(sig_A1^2)-1)*exp(2*mu_A2 + (sig_A1^2)/2))*qt(c, Inf)

points(sig_A1, mu_A2, col = "black", pch = 18, cex = 1.5)
points(sig_A1, mu_A1, col = "red", pch = 9, cex = 1.5)


#### Back test ####
start <- as.Date("2022/02/04")
end <- as.Date("2022/04/29")
#Rb <- R[as.Date(start:end), ]
#realR <- apply(Rb, 2, sum)
#R.s <- exp(realR)-1

pp <- adPrice[as.Date(start:end), ]
rr <- as.vector(tail(pp, 1))/as.vector(head(pp, 1)) - 1
names(rr) <- names(pp)
write.csv(rr, file = "actual returns.csv")

##### Test for Optimal portfolio ####
sum(w_p*rr[-51])


##### Test for Optimal portfolio at the same s d####
sum(w_A2*rr[-51])



