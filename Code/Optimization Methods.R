# Optimization Methods

x = seq(-2,2,0.001)
ex = exp(-x^2)*sin(pi*x^3)
plot(y = ex, x = x, xlim = c(-2, 2),col="blue",cex=0.2,
     main="Cool Looking Function")
abline(h=0,col="gray2")
abline(v=0,col="gray2")

h <- function(x) exp(-x^2)*sin(pi*x^3)

# Newton-Raphson Method 

x = seq(-2,2,0.001)
ex = exp(-x^2)*sin(pi*x^3)
plot(y = -ex, x = x, xlim = c(-2, 2),col="blue",cex=0.2,
     main="Reversed Function")
abline(h=0,col="gray2")
abline(v=0,col="gray2")

h_minus <- function(x) -(exp(-x^2)*sin(pi*x^3))

nlm_simulation <- function(iterlim){
  u <- runif(1,-1,1)
  result_max_nlm <- nlm(h_minus,p=u,iterlim = iterlim)
  nlm_max <- -result_max_nlm$minimum
  nlm_max_x <- result_max_nlm$estimate 
  return(list("Max Nlm"=nlm_max,"Max Nlm X" = nlm_max_x,"U"=u))
}

nlm_simulation(100)

nelm <- nlm_simulation(100)
nlm_max <- nelm$`Max Nlm`
nlm_max_x <- nelm$`Max Nlm X`

par(mfrow=c(1,2))
x = seq(-2,2,0.001)
ex = exp(-x^2)*sin(pi*x^3)
plot(y = ex, x = x, xlim = c(-2, 2),col="blue",cex=0.2,
     main=paste0("Newton-Raphson Method with U =",nelm$U))
abline(h=0,col="gray2")
abline(v=0,col="gray2")
points(x = nlm_max_x, y = nlm_max, col = "purple", bg = "purple", cex = 1, type = "p", pch = 19)
points(x=nelm$U,y = h(nelm$U),col="red",pch=19,bg="red")
text(x=nelm$U,y = h(nelm$U),labels = "Starting Point",pos = 4)
text(x = nlm_max_x, y = nlm_max, labels = paste0("(",round(nlm_max_x,4),",",round(nlm_max,4),")"), pos = 4, offset = 0.5)

# Basic-Monte Carlo Optimization 

bmco_unif <- function(n){
  u <- runif(n,-2,2)
  basic_mc_max <- max(h(u))
  basic_mc_max_x <- u[which.max(h(u))]
  return(list("Maximum Unif Dist" = basic_mc_max, "Maximum X-value Unif Dist"=basic_mc_max_x))
}

bmco_unif(100)

bmco_u <- bmco_unif(100)
bmco_u_max <- bmco_u$`Maximum Unif Dist`
bmco_u_max_x <- bmco_u$`Maximum X-value Unif Dist`

x = seq(-2,2,0.001)
ex = exp(-x^2)*sin(pi*x^3)
plot(y = ex, x = x, xlim = c(-2, 2),col="blue",cex=0.2,
     main="Basic Monte Carlo - Uniform Distribution")
abline(h=0,col="gray2")
abline(v=0,col="gray2")
points(x = bmco_u_max_x, y = bmco_u_max, col = "orange", bg = "orange", cex = 1, type = "p", pch = 19)
text(x = bmco_u_max_x, y = bmco_u_max, labels = paste0("(",round(bmco_u_max_x,4),",",round(bmco_u_max,4),")"), pos = 4, offset = 0.5)


# Basic-Monte Carlo Optimization (Something other than uniform) 

bmco_beta <- function(n){
  alpha=25
  beta=10
  b <- rbeta(n,alpha,beta)
  h_b <- h(b)
  max_h_b <- max(h_b)
  max_h_b_x <- b[which.max(h_b)]
  return(list("Maximum Beta Dist" = max_h_b, "Maximum X-value Beta Dist"=max_h_b_x))
}

bmco_beta(100)

bmco_b <- bmco_beta(n)
bmco_b_max <- bmco_b$`Maximum Beta Dist`
bmco_b_max_x <- bmco_b$`Maximum X-value Beta Dist`

x = seq(-2,2,0.001)
ex = exp(-x^2)*sin(pi*x^3)
plot(y = ex, x = x, xlim = c(-2, 2),col="blue",cex=0.2,
     main="Basic Monte Carlo - Beta Distribution")
abline(h=0,col="gray2")
abline(v=0,col="gray2")
points(x = bmco_b_max_x, y = bmco_b_max, col = "green", bg = "green", cex = 1, type = "p", pch = 19)
text(x = bmco_u_max_x, y = bmco_u_max, labels = paste0("(",round(bmco_b_max_x,4),",",round(bmco_b_max,4),")"), pos = 4, offset = 0.5)
curve(dbeta(x,25,10),add = T,col="red")
legend("topleft", legend = c("Function", "Beta Density"),
       col = c("blue", "red"), lty = 1,cex = 0.5)

# Simulated Annealing 

Simulated_annealing <- function(min_iter,max_iter){
  
  tol <- 1e-20
  temp <- 0.05/log(1+1:max.iter)
  
  x <- runif(1,-2,2)
  hval <- h(x)
  hcur <- hval
  rho_vec <- numeric()
  props <- x
  
  iter <- 1
  diff <- 1
  
  for (i in 1:max.iter) {
    prop <- x[i]+runif(1,-5,5)
    props <- c(props,prop)
    rho <- exp((h(prop)-hcur)/temp[i])
    rho <- min(rho,1)
    rho_vec <- c(rho_vec,rho)
    u <- runif(1)
    
    if((u>rho)||(prop>2)|| (prop< -2)) prop <- x[i]
    x <- c(x,prop)
    hcur <- h(prop)
    hval <- c(hval,hcur)
    
    diff <- max(hval)-max(hval[1:(i/2)])
    if((iter>min.iter)&&(length(unique(x[(i/2):i]))>1)&&(diff<tol)) break
    
    iter <- iter + 1
  }
  
  max_sa_x <- x[which.max(hval)]
  max_sa <- max(hval)
  return(list("Max x"=max_sa_x,"Max Value"=max_sa,"x"=x,"hval"=hval))
}

sa_result <- Simulated_annealing(1000,3000)
x = seq(-2,2,0.001)
ex = exp(-x^2)*sin(pi*x^3)
plot(y = ex, x = x, xlim = c(-2, 2),col="blue",cex=0.2,
     main="Simulated Annealing")
abline(h=0,col="gray2")
abline(v=0,col="gray2")
points(sa_result$x[1],sa_result$hval[1],col="red",pch=19)
lines(sa_result$x,sa_result$hval,col="red",type = "b")
points(x = sa_result$`Max x`,y = sa_result$`Max Value`,pch=19,col="aquamarine")
text(x = sa_result$`Max x`, y = sa_result$`Max Value`, labels = paste0("(",round(sa_result$`Max x`,4),",",round(sa_result$`Max Value`,4),")"), pos = 3, offset = 0.5)


# Comparison: Iter Lim = 100 | n=100 | miniter = 100 & maxiter = 300 

n <- 100

nelm <- nlm_simulation(n)
nlm_max <- nelm$`Max Nlm`
nlm_max_x <- nelm$`Max Nlm X`

bmco_u <- bmco_unif(n)
bmco_u_max <- bmco_u$`Maximum Unif Dist`
bmco_u_max_x <- bmco_u$`Maximum X-value Unif Dist`

bmco_b <- bmco_beta(n)
bmco_b_max <- bmco_b$`Maximum Beta Dist`
bmco_b_max_x <- bmco_b$`Maximum X-value Beta Dist`

sa_result <- Simulated_annealing(n,300)
sa_max <- sa_result$`Max Value`
sa_max_x <- sa_result$`Max x`

x = seq(-2,2,0.001)
ex = exp(-x^2)*sin(pi*x^3)
plot(y = ex, x = x, xlim = c(-2, 2),col="blue",cex=0.2,
     main="Cool Looking Function")
abline(h=0,col="gray2")
abline(v=0,col="gray2")
points(x = nlm_max_x, y = nlm_max, col = "purple", bg = "purple", cex = 1, type = "p", pch = 19)
points(x = bmco_u_max_x, y = bmco_u_max, col = "orange", bg = "orange", cex = 1, type = "p", pch = 19)
points(x = bmco_b_max_x, y = bmco_b_max, col = "green", bg = "green", cex = 1, type = "p", pch = 19)
points(x = max_sa_x, y = max_sa, col = "aquamarine", bg = "aquamarine", cex = 1, type = "p", pch = 19)
text(x = nlm_max_x, y = nlm_max, labels = "nlm_max", pos = 3, offset = 0.5)
text(x = bmco_u_max_x, y = bmco_u_max, labels = "bmco_u_max", pos = 3, offset = 0.5)
text(x = bmco_b_max_x, y = bmco_b_max, labels = "bmco_b_max", pos = 3, offset = 0.5)
text(x = sa_max_x, y = sa_max, labels = "max_sa", pos = 3, offset = 0.5)

# Table Format 

max_df <- data.frame("X" = round(c(nlm_max_x,bmco_u_max_x,bmco_b_max_x,sa_max_x),7),
                     "Maximum Result" = round(c(nlm_max,bmco_u_max,bmco_b_max,sa_max),8))

rownames(max_df) <- c("Newton Raphson","Basic-Uniform","Basic-Beta","Simulated Annealing")

max_df

# Simulation | n = 100

n <- 100
B <- 1000
res_nlm <- rep(0,B)
res_bmco_u <- rep(0,B)
res_bmco_b <- rep(0,B)
res_sa <- rep(0,B)

for (i in 1:B) {
  res_nlm[i] <- nlm_simulation(n)$`Max Nlm`
  res_bmco_u[i] <- bmco_unif(n)$`Maximum Unif Dist`
  res_bmco_b[i] <- bmco_beta(n)$`Maximum Beta Dist`
  res_sa[i] <- Simulated_annealing(1000,3000)$`Max Value`
}

df <- data.frame("1"=res_nlm,"2"=res_bmco_u,"3"=res_bmco_b,"4"=res_sa)
colnames(df) <- c("Newton-Raphson","Basic-MC-Uniform",
                  "Basic-MC-Beta","Simulated Annealing")
boxplot(df,range=0,col=c("orange","yellow","blue","lightblue")
        ,main="Simulation For Optimization Methods")


