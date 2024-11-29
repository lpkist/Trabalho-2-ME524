library(tidyverse)
newton_raphson <- function(y, x, beta_init, f, df_dbeta, d2f_dbeta2, tol,
                           tipo_line_search = 'nr', somar_diag_hess, 
                           gamma_init = 1, Q, exp_max, exp_min){
  av_grad <- 0
  av_hess <- 0
  av_Q <- 0
  n_it <- 0
  prod_vetores <- 0
  prod_matriz_vetor <- 0; gamma_init_aux <- gamma_init; dif <- Inf
  max_gamma <- -Inf
  min_gamma <- Inf
  calcula_grad_hess_Q <- function(beta){
    grad_Q <- matrix(0, nrow = length(beta))
    hess_Q <- matrix(0, nrow = length(beta), ncol = length(beta))
    for(i in 1:length(x)){
      x_i <- x[i]
      y_i <- y[i]
      f_x_beta <- f(x_i, beta)
      df_x_beta <- df_dbeta(x_i, beta)
      d2f_x_beta <- d2f_dbeta2(x_i, beta)
      grad_Q <- grad_Q + -2*(t((y_i - f_x_beta)*df_x_beta))
      hess_Q <- hess_Q + 2*((y_i - f_x_beta)*d2f_x_beta+
                              t(df_x_beta) %*% df_x_beta)
      h <<- hess_Q
    }
    return(list(grad_Q, hess_Q))
  }
  
  while(dif > tol){
    ders_Q <- calcula_grad_hess_Q(beta_init)
    grad_Q <- ders_Q[[1]]
    av_grad <- av_grad + 1
    hess_Q <- ders_Q[[2]]
    av_hess <- av_hess + 1
    if(tipo_line_search == 'nr'){
      Q_init <- Q(beta_init)
      av_Q <- av_Q + 1
      diag(hess_Q) <- diag(hess_Q) + somar_diag_hess
      p_k_1 <- solve(hess_Q, grad_Q)
      prod_matriz_vetor <- prod_matriz_vetor + 1
      q_primes <- function(gamma, p_k_1){
        beta_k <- beta_init - gamma*p_k_1
        ders_Q <- calcula_grad_hess_Q(beta_k)
        grad_Q <- ders_Q[[1]]
        hess_Q <- ders_Q[[2]]
        diag(hess_Q) <- diag(hess_Q) + somar_diag_hess
        q_prime <- -t(grad_Q) %*% p_k_1
        q_prime2 <- t(p_k_1) %*% hess_Q %*% p_k_1
        return(list(q_prime, q_prime2))
      }
      gamma_init <- gamma_init_aux
      dif_gamma <- Inf
      max_it <- 100
      n_it_gamma <- 0
      while((dif_gamma > 1e-4)&(n_it_gamma < max_it)){
        q_primes_i <- q_primes(gamma_init, p_k_1)
        av_grad <- av_grad + 1
        av_hess <- av_hess + 1
        prod_matriz_vetor <- prod_matriz_vetor + 1
        prod_vetores <- prod_vetores + 2
        q_prime <- q_primes_i[[1]]
        q_prime2 <- q_primes_i[[2]]
        gamma <- gamma_init - gamma_init_aux*q_prime/q_prime2
        gamma <- gamma[[1]]
        dif_gamma <- abs(gamma - gamma_init)
        gamma_init <- gamma
        n_it_gamma <- n_it_gamma+1
      }
      av_Q <- av_Q + 1
      if(Q(beta_init - gamma*p_k_1) > Q_init) gamma <- gamma_init_aux^0.1
    } else if(tipo_line_search == F){
      gamma <- gamma_init
      diag(hess_Q) <- diag(hess_Q) + somar_diag_hess
    }else if(tipo_line_search == 'Grid'){
      gamma <- 0
      if(!exists("exp_s")) exp_s <- NA
      Q_init <- Q(beta_init)
      av_Q <- av_Q + 1
      df <- data.frame(gamma = gamma, exp_s = exp_s, q = Q_init)
      for(exp_s in -12:20){
        hess_aux <- hess_Q
        diag(hess_aux) <- diag(hess_aux) + 3.1^exp_s
        try({
          p_k_1 <- solve(hess_aux, grad_Q)
          prod_matriz_vetor <- prod_matriz_vetor + 1
          arroz <- "feijao"
        })
        if(!exists("arroz")){
          next
        }
        rm(arroz)
        for(expoente in exp_max:exp_min){
          gamma <- 2^expoente
          df <- rbind(df, data.frame(gamma = gamma, exp_s = exp_s, q = Q(beta_init - gamma*p_k_1)))
          av_Q <- av_Q + 1
        }
      }
      df <- na.omit(df)
      gamma <- as.numeric(df[which.min(df$q), "gamma"])
      exp_s <- as.numeric(df[which.min(df$q), "exp_s"])
      hess_aux <- hess_Q
      diag(hess_aux) <- diag(hess_aux) + 3.1^exp_s
      hess_Q <- hess_aux
    }
    beta <- beta_init - gamma*solve(hess_Q, grad_Q)
    prod_matriz_vetor <- prod_matriz_vetor + 1
    dif <- sqrt(sum((beta-beta_init)^2)) # norma 2
    beta_init <- beta
    n_it <- n_it + 1
    ultimo_beta <<- beta
  }
  return(list(beta = beta, dif = dif, n_it = n_it, grad = grad_Q,
              hess = hess_Q, av_hess = av_hess,
              prod_matriz_vetor = prod_matriz_vetor,
              av_Q = av_Q, av_grad = av_grad, prod_vetores = prod_vetores))
}

f1 <- function(x, beta){
  return(beta[1]*(1-exp(-beta[2]*x)))
}
df1_dbeta <- function(x, beta){
  return(matrix(c(1-exp(-beta[2]*x), beta[1]*x*exp(-beta[2]*x)), ncol = 2))
}
d2f1_dbeta2 <- function(x, beta){
  return(matrix(c(0, x*exp(-beta[2]*x),
                  x*exp(-beta[2]*x), -beta[1]*x^2*exp(-beta[2]*x)), nrow = 2))
}

Q1_beta <- function(beta){
  Q <- 0
  for(i in 1:nrow(Misra1a)){
    Q <- Q + (Misra1a$y[i] - f1(Misra1a$x[i], beta))^2
  }
  Q
}

Misra1a <- data.frame(y = c(10.07, 14.73, 17.94, 23.93, 29.61,
                            35.18, 40.02, 44.82, 50.76, 55.05,
                            61.01, 66.4, 75.47, 81.78),
                      x = c(77.6, 114.9, 141.1, 190.8, 239.9,
                            289, 332.8, 378.4, 434.8, 477.3,
                            536.8, 593.1, 689.1, 760))


# ponto inicial 1
res_nr_1_1_Grid <- data.frame(tol = NA, dif_otimo = NA, dif = NA, 
                              norma2_grad = NA, av_Q = NA, n_it = NA,
                              av_grad = NA, tempo = NA, dif_q_ot = NA, 
                              av_hess = NA,
                              ls = "Grid")
beta_init <- c(500, 0.0001)
gabarito <- c(2.3894212918E+02, 5.5015643181E-04)
av_Q <- 0
av_grad <- 0
n_it <- 0
tempo <- 0
av_hess <- 0
for(tol in 2:-4){
  print(tol)
  tempo_i <- as.numeric(system.time({
    aux <- newton_raphson(y = Misra1a$y,
                          x = Misra1a$x,
                          beta_init = beta_init,
                          f =  f1,
                          df_dbeta = df1_dbeta,
                          d2f_dbeta2 = d2f1_dbeta2,
                          tol = 10^tol,
                          tipo_line_search = "Grid",
                          Q = Q1_beta, exp_max = 5, exp_min = -1)
  })['elapsed'])
  beta_init <- aux$beta
  av_Q <- av_Q + aux$av_Q
  av_grad <- av_grad + aux$av_grad
  av_hess <- av_hess + aux$av_hess
  n_it <- n_it + aux$n_it
  tempo <- tempo + tempo_i
  if(tol == -1){ nr11_ot <- aux; min(eigen(nr11_ot$hess)$values)}
  res_nr_1_1_Grid <- rbind(res_nr_1_1_Grid,
                           data.frame(tol = tol,
                                      dif_otimo = sqrt(sum((gabarito - beta_init)^2)), 
                                      dif = aux$dif, 
                                      norma2_grad = sqrt(sum(aux$grad^2)),
                                      av_Q = av_Q, 
                                      n_it = n_it,
                                      av_grad = av_grad,
                                      av_hess = av_hess,
                                      tempo = tempo,
                                      dif_q_ot = Q1_beta(beta_init) - Q1_beta(gabarito),
                                      ls = "Grid"))
}

res_nr_1_1_NR <- data.frame(tol = NA, dif_otimo = NA, dif = NA, 
                            norma2_grad = NA, av_Q = NA, n_it = NA,
                            av_grad = NA, tempo = NA, dif_q_ot = NA, 
                            av_hess = NA,
                            ls = "NR")
beta_init <- c(500, 0.0001)
gabarito <- c(2.3894212918E+02, 5.5015643181E-04)
av_Q <- 0
av_grad <- 0
n_it <- 0
tempo <- 0
av_hess <- 0
for(tol in 2:-4){
  print(tol)
  tempo_i <- as.numeric(system.time({
    aux <- newton_raphson(y = Misra1a$y,
                          x = Misra1a$x,
                          beta_init = beta_init,
                          f =  f1,
                          df_dbeta = df1_dbeta,
                          d2f_dbeta2 = d2f1_dbeta2,
                          tol = 10^tol,
                          tipo_line_search = "nr",
                          gamma_init = 1,
                          somar_diag_hess = 10^-2,
                          Q = Q1_beta, exp_max = 5, exp_min = -1)
  })['elapsed'])
  beta_init <- aux$beta
  av_Q <- av_Q + aux$av_Q
  av_grad <- av_grad + aux$av_grad
  av_hess <- av_hess + aux$av_hess
  n_it <- n_it + aux$n_it
  tempo <- tempo + tempo_i
  
  res_nr_1_1_NR <- rbind(res_nr_1_1_NR,
                         data.frame(tol = tol,
                                    dif_otimo = sqrt(sum((gabarito - beta_init)^2)), 
                                    dif = aux$dif, 
                                    norma2_grad = sqrt(sum(aux$grad^2)),
                                    av_Q = av_Q, 
                                    n_it = n_it,
                                    av_grad = av_grad,
                                    av_hess = av_hess,
                                    tempo = tempo,
                                    dif_q_ot = Q1_beta(beta_init) - Q1_beta(gabarito),
                                    ls = "NR"))
}

res_nr_1_1_Sem <- data.frame(tol = NA, dif_otimo = NA, dif = NA, 
                             norma2_grad = NA, av_Q = NA, n_it = NA,
                             av_grad = NA, tempo = NA, dif_q_ot = NA, 
                             av_hess = NA,
                             ls = "Sem")
beta_init <- c(500, 0.0001)
gabarito <- c(2.3894212918E+02, 5.5015643181E-04)
av_Q <- 0
av_grad <- 0
n_it <- 0
tempo <- 0
av_hess <- 0
for(tol in 2:-4){
  print(tol)
  tempo_i <- as.numeric(system.time({
    aux <- newton_raphson(y = Misra1a$y,
                          x = Misra1a$x,
                          beta_init = beta_init,
                          f =  f1,
                          df_dbeta = df1_dbeta,
                          d2f_dbeta2 = d2f1_dbeta2,
                          tol = 10^tol,
                          tipo_line_search = F,
                          gamma_init = 1,
                          somar_diag_hess = 0.1,
                          Q = Q1_beta, exp_max = 5, exp_min = -1)
  })['elapsed'])
  beta_init <- aux$beta
  av_Q <- av_Q + aux$av_Q
  av_grad <- av_grad + aux$av_grad
  av_hess <- av_hess + aux$av_hess
  n_it <- n_it + aux$n_it
  tempo <- tempo + tempo_i
  
  res_nr_1_1_Sem <- rbind(res_nr_1_1_Sem,
                          data.frame(tol = tol,
                                     dif_otimo = sqrt(sum((gabarito - beta_init)^2)), 
                                     dif = aux$dif, 
                                     norma2_grad = sqrt(sum(aux$grad^2)),
                                     av_Q = av_Q, 
                                     n_it = n_it,
                                     av_grad = av_grad,
                                     av_hess = av_hess,
                                     tempo = tempo,
                                     dif_q_ot = Q1_beta(beta_init) - Q1_beta(gabarito),
                                     ls = "Sem"))
}



nr11 <- 
  rbind(res_nr_1_1_Grid,
        res_nr_1_1_NR,
        res_nr_1_1_Sem) %>% 
  transmute("Tolerância" = tol,
            "Log10 do nº de avaliações do gradiente" = log(av_grad, 10),
            "Log10 do nº de iterações" = log(n_it, 10),
            "Log10 do nº de avaliações de Q1 + 1" =log(av_Q + 1, 10),
            "Log10 do nº de avaliações de H1 + 1" = log(av_hess + 1, 10),
            "Log10 da distância ao mínimo global" = log(dif_otimo, 10),
            "Log10 da norma 2 do gradiente" = log(norma2_grad, 10),
            "Log10 da diferença para a perda ótima" = log(ifelse(dif_q_ot < 0, 10^-15, dif_q_ot), 10),
            "Tempo (s) até a convergência" = tempo,
            "Método de line-search" = ls,
  ) %>% 
  pivot_longer(cols = 2:9, names_to = "kpi", values_to = "Valor") %>% 
  ggplot(aes(x = `Tolerância`, y = Valor, color = `Método de line-search`))+
  geom_point()+
  geom_line()+
  facet_wrap(~kpi, scales = "free")+
  scale_x_reverse(breaks = -4:2, labels = paste0("10^", -4:2))+
  theme_bw()+
  ggtitle('Ponto inicial: (500, 0.0001)')+
  guides(color = guide_legend(position = "top", 
                              title = "Método de line-search",
                              title.position="top"))+
  theme(plot.title = element_text(size=18, hjust = 0.5),
        axis.title = element_text(size=14),
        legend.title = element_text(size=14))

# ponto inicial 2
res_nr_1_2_Grid <- data.frame(tol = NA, dif_otimo = NA, dif = NA, 
                              norma2_grad = NA, av_Q = NA, n_it = NA,
                              av_grad = NA, tempo = NA, dif_q_ot = NA, 
                              av_hess = NA,
                              ls = "Grid")
beta_init <- c(250,0.0005)
gabarito <- c(2.3894212918E+02, 5.5015643181E-04)
av_Q <- 0
av_grad <- 0
n_it <- 0
tempo <- 0
av_hess <- 0
for(tol in 2:-4){
  print(tol)
  tempo_i <- as.numeric(system.time({
    aux <- newton_raphson(y = Misra1a$y,
                          x = Misra1a$x,
                          beta_init = beta_init,
                          f =  f1,
                          df_dbeta = df1_dbeta,
                          d2f_dbeta2 = d2f1_dbeta2,
                          tol = 10^tol,
                          tipo_line_search = "Grid",
                          Q = Q1_beta, exp_max = 5, exp_min = -1)
  })['elapsed'])
  beta_init <- aux$beta
  av_Q <- av_Q + aux$av_Q
  av_grad <- av_grad + aux$av_grad
  av_hess <- av_hess + aux$av_hess
  n_it <- n_it + aux$n_it
  tempo <- tempo + tempo_i
  if(tol == -1){ nr12_ot <- aux; min(eigen(nr12_ot$hess)$values)}
  res_nr_1_2_Grid <- rbind(res_nr_1_2_Grid,
                           data.frame(tol = tol,
                                      dif_otimo = sqrt(sum((gabarito - beta_init)^2)), 
                                      dif = aux$dif, 
                                      norma2_grad = sqrt(sum(aux$grad^2)),
                                      av_Q = av_Q, 
                                      n_it = n_it,
                                      av_grad = av_grad,
                                      av_hess = av_hess,
                                      tempo = tempo,
                                      dif_q_ot = Q1_beta(beta_init) - Q1_beta(gabarito),
                                      ls = "Grid"))
}

res_nr_1_2_NR <- data.frame(tol = NA, dif_otimo = NA, dif = NA, 
                            norma2_grad = NA, av_Q = NA, n_it = NA,
                            av_grad = NA, tempo = NA, dif_q_ot = NA, 
                            av_hess = NA,
                            ls = "NR")
beta_init <- c(250,0.0005)
gabarito <- c(2.3894212918E+02, 5.5015643181E-04)
av_Q <- 0
av_grad <- 0
n_it <- 0
tempo <- 0
av_hess <- 0
for(tol in 2:-4){
  print(tol)
  tempo_i <- as.numeric(system.time({
    aux <- newton_raphson(y = Misra1a$y,
                          x = Misra1a$x,
                          beta_init = beta_init,
                          f =  f1,
                          df_dbeta = df1_dbeta,
                          d2f_dbeta2 = d2f1_dbeta2,
                          tol = 10^tol,
                          tipo_line_search = "nr",
                          gamma_init = 1,
                          somar_diag_hess = 10^-2,
                          Q = Q1_beta, exp_max = 5, exp_min = -1)
  })['elapsed'])
  beta_init <- aux$beta
  av_Q <- av_Q + aux$av_Q
  av_grad <- av_grad + aux$av_grad
  av_hess <- av_hess + aux$av_hess
  n_it <- n_it + aux$n_it
  tempo <- tempo + tempo_i
  
  res_nr_1_2_NR <- rbind(res_nr_1_2_NR,
                         data.frame(tol = tol,
                                    dif_otimo = sqrt(sum((gabarito - beta_init)^2)), 
                                    dif = aux$dif, 
                                    norma2_grad = sqrt(sum(aux$grad^2)),
                                    av_Q = av_Q, 
                                    n_it = n_it,
                                    av_grad = av_grad,
                                    av_hess = av_hess,
                                    tempo = tempo,
                                    dif_q_ot = Q1_beta(beta_init) - Q1_beta(gabarito),
                                    ls = "NR"))
}

res_nr_1_2_Sem <- data.frame(tol = NA, dif_otimo = NA, dif = NA, 
                             norma2_grad = NA, av_Q = NA, n_it = NA,
                             av_grad = NA, tempo = NA, dif_q_ot = NA, 
                             av_hess = NA,
                             ls = "Sem")
beta_init <- c(250,0.0005)
gabarito <- c(2.3894212918E+02, 5.5015643181E-04)
av_Q <- 0
av_grad <- 0
n_it <- 0
tempo <- 0
av_hess <- 0
for(tol in 2:-4){
  print(tol)
  tempo_i <- as.numeric(system.time({
    aux <- newton_raphson(y = Misra1a$y,
                          x = Misra1a$x,
                          beta_init = beta_init,
                          f =  f1,
                          df_dbeta = df1_dbeta,
                          d2f_dbeta2 = d2f1_dbeta2,
                          tol = 10^tol,
                          tipo_line_search = F,
                          gamma_init = 1,
                          somar_diag_hess = 0.01,
                          Q = Q1_beta, exp_max = 5, exp_min = -1)
  })['elapsed'])
  beta_init <- aux$beta
  av_Q <- av_Q + aux$av_Q
  av_grad <- av_grad + aux$av_grad
  av_hess <- av_hess + aux$av_hess
  n_it <- n_it + aux$n_it
  tempo <- tempo + tempo_i
  
  res_nr_1_2_Sem <- rbind(res_nr_1_2_Sem,
                          data.frame(tol = tol,
                                     dif_otimo = sqrt(sum((gabarito - beta_init)^2)), 
                                     dif = aux$dif, 
                                     norma2_grad = sqrt(sum(aux$grad^2)),
                                     av_Q = av_Q, 
                                     n_it = n_it,
                                     av_grad = av_grad,
                                     av_hess = av_hess,
                                     tempo = tempo,
                                     dif_q_ot = Q1_beta(beta_init) - Q1_beta(gabarito),
                                     ls = "Sem"))
}



nr12 <- rbind(res_nr_1_2_Grid,
              res_nr_1_2_NR,
              res_nr_1_2_Sem) %>% 
  transmute("Tolerância" = tol,
            "Log10 do nº de avaliações do gradiente" = log(av_grad, 10),
            "Log10 do nº de iterações" = log(n_it, 10),
            "Log10 do nº de avaliações de Q2 + 1" =log(av_Q + 1, 10),
            "Log10 do nº de avaliações de H2 + 1" = log(av_hess + 1, 10),
            "Log10 da distância ao mínimo global" = log(dif_otimo, 10),
            "Log10 da norma 2 do gradiente" = log(norma2_grad, 10),
            "Log10 da diferença para a perda ótima" = log(ifelse(dif_q_ot < 0, 10^-15, dif_q_ot), 10),
            "Tempo (s) até a convergência" = tempo,
            "Método de line-search" = ls,
  ) %>% 
  pivot_longer(cols = 2:9, names_to = "kpi", values_to = "Valor") %>% 
  ggplot(aes(x = `Tolerância`, y = Valor, color = `Método de line-search`))+
  geom_point()+
  geom_line()+
  facet_wrap(~kpi, scales = "free")+
  scale_x_reverse(breaks = -4:2, labels = paste0("10^", -4:2))+
  theme_bw()+
  ggtitle('Ponto inicial: (250,0.0005)')+
  guides(color = guide_legend(position = "top", 
                              title = "Método de line-search",
                              title.position="top"))+
  theme(plot.title = element_text(size=18, hjust = 0.5),
        axis.title = element_text(size=14),
        legend.title = element_text(size=14))

# problema 2
Q2_beta <- function(beta){
  Q <- 0
  for(i in 1:nrow(Thurber)){
    Q <- Q + (Thurber$y[i] - f2(Thurber$x[i], beta))^2
  }
  Q
}

f2 <- function(x, beta){
  betas_1 <- matrix(beta[1:4], nrow = 4)
  X1 <- matrix(c(1, x, x^2, x^3), ncol = 4)
  
  betas_2 <- matrix(beta[5:7], nrow = 3)
  X2 <-  matrix(c(x, x^2, x^3), ncol = 3)
  f2 <- as.numeric((1 + X2 %*% betas_2)^(-1)*X1 %*% betas_1)
  if(is.finite(f2)){
    return(f2)
  } else{
    return(as.numeric((0.0001)^(-1)*X1 %*% betas_1))
  }
}
df2_dbeta <- function(x, beta){
  betas_1 <- matrix(beta[1:4], nrow = 4)
  X1 <- matrix(c(1, x, x^2, x^3), ncol = 4)
  
  betas_2 <- matrix(beta[5:7], nrow = 3)
  X2 <-  matrix(c(x, x^2, x^3), ncol = 3)
  
  der <- matrix(c(as.numeric((1 + X2 %*% betas_2)^(-1))*X1, 
                  -as.numeric((1 + X2 %*% betas_2))^(-2)*X1 %*% betas_1 %*% X2), ncol = 7)
  if(sum(!is.finite(der))){
    return(matrix(c(as.numeric((0.0001)^(-1))*X1, 
                    -as.numeric((0.0001))^(-2)*X1 %*% betas_1 %*% X2), ncol = 7))
  } else{
    return(der)
  }
  
}

d2f2_dbeta2 <- function(x, beta){
  betas_1 <- matrix(beta[1:4], nrow = 4)
  X1 <- matrix(c(1, x, x^2, x^3), ncol = 4)
  
  betas_2 <- matrix(beta[5:7], nrow = 3)
  X2 <-  matrix(c(x, x^2, x^3), ncol = 3)
  
  hess <- matrix(0, nrow = 7, ncol = 7)
  
  hess[5:7, 1:4] <- as.numeric(-(1 + X2 %*% betas_2)^(-2))*t(X2) %*% X1 
  if(sum(!is.finite(hess[5:7, 1:4]))){
    hess[5:7, 1:4] <- as.numeric((0.0001)^(-2))*t(X2) %*% X1 
  }
  hess[1:4, 5:7] <- as.numeric(-(1 + X2 %*% betas_2)^(-2))*t(X1) %*% X2
  if(sum(!is.finite(hess[1:4, 5:7]))){
    hess[1:4, 5:7] <- as.numeric((0.0001)^(-2))*t(X1) %*% X2
  }
  hess[5:7, 5:7] <- as.numeric(2*(1 + X2 %*% betas_2)^(-3))*t(X2) %*% t(betas_1) %*% t(X1) %*% X2
  if(sum(!is.finite(hess[5:7, 5:7]))){
    hess[5:7, 5:7] <- as.numeric(-2*(0.0001)^(-3))*t(X2) %*% t(betas_1) %*% t(X1) %*% X2
  }  
  return(hess)
}

Thurber <- data.frame(y = c(80.574, 84.248, 87.264, 87.195, 89.076, 89.608,
                            89.868, 90.101, 92.405, 95.854, 100.696, 101.06,
                            401.672, 390.724, 567.534, 635.316, 733.054,
                            759.087, 894.206, 990.785, 1090.109, 1080.914,
                            1122.643, 1178.351, 1260.531, 1273.514, 1288.339,
                            1327.543, 1353.863, 1414.509, 1425.208, 1421.384,
                            1442.962, 1464.35, 1468.705, 1447.894, 1457.628),
                      x = c(-3.067,-2.981,-2.921,-2.912,-2.84,-2.797,-2.702,-2.699,-2.633,-2.481,-2.363,-2.322,-1.501,-1.46,-1.274,-1.212,-1.1,-1.046,-0.915,-0.714,-0.566,-0.545,-0.4,-0.309,-0.109,-0.103, 0.01, 0.119, 0.377, 0.79, 0.963,
                            1.006, 1.115, 1.572, 1.841, 2.047, 2.2))


# ponto inicial 1
res_nr_2_1_Grid <- data.frame(tol = NA, dif_otimo = NA, dif = NA, 
                              norma2_grad = NA, av_Q = NA, n_it = NA,
                              av_grad = NA, tempo = NA, dif_q_ot = NA, 
                              av_hess = NA,
                              ls = "Grid")
beta_init <- c(1300,1500,500,75,1,0.4,0.05)
gabarito <- c(1.2881396800E+03,
              1.4910792535E+03,
              5.8323836877E+02,
              7.5416644291E+01,
              9.6629502864E-01,
              3.9797285797E-01,
              4.9727297349E-02)
av_Q <- 0
av_grad <- 0
n_it <- 0
tempo <- 0
av_hess <- 0
for(tol in 2:-4){
  print(tol)
  tempo_i <- as.numeric(system.time({
    aux <- newton_raphson(y = Thurber$y,
                          x = Thurber$x,
                          beta_init = beta_init,
                          f =  f2,
                          df_dbeta = df2_dbeta,
                          d2f_dbeta2 = d2f2_dbeta2,
                          tol = 10^tol,
                          tipo_line_search = "Grid",
                          Q = Q2_beta, exp_max = 20, exp_min = -10)
  })['elapsed'])
  beta_init <- aux$beta
  av_Q <- av_Q + aux$av_Q
  av_grad <- av_grad + aux$av_grad
  av_hess <- av_hess + aux$av_hess
  n_it <- n_it + aux$n_it
  tempo <- tempo + tempo_i
  res_nr_2_1_Grid <- rbind(res_nr_2_1_Grid,
                           data.frame(tol = tol,
                                      dif_otimo = sqrt(sum((gabarito - beta_init)^2)), 
                                      dif = aux$dif, 
                                      norma2_grad = sqrt(sum(aux$grad^2)),
                                      av_Q = av_Q, 
                                      n_it = n_it,
                                      av_grad = av_grad,
                                      av_hess = av_hess,
                                      tempo = tempo,
                                      dif_q_ot = Q2_beta(beta_init) - Q2_beta(gabarito),
                                      ls = "Grid"))
}
nr21_ot <- aux; min(eigen(nr21_ot$hess)$values)
res_nr_2_1_NR <- data.frame(tol = NA, dif_otimo = NA, dif = NA, 
                            norma2_grad = NA, av_Q = NA, n_it = NA,
                            av_grad = NA, tempo = NA, dif_q_ot = NA, 
                            av_hess = NA,
                            ls = "NR")
beta_init <- c(1300,1500,500,75,1,0.4,0.05)
gabarito <- c(1.2881396800E+03,
              1.4910792535E+03,
              5.8323836877E+02,
              7.5416644291E+01,
              9.6629502864E-01,
              3.9797285797E-01,
              4.9727297349E-02)
av_Q <- 0
av_grad <- 0
n_it <- 0
tempo <- 0
av_hess <- 0
for(tol in 2:-3){
  print(tol)
  tempo_i <- as.numeric(system.time({
    aux <- newton_raphson(y = Thurber$y,
                          x = Thurber$x,
                          beta_init = beta_init,
                          f =  f2,
                          df_dbeta = df2_dbeta,
                          d2f_dbeta2 = d2f2_dbeta2,
                          tol = 10^tol,
                          tipo_line_search = "nr",
                          gamma_init = 1e-5,
                          somar_diag_hess = 10^4,
                          Q = Q2_beta, exp_max = 5, exp_min = -1)
  })['elapsed'])
  beta_init <- aux$beta
  av_Q <- av_Q + aux$av_Q
  av_grad <- av_grad + aux$av_grad
  av_hess <- av_hess + aux$av_hess
  n_it <- n_it + aux$n_it
  tempo <- tempo + tempo_i
  
  res_nr_2_1_NR <- rbind(res_nr_2_1_NR,
                         data.frame(tol = tol,
                                    dif_otimo = sqrt(sum((gabarito - beta_init)^2)), 
                                    dif = aux$dif, 
                                    norma2_grad = sqrt(sum(aux$grad^2)),
                                    av_Q = av_Q, 
                                    n_it = n_it,
                                    av_grad = av_grad,
                                    av_hess = av_hess,
                                    tempo = tempo,
                                    dif_q_ot = Q2_beta(beta_init) - Q2_beta(gabarito),
                                    ls = "NR"))
}

res_nr_2_1_Sem <- data.frame(tol = NA, dif_otimo = NA, dif = NA, 
                             norma2_grad = NA, av_Q = NA, n_it = NA,
                             av_grad = NA, tempo = NA, dif_q_ot = NA, 
                             av_hess = NA,
                             ls = "Sem")
beta_init <- c(1300,1500,500,75,1,0.4,0.05)
gabarito <- c(1.2881396800E+03,
              1.4910792535E+03,
              5.8323836877E+02,
              7.5416644291E+01,
              9.6629502864E-01,
              3.9797285797E-01,
              4.9727297349E-02)
av_Q <- 0
av_grad <- 0
n_it <- 0
tempo <- 0
av_hess <- 0
for(tol in 2:-3){
  print(tol)
  tempo_i <- as.numeric(system.time({
    aux <- newton_raphson(y = Thurber$y,
                          x = Thurber$x,
                          beta_init = beta_init,
                          f =  f2,
                          df_dbeta = df2_dbeta,
                          d2f_dbeta2 = d2f2_dbeta2,
                          tol = 10^tol,
                          tipo_line_search = F,
                          gamma_init = 10^-3,
                          somar_diag_hess = 10^12,
                          Q = Q2_beta, exp_max = 5, exp_min = -1)
  })['elapsed'])
  beta_init <- aux$beta
  av_Q <- av_Q + aux$av_Q
  av_grad <- av_grad + aux$av_grad
  av_hess <- av_hess + aux$av_hess
  n_it <- n_it + aux$n_it
  tempo <- tempo + tempo_i
  
  res_nr_2_1_Sem <- rbind(res_nr_2_1_Sem,
                          data.frame(tol = tol,
                                     dif_otimo = sqrt(sum((gabarito - beta_init)^2)), 
                                     dif = aux$dif, 
                                     norma2_grad = sqrt(sum(aux$grad^2)),
                                     av_Q = av_Q, 
                                     n_it = n_it,
                                     av_grad = av_grad,
                                     av_hess = av_hess,
                                     tempo = tempo,
                                     dif_q_ot = Q2_beta(beta_init) - Q2_beta(gabarito),
                                     ls = "Sem"))
}

nr21 <- res_nr_2_1_Grid %>% 
  transmute("Tolerância" = tol,
            "Log10 do nº de avaliações do gradiente" = log(av_grad, 10),
            "Log10 do nº de iterações" = log(n_it, 10),
            "Log10 do nº de avaliações de Q2 + 1" =log(av_Q + 1, 10),
            "Log10 do nº de avaliações de H2 + 1" = log(av_hess + 1, 10),
            "Log10 da distância ao mínimo global" = log(dif_otimo, 10),
            "Log10 da norma 2 do gradiente" = log(norma2_grad, 10),
            "Log10 da diferença para a perda ótima" = log(ifelse(dif_q_ot < 0, 10^-15, dif_q_ot), 10),
            "Tempo (s) até a convergência" = tempo,
            "Método de line-search" = ls,
  ) %>% 
  pivot_longer(cols = 2:9, names_to = "kpi", values_to = "Valor") %>% 
  ggplot(aes(x = `Tolerância`, y = Valor, color = `Método de line-search`))+
  geom_point()+
  geom_line()+
  facet_wrap(~kpi, scales = "free")+
  scale_x_reverse(breaks = -4:2, labels = paste0("10^", -4:2))+
  theme_bw()+
  ggtitle('Ponto inicial: (1300,1500,500,75,1,0.4,0.05)')+
  guides(color = guide_legend(position = "top", 
                              title = "Método de line-search",
                              title.position="top"))+
  theme(plot.title = element_text(size=18, hjust = 0.5),
        axis.title = element_text(size=14),
        legend.title = element_text(size=14))


# ponto inicial 2
res_nr_2_2_Grid <- data.frame(tol = NA, dif_otimo = NA, dif = NA, 
                              norma2_grad = NA, av_Q = NA, n_it = NA,
                              av_grad = NA, tempo = NA, dif_q_ot = NA, 
                              av_hess = NA,
                              ls = "Grid")
beta_init <- c(1000, 1000, 400, 40, 0.7, 0.3, 0.03)
gabarito <- c(1.2881396800E+03,
              1.4910792535E+03,
              5.8323836877E+02,
              7.5416644291E+01,
              9.6629502864E-01,
              3.9797285797E-01,
              4.9727297349E-02)
av_Q <- 0
av_grad <- 0
n_it <- 0
tempo <- 0
av_hess <- 0
for(tol in 2:-4){
  print(tol)
  tempo_i <- as.numeric(system.time({
    aux <- newton_raphson(y = Thurber$y,
                          x = Thurber$x,
                          beta_init = beta_init,
                          f =  f2,
                          df_dbeta = df2_dbeta,
                          d2f_dbeta2 = d2f2_dbeta2,
                          tol = 10^tol,
                          tipo_line_search = "Grid",
                          Q = Q2_beta, exp_max = 20, exp_min = -80)
  })['elapsed'])
  beta_init <- aux$beta
  av_Q <- av_Q + aux$av_Q
  av_grad <- av_grad + aux$av_grad
  av_hess <- av_hess + aux$av_hess
  n_it <- n_it + aux$n_it
  tempo <- tempo + tempo_i
  
  res_nr_2_2_Grid <- rbind(res_nr_2_2_Grid,
                           data.frame(tol = tol,
                                      dif_otimo = sqrt(sum((gabarito - beta_init)^2)), 
                                      dif = aux$dif, 
                                      norma2_grad = sqrt(sum(aux$grad^2)),
                                      av_Q = av_Q, 
                                      n_it = n_it,
                                      av_grad = av_grad,
                                      av_hess = av_hess,
                                      tempo = tempo,
                                      dif_q_ot = Q2_beta(beta_init) - Q2_beta(gabarito),
                                      ls = "Grid"))
}
nr22_ot <- aux; min(eigen(nr22_ot$hess)$values)
res_nr_2_2_NR <- data.frame(tol = NA, dif_otimo = NA, dif = NA, 
                            norma2_grad = NA, av_Q = NA, n_it = NA,
                            av_grad = NA, tempo = NA, dif_q_ot = NA, 
                            av_hess = NA,
                            ls = "NR")
beta_init <- c(1000, 1000, 400, 40, 0.7, 0.3, 0.03)
gabarito <- c(1.2881396800E+03,
              1.4910792535E+03,
              5.8323836877E+02,
              7.5416644291E+01,
              9.6629502864E-01,
              3.9797285797E-01,
              4.9727297349E-02)
av_Q <- 0
av_grad <- 0
n_it <- 0
tempo <- 0
av_hess <- 0
for(tol in 2:-4){
  print(tol)
  tempo_i <- as.numeric(system.time({
    aux <- newton_raphson(y = Thurber$y,
                          x = Thurber$x,
                          beta_init = beta_init,
                          f =  f2,
                          df_dbeta = df2_dbeta,
                          d2f_dbeta2 = d2f2_dbeta2,
                          tol = 10^tol,
                          tipo_line_search = "nr",
                          gamma_init = 1e-3,
                          somar_diag_hess = 10^10,
                          Q = Q2_beta, exp_max = 5, exp_min = -1)
  })['elapsed'])
  beta_init <- aux$beta
  av_Q <- av_Q + aux$av_Q
  av_grad <- av_grad + aux$av_grad
  av_hess <- av_hess + aux$av_hess
  n_it <- n_it + aux$n_it
  tempo <- tempo + tempo_i
  
  res_nr_2_2_NR <- rbind(res_nr_2_2_NR,
                         data.frame(tol = tol,
                                    dif_otimo = sqrt(sum((gabarito - beta_init)^2)), 
                                    dif = aux$dif, 
                                    norma2_grad = sqrt(sum(aux$grad^2)),
                                    av_Q = av_Q, 
                                    n_it = n_it,
                                    av_grad = av_grad,
                                    av_hess = av_hess,
                                    tempo = tempo,
                                    dif_q_ot = Q2_beta(beta_init) - Q2_beta(gabarito),
                                    ls = "NR"))
}

res_nr_2_2_Sem <- data.frame(tol = NA, dif_otimo = NA, dif = NA, 
                             norma2_grad = NA, av_Q = NA, n_it = NA,
                             av_grad = NA, tempo = NA, dif_q_ot = NA, 
                             av_hess = NA,
                             ls = "Sem")
beta_init <- c(1000, 1000, 400, 40, 0.7, 0.3, 0.03)
gabarito <- c(1.2881396800E+03,
              1.4910792535E+03,
              5.8323836877E+02,
              7.5416644291E+01,
              9.6629502864E-01,
              3.9797285797E-01,
              4.9727297349E-02)
av_Q <- 0
av_grad <- 0
n_it <- 0
tempo <- 0
av_hess <- 0
for(tol in 2:-4){
  print(tol)
  tempo_i <- as.numeric(system.time({
    aux <- newton_raphson(y = Thurber$y,
                          x = Thurber$x,
                          beta_init = beta_init,
                          f =  f2,
                          df_dbeta = df2_dbeta,
                          d2f_dbeta2 = d2f2_dbeta2,
                          tol = 10^tol,
                          tipo_line_search = F,
                          gamma_init = 10^-6,
                          somar_diag_hess = 10^6,
                          Q = Q2_beta, exp_max = 5, exp_min = -1)
  })['elapsed'])
  beta_init <- aux$beta
  av_Q <- av_Q + aux$av_Q
  av_grad <- av_grad + aux$av_grad
  av_hess <- av_hess + aux$av_hess
  n_it <- n_it + aux$n_it
  tempo <- tempo + tempo_i
  
  res_nr_2_2_Sem <- rbind(res_nr_2_2_Sem,
                          data.frame(tol = tol,
                                     dif_otimo = sqrt(sum((gabarito - beta_init)^2)), 
                                     dif = aux$dif, 
                                     norma2_grad = sqrt(sum(aux$grad^2)),
                                     av_Q = av_Q, 
                                     n_it = n_it,
                                     av_grad = av_grad,
                                     av_hess = av_hess,
                                     tempo = tempo,
                                     dif_q_ot = Q2_beta(beta_init) - Q2_beta(gabarito),
                                     ls = "Sem"))
}



nr22 <- res_nr_2_2_Grid %>% 
  transmute("Tolerância" = tol,
            "Log10 do nº de avaliações do gradiente" = log(av_grad, 10),
            "Log10 do nº de iterações" = log(n_it, 10),
            "Log10 do nº de avaliações de Q2 + 1" =log(av_Q + 1, 10),
            "Log10 do nº de avaliações de H2 + 1" = log(av_hess + 1, 10),
            "Log10 da distância ao mínimo global" = log(dif_otimo, 10),
            "Log10 da norma 2 do gradiente" = log(norma2_grad, 10),
            "Log10 da diferença para a perda ótima" = log(ifelse(dif_q_ot < 0, 10^-15, dif_q_ot), 10),
            "Tempo (s) até a convergência" = tempo,
            "Método de line-search" = ls,
  ) %>% 
  pivot_longer(cols = 2:9, names_to = "kpi", values_to = "Valor") %>% 
  ggplot(aes(x = `Tolerância`, y = Valor, color = `Método de line-search`))+
  geom_point()+
  geom_line()+
  facet_wrap(~kpi, scales = "free")+
  scale_x_reverse(breaks = -4:2, labels = paste0("10^", -4:2))+
  theme_bw()+
  ggtitle('Ponto inicial: (1000, 1000, 400, 40, 0.7, 0.3, 0.03)')+
  guides(color = guide_legend(position = "top", 
                              title = "Método de line-search",
                              title.position="top"))+
  theme(plot.title = element_text(size=18, hjust = 0.5),
        axis.title = element_text(size=14),
        legend.title = element_text(size=14))

# fim do problema 2
# ggsave("nr1.pdf", nr11+nr12, width = 1900, height = 850, units = "px")
# ggsave("nr2.pdf", nr21+nr22, width = 1900, height = 850, units = "px")


nr11_ot <- newton_raphson(y = Misra1a$y,
                          x = Misra1a$x,
                          beta_init = c(500, 0.0001),
                          f =  f1,
                          df_dbeta = df1_dbeta,
                          d2f_dbeta2 = d2f1_dbeta2,
                          tol = 10^-1,
                          tipo_line_search = "Grid",
                          Q = Q1_beta, exp_max = 5, exp_min = -1)
nr11_ot$beta
min(eigen(nr11_ot$hess)$values)

library(tidyverse)
## grad desc
gradient_descent <-  function(y, x, beta_init, f, df_dbeta, tol,
                              tipo_line_search = 'gd', gamma_init = 1, Q,
                              exp_max, exp_min){
  av_grad <- 0; av_Q <- 0; n_it <- 0; prod_grad <- 0
  gamma_init_or <- gamma_init
  dif <- Inf; max_gamma <- -Inf; min_gamma <- Inf
  calcula_grad_Q <- function(beta){ # um cálculo de gradiente
    grad_Q <- matrix(0, nrow = length(beta))
    for(i in 1:length(x)){
      x_i <- x[i]
      y_i <- y[i]
      f_x_beta <- f(x_i, beta)
      df_x_beta <- df_dbeta(x_i, beta)
      grad_Q <- grad_Q -2*(t((y_i - f_x_beta)*df_x_beta))
    }
    return(grad_Q)
  }
  
  while(dif > tol){
    grad_Q <- calcula_grad_Q(beta_init)
    av_grad <- av_grad + 1
    p_k_1 <- grad_Q
    
    if(tipo_line_search == 'gd'){
      q_primes <- function(gamma, grad_Q, beta_init){ # um produto de gradientes
        beta_k <- beta_init - gamma*grad_Q
        grad_at <- calcula_grad_Q(beta_k)
        q_prime <- -t(grad_at) %*% grad_Q
        return(q_prime)
      }
      dif_gamma <- Inf
      max_it <- 100
      n_it_gamma <- 0
      gamma_init_aux <- gamma_init_or
      Q_init <- Q(beta_init)
      while((dif_gamma > min(0.1, tol))&(n_it_gamma < max_it)){
        # aqui é o while do ls
        q_prime <- q_primes(gamma_init, grad_Q, beta_init)
        av_grad <- av_grad + 1
        prod_grad <- prod_grad + 1
        
        gamma <- as.numeric(gamma_init - gamma_init_aux*q_prime)
        if(Q(beta_init - gamma*grad_Q) > Q_init){
          n_it_gamma <- n_it_gamma + 1
          gamma_init_aux <- gamma_init_aux/2
          if(n_it_gamma == max_it){ 
            if(exists("gamma_backup")) {gamma <- gamma_backup} else gamma <- gamma_init_or^0.1
          }
          next
        } else{
          gamma_backup <- gamma
          gamma_init_aux <- gamma_init_or
        }
        if(!is.finite(gamma)) gamma <- gamma_init_aux
        gamma <- gamma[[1]]
        dif_gamma <- abs(gamma - gamma_init)
        gamma_init <- gamma
        n_it_gamma <- n_it_gamma+1
      }
      if(exists("gamma_backup")) rm(gamma_backup)
    }else if(tipo_line_search == 'Grid'){
      avalia_expoente <- function(expoente){
        gamma <- 2^expoente
        data.frame(gamma = gamma, q = Q(beta_init - gamma*grad_Q))
      }
      
      df <- list_rbind(map(c(0, exp_max:exp_min), avalia_expoente))
      av_Q <- av_Q + length(exp_max:exp_min) + 1
      
      gamma <- as.numeric(df[which.min(df$q), "gamma"])
      
    } else if(tipo_line_search == F){
      gamma <- gamma_init
    } 
    if(gamma < min_gamma) min_gamma <- gamma
    if(gamma > max_gamma) max_gamma <- gamma
    
    beta <- beta_init - gamma*grad_Q
    #n_op <- n_op + p*(2*p-1)
    dif <- sqrt(sum((beta-beta_init)^2)) # norma 2
    beta_init <- beta
    n_it <- n_it + 1
    
  }
  return(list(beta = beta, dif = dif, n_it = n_it, grad = grad_Q,
              av_Q = av_Q, av_grad = av_grad, prod_grad = prod_grad))
}

# começo do problema 1
f1 <- function(x, beta){
  return(beta[1]*(1-exp(-beta[2]*x)))
}
df1_dbeta <- function(x, beta){
  return(matrix(c(1-exp(-beta[2]*x), beta[1]*x*exp(-beta[2]*x)), ncol = 2))
}

Misra1a <- data.frame(y = c(10.07, 14.73, 17.94, 23.93, 29.61,
                            35.18, 40.02, 44.82, 50.76, 55.05,
                            61.01, 66.4, 75.47, 81.78),
                      x = c(77.6, 114.9, 141.1, 190.8, 239.9,
                            289, 332.8, 378.4, 434.8, 477.3,
                            536.8, 593.1, 689.1, 760))

Q1_beta <- function(beta){
  Q <- 0
  for(i in 1:nrow(Misra1a)){
    Q <- Q + (Misra1a$y[i] - f1(Misra1a$x[i]/1e6, beta))^2
  }
  Q
}

# ponto inicial 1
res_gd_1_1_Grid <- data.frame(tol = NA, dif_otimo = NA, dif = NA, 
                              norma2_grad = NA, av_Q = NA, n_it = NA,
                              av_grad = NA, tempo = NA, dif_q_ot = NA, 
                              prod_grad = NA,
                              ls = "Grid")
beta_init <- c(500, 100)
gabarito <- c(2.3894212918E+02, 5.5015643181E-04)
av_Q <- 0
av_grad <- 0
n_it <- 0
tempo <- 0
prod_grad <- 0
for(tol in 2:-4){
  print(tol)
  tempo_i <- as.numeric(system.time({
    aux <- gradient_descent(y = Misra1a$y,
                            x =  Misra1a$x/1e6,
                            beta_init = beta_init,
                            f =  f1,
                            df_dbeta = df1_dbeta,
                            tol = 10^tol, tipo_line_search = "Grid",
                            Q = Q1_beta, exp_max = 5, exp_min = -1)
  })['elapsed'])
  beta_init <- aux$beta
  av_Q <- av_Q + aux$av_Q
  av_grad <- av_grad + aux$av_grad
  prod_grad <- prod_grad + aux$prod_grad
  n_it <- n_it + aux$n_it
  tempo <- tempo + tempo_i
  
  res_gd_1_1_Grid <- rbind(res_gd_1_1_Grid,
                           data.frame(tol = tol,
                                      dif_otimo = sqrt(sum((gabarito - beta_init*c(1, 1e-6))^2)), 
                                      dif = aux$dif, 
                                      norma2_grad = sqrt(sum(aux$grad^2)),
                                      av_Q = av_Q, 
                                      n_it = n_it,
                                      av_grad = av_grad,
                                      prod_grad = prod_grad,
                                      tempo = tempo,
                                      dif_q_ot = Q1_beta(beta_init) - Q1_beta(gabarito*c(1, 1e6)),
                                      ls = "Grid"))
}

res_gd_1_1_gd <- data.frame(tol = NA, dif_otimo = NA, dif = NA, 
                            norma2_grad = NA, av_Q = NA, n_it = NA,
                            av_grad = NA, tempo = NA, dif_q_ot = NA, 
                            prod_grad = NA,
                            ls = "Gradient Descent")
beta_init <- c(500, 100)
gabarito <- c(2.3894212918E+02, 5.5015643181E-04)
av_Q <- 0
av_grad <- 0
n_it <- 0
tempo <- 0
prod_grad <- 0
for(tol in 2:-4){
  print(tol)
  tempo_i <- as.numeric(system.time({
    aux <- gradient_descent(y = Misra1a$y,
                            x =  Misra1a$x/1e6,
                            beta_init = beta_init,
                            f =  f1,
                            df_dbeta = df1_dbeta,
                            tol = 10^tol, tipo_line_search = "gd",
                            gamma_init = 10^-3,
                            Q = Q1_beta)
  })['elapsed'])
  beta_init <- aux$beta
  av_Q <- av_Q + aux$av_Q
  av_grad <- av_grad + aux$av_grad
  prod_grad <- prod_grad + aux$prod_grad
  n_it <- n_it + aux$n_it
  tempo <- tempo + tempo_i
  
  res_gd_1_1_gd <- rbind(res_gd_1_1_gd,
                         data.frame(tol = tol,
                                    dif_otimo = sqrt(sum((gabarito - beta_init*c(1, 1e-6))^2)), 
                                    dif = aux$dif, 
                                    norma2_grad = sqrt(sum(aux$grad^2)),
                                    av_Q = av_Q, 
                                    n_it = n_it,
                                    av_grad = av_grad,
                                    prod_grad = prod_grad,
                                    tempo = tempo,
                                    dif_q_ot = Q1_beta(beta_init) - Q1_beta(gabarito*c(1, 1e6)),
                                    ls = "Gradient Descent"))
}

res_gd_1_1_sem <- data.frame(tol = NA, dif_otimo = NA, dif = NA, 
                             norma2_grad = NA, av_Q = NA, n_it = NA,
                             av_grad = NA, tempo = NA, dif_q_ot = NA,
                             prod_grad = NA,
                             ls = "Sem (gamma = 1)")
beta_init <- c(500, 100)
gabarito <- c(2.3894212918E+02, 5.5015643181E-04)
av_Q <- 0
av_grad <- 0
n_it <- 0
tempo <- 0
prod_grad <- 0

for(tol in 2:-4){
  print(tol)
  tempo_i <- as.numeric(system.time({
    aux <- gradient_descent(y = Misra1a$y,
                            x =  Misra1a$x/1e6,
                            beta_init = beta_init,
                            f =  f1,
                            df_dbeta = df1_dbeta,
                            tol = 10^tol, tipo_line_search = F,
                            gamma_init = 1,
                            Q = Q1_beta)
  })['elapsed'])
  beta_init <- aux$beta
  av_Q <- av_Q + aux$av_Q
  av_grad <- av_grad + aux$av_grad
  n_it <- n_it + aux$n_it
  tempo <- tempo + tempo_i
  prod_grad <- prod_grad + aux$prod_grad
  
  res_gd_1_1_sem <- rbind(res_gd_1_1_sem,
                          data.frame(tol = tol,
                                     dif_otimo = sqrt(sum((gabarito - beta_init*c(1, 1e-6))^2)), 
                                     dif = aux$dif, 
                                     norma2_grad = sqrt(sum(aux$grad^2)),
                                     av_Q = av_Q, 
                                     n_it = n_it,
                                     av_grad = av_grad,
                                     tempo = tempo,
                                     dif_q_ot = Q1_beta(beta_init) - Q1_beta(gabarito*c(1, 1e6)),
                                     prod_grad = prod_grad,
                                     ls = "Sem (gamma = 1)"))
}

gd11 <- rbind(res_gd_1_1_Grid,
              res_gd_1_1_gd,
              res_gd_1_1_sem) %>% 
  transmute("Tolerância" = tol,
            "Log10 do nº de avaliações do gradiente" = log(av_grad, 10),
            "Log10 do nº de iterações" = log(n_it, 10),
            "Log10 do nº de avaliações de Q1 + 1" =log(av_Q + 1, 10),
            
            "Log10 da distância ao mínimo global" = log(dif_otimo, 10),
            "Log10 da norma 2 do gradiente" = log(norma2_grad, 10),
            "Log10 da diferença para a perda ótima" = log(dif_q_ot, 10),
            "Tempo (s) até a convergência" = tempo,
            "Método de line-search" = ls,
  ) %>% 
  pivot_longer(cols = 2:8, names_to = "kpi", values_to = "valor") %>% 
  ggplot(aes(x = `Tolerância`, y = valor, color = `Método de line-search`,
             shape = `Método de line-search`))+
  geom_point()+
  geom_line()+
  facet_wrap(~kpi, scales = "free")+
  scale_x_reverse(breaks = -4:2, labels = paste0("10^", -4:2))+
  theme_bw()+
  ggtitle('Ponto inicial: (500, 0.0001)')+
  guides(color = guide_legend(position = "top", 
                              title = "Método de line-search",
                              title.position="top"))+
  theme(plot.title = element_text(size=18, hjust = 0.5),
        axis.title = element_text(size=14),
        legend.title = element_text(size=14))

# ponto inicial 2
res_gd_1_2_Grid <- data.frame(tol = NA, dif_otimo = NA, dif = NA, 
                              norma2_grad = NA, av_Q = NA, n_it = NA,
                              av_grad = NA, tempo = NA, dif_q_ot = NA, 
                              ls = "Grid")
beta_init <- c(250,500)
gabarito <- c(2.3894212918E+02, 5.5015643181E-04)
av_Q <- 0
av_grad <- 0
n_it <- 0
tempo <- 0
for(tol in 2:-4){
  print(tol)
  tempo_i <- as.numeric(system.time({
    aux <- gradient_descent(y = Misra1a$y,
                            x =  Misra1a$x/1e6,
                            beta_init = beta_init,
                            f =  f1,
                            df_dbeta = df1_dbeta,
                            tol = 10^tol, tipo_line_search = "Grid",
                            Q = Q1_beta, exp_max = 5, exp_min = -1)
  })['elapsed'])
  beta_init <- aux$beta
  av_Q <- av_Q + aux$av_Q
  av_grad <- av_grad + aux$av_grad
  n_it <- n_it + aux$n_it
  tempo <- tempo + tempo_i
  
  res_gd_1_2_Grid <- rbind(res_gd_1_2_Grid,
                           data.frame(tol = tol,
                                      dif_otimo = sqrt(sum(gabarito - beta_init*c(1, 1e-6))^2), 
                                      dif = aux$dif, 
                                      norma2_grad = sqrt(sum(aux$grad)^2),
                                      av_Q = av_Q, 
                                      n_it = n_it,
                                      av_grad = av_grad,
                                      tempo = tempo,
                                      dif_q_ot = Q1_beta(beta_init) - Q1_beta(gabarito*c(1, 1e6)),
                                      ls = "Grid"))
}

res_gd_1_2_gd <- data.frame(tol = NA, dif_otimo = NA, dif = NA, 
                            norma2_grad = NA, av_Q = NA, n_it = NA,
                            av_grad = NA, tempo = NA, dif_q_ot = NA, 
                            ls = "Gradient Descent")
beta_init <- c(250,500)
gabarito <- c(2.3894212918E+02, 5.5015643181E-04)
av_Q <- 0
av_grad <- 0
n_it <- 0
tempo <- 0
for(tol in 2:-4){
  print(tol)
  tempo_i <- as.numeric(system.time({
    aux <- gradient_descent(y = Misra1a$y,
                            x =  Misra1a$x/1e6,
                            beta_init = beta_init,
                            f =  f1,
                            df_dbeta = df1_dbeta,
                            tol = 10^tol, tipo_line_search = "gd",
                            gamma_init = 10^-3,
                            Q = Q1_beta)
  })['elapsed'])
  beta_init <- aux$beta
  av_Q <- av_Q + aux$av_Q
  av_grad <- av_grad + aux$av_grad
  n_it <- n_it + aux$n_it
  tempo <- tempo + tempo_i
  
  res_gd_1_2_gd <- rbind(res_gd_1_2_gd,
                         data.frame(tol = tol,
                                    dif_otimo = sqrt(sum(gabarito - beta_init*c(1, 1e-6))^2), 
                                    dif = aux$dif, 
                                    norma2_grad = sqrt(sum(aux$grad)^2),
                                    av_Q = av_Q, 
                                    n_it = n_it,
                                    av_grad = av_grad,
                                    tempo = tempo,
                                    dif_q_ot = Q1_beta(beta_init) - Q1_beta(gabarito*c(1, 1e6)),
                                    ls = "Gradient Descent"))
}

res_gd_1_2_sem <- data.frame(tol = NA, dif_otimo = NA, dif = NA, 
                             norma2_grad = NA, av_Q = NA, n_it = NA,
                             av_grad = NA, tempo = NA, dif_q_ot = NA,
                             ls = NA)
beta_init <- c(250,500)
gabarito <- c(2.3894212918E+02, 5.5015643181E-04)
av_Q <- 0
av_grad <- 0
n_it <- 0
tempo <- 0
for(tol in 2:-4){
  print(tol)
  tempo_i <- as.numeric(system.time({
    aux <- gradient_descent(y = Misra1a$y,
                            x =  Misra1a$x/1e6,
                            beta_init = beta_init,
                            f =  f1,
                            df_dbeta = df1_dbeta,
                            tol = 10^tol, tipo_line_search = F,
                            gamma_init = 1,
                            Q = Q1_beta)
  })['elapsed'])
  beta_init <- aux$beta
  av_Q <- av_Q + aux$av_Q
  av_grad <- av_grad + aux$av_grad
  n_it <- n_it + aux$n_it
  tempo <- tempo + tempo_i
  
  res_gd_1_2_sem <- rbind(res_gd_1_2_sem,
                          data.frame(tol = tol,
                                     dif_otimo = sqrt(sum(gabarito - beta_init*c(1, 1e-6))^2), 
                                     dif = aux$dif, 
                                     norma2_grad = sqrt(sum(aux$grad)^2),
                                     av_Q = av_Q, 
                                     n_it = n_it,
                                     av_grad = av_grad,
                                     tempo = tempo,
                                     dif_q_ot = Q1_beta(beta_init) - Q1_beta(gabarito*c(1, 1e6)),
                                     ls = "Sem (gamma = 1)"))
}

gd12 <- rbind(res_gd_1_2_Grid,
              res_gd_1_2_gd,
              res_gd_1_2_sem) %>% 
  drop_na() %>% 
  transmute("Tolerância" = tol,
            "Log10 do nº de avaliações do gradiente" = log(av_grad, 10),
            "Log10 do nº de iterações" = log(n_it, 10),
            "Log10 do nº de avaliações de Q1 + 1" =log(av_Q + 1, 10),
            
            "Log10 da distância ao mínimo global" = log(dif_otimo, 10),
            "Log10 da norma 2 do gradiente" = log(norma2_grad, 10),
            "Log10 da diferença para a perda ótima" = log(dif_q_ot, 10),
            "Tempo (s) até a convergência" = tempo,
            "Método de line-search" = ls,
  ) %>% 
  pivot_longer(cols = 2:8, names_to = "kpi", values_to = "valor") %>% 
  ggplot(aes(x = `Tolerância`, y = valor, color = `Método de line-search`,
             shape = `Método de line-search`))+
  geom_point()+
  geom_line()+
  facet_wrap(~kpi, scales = "free")+
  scale_x_reverse(breaks = -4:2, labels = paste0("10^", -4:2))+
  theme_bw()+
  ggtitle('Ponto inicial: (500, 0.0005)')+
  guides(color = guide_legend(position = "top", 
                              title = "Método de line-search",
                              title.position="top"))+
  theme(plot.title = element_text(size=18, hjust = 0.5),
        axis.title = element_text(size=14),
        legend.title = element_text(size=14))

# fim do problema 1

# começo do problema 2
Thurber <- data.frame(y = c(80.574, 84.248, 87.264, 87.195, 89.076, 89.608,
                            89.868, 90.101, 92.405, 95.854, 100.696, 101.06,
                            401.672, 390.724, 567.534, 635.316, 733.054,
                            759.087, 894.206, 990.785, 1090.109, 1080.914,
                            1122.643, 1178.351, 1260.531, 1273.514, 1288.339,
                            1327.543, 1353.863, 1414.509, 1425.208, 1421.384,
                            1442.962, 1464.35, 1468.705, 1447.894, 1457.628),
                      x = c(-3.067,-2.981,-2.921,-2.912,-2.84,-2.797,-2.702,-2.699,-2.633,-2.481,-2.363,-2.322,-1.501,-1.46,-1.274,-1.212,-1.1,-1.046,-0.915,-0.714,-0.566,-0.545,-0.4,-0.309,-0.109,-0.103, 0.01, 0.119, 0.377, 0.79, 0.963,
                            1.006, 1.115, 1.572, 1.841, 2.047, 2.2))

Q2_beta <- function(beta){
  Q <- 0
  for(i in 1:nrow(Thurber)){
    Q <- Q + (Thurber$y[i] - f2(Thurber$x[i], beta))^2
  }
  Q
}

f2 <- function(x, beta){
  betas_1 <- matrix(beta[1:4], nrow = 4)
  X1 <- matrix(c(1*1000, x*1000, x^2*400, x^3*40), ncol = 4)/100
  
  betas_2 <- matrix(beta[5:7], nrow = 3)
  X2 <-  matrix(c(x*0.7, x^2*0.3, x^3*0.03), ncol = 3)/100
  
  f2 <- as.numeric((1 + X2 %*% betas_2)^(-1)*X1 %*% betas_1)
  if(is.finite(f2)){
    return(f2)
  } else{
    return(as.numeric((0.0001)^(-1)*X1 %*% betas_1))
  }
}

df2_dbeta <- function(x, beta){
  betas_1 <- matrix(beta[1:4], nrow = 4)
  X1 <- matrix(c(1*1000, x*1000, x^2*400, x^3*40), ncol = 4)/100
  
  betas_2 <- matrix(beta[5:7], nrow = 3)
  X2 <-  matrix(c(x*0.7, x^2*0.3, x^3*0.03), ncol = 3)/100
  
  der <- matrix(c(as.numeric((1 + X2 %*% betas_2)^(-1))*X1, 
                  -as.numeric((1 + X2 %*% betas_2))^(-2)*X1 %*% betas_1 %*% X2), ncol = 7)
  if(sum(!is.finite(der))){
    return(matrix(c(as.numeric((0.0001)^(-1))*X1, 
                    -as.numeric((0.0001))^(-2)*X1 %*% betas_1 %*% X2), ncol = 7))
  } else{
    return(der)
  }
  
}
Q2_beta <- function(beta){
  Q <- 0
  for(i in 1:nrow(Thurber)){
    Q <- Q + (Thurber$y[i] - f2(Thurber$x[i], beta))^2
  }
  Q
}

res_gd_2_1_Grid <- data.frame(tol = NA, dif_otimo = NA, dif = NA, 
                              norma2_grad = NA, av_Q = NA, n_it = NA,
                              av_grad = NA, tempo = NA, dif_q_ot = NA, 
                              ls = "Grid")
beta_init <- rep(100, 7)
gabarito <- c(1.2881396800E+03,
              1.4910792535E+03,
              5.8323836877E+02,
              7.5416644291E+01,
              9.6629502864E-01,
              3.9797285797E-01,
              4.9727297349E-02)
av_Q <- 0
av_grad <- 0
n_it <- 0
tempo <- 0
for(tol in 2:-3){
  print("Tolerância")
  print(tol)
  tempo_i <- as.numeric(system.time({
    aux <- gradient_descent(y = Thurber$y,
                            x = Thurber$x,
                            beta_init = beta_init,
                            f =  f2,
                            df_dbeta = df2_dbeta,
                            tol = 10^tol, tipo_line_search = "Grid",
                            Q = Q2_beta, exp_max = -10, exp_min = -20)
  })['elapsed'])
  beta_init <- aux$beta
  av_Q <- av_Q + aux$av_Q
  av_grad <- av_grad + aux$av_grad
  n_it <- n_it + aux$n_it
  tempo <- tempo + tempo_i
  
  res_gd_2_1_Grid <- rbind(res_gd_2_1_Grid,
                           data.frame(tol = tol,
                                      dif_otimo = sqrt(sum((gabarito - beta_init*c(1000,1000,400,40,0.7,0.3,0.03))^2)), 
                                      dif = aux$dif, 
                                      norma2_grad = sqrt(sum(aux$grad^2)),
                                      av_Q = av_Q, 
                                      n_it = n_it,
                                      av_grad = av_grad,
                                      tempo = tempo,
                                      dif_q_ot = Q2_beta(beta_init) -  Q2_beta(gabarito/(c(1000,1000,400,40,0.7,0.3,0.03))/0.01),
                                      ls = "Grid"))
}

res_gd_2_1_gd <- data.frame(tol = NA, dif_otimo = NA, dif = NA, 
                            norma2_grad = NA, av_Q = NA, n_it = NA,
                            av_grad = NA, tempo = NA, dif_q_ot = NA, 
                            ls = "Gradient Descent")
beta_init <- rep(100, 7)
gabarito <- c(1.2881396800E+03,
              1.4910792535E+03,
              5.8323836877E+02,
              7.5416644291E+01,
              9.6629502864E-01,
              3.9797285797E-01,
              4.9727297349E-02)
av_Q <- 0
av_grad <- 0
n_it <- 0
tempo <- 0
for(tol in 2:-3){
  print("Tolerância")
  print(tol)
  tempo_i <- as.numeric(system.time({
    aux <- gradient_descent(y = Thurber$y,
                            x = Thurber$x,
                            beta_init = beta_init,
                            f =  f2,
                            df_dbeta = df2_dbeta,
                            tol = 10^tol, tipo_line_search = 'gd',
                            gamma_init = 10^-17,
                            Q = Q2_beta)
  })['elapsed'])
  beta_init <- aux$beta
  av_Q <- av_Q + aux$av_Q
  av_grad <- av_grad + aux$av_grad
  n_it <- n_it + aux$n_it
  tempo <- tempo + tempo_i
  
  res_gd_2_1_gd <- rbind(res_gd_2_1_gd,
                         data.frame(tol = tol,
                                    dif_otimo = sqrt(sum((gabarito - beta_init*c(1000,1000,400,40,0.7,0.3,0.03))^2)), 
                                    dif = aux$dif, 
                                    norma2_grad = sqrt(sum(aux$grad^2)),
                                    av_Q = av_Q, 
                                    n_it = n_it,
                                    av_grad = av_grad,
                                    tempo = tempo,
                                    dif_q_ot = Q2_beta(beta_init) -  Q2_beta(gabarito/(c(1000,1000,400,40,0.7,0.3,0.03))/0.01),
                                    ls = "Gradient Descent"))
}

res_gd_2_1_sem <- data.frame(tol = NA, dif_otimo = NA, dif = NA, 
                             norma2_grad = NA, av_Q = NA, n_it = NA,
                             av_grad = NA, tempo = NA, dif_q_ot = NA, 
                             ls = "Sem (Gamma = 10^-6)")
beta_init <- rep(100, 7)
gabarito <- c(1.2881396800E+03,
              1.4910792535E+03,
              5.8323836877E+02,
              7.5416644291E+01,
              9.6629502864E-01,
              3.9797285797E-01,
              4.9727297349E-02)
av_Q <- 0
av_grad <- 0
n_it <- 0
tempo <- 0
for(tol in 2:-3){
  print("Tolerância")
  print(tol)
  tempo_i <- as.numeric(system.time({
    aux <- gradient_descent(y = Thurber$y,
                            x = Thurber$x,
                            beta_init = beta_init,
                            f =  f2,
                            df_dbeta = df2_dbeta,
                            tol = 10^tol, tipo_line_search = F,
                            gamma_init = 1e-6,
                            Q = Q2_beta)
  })['elapsed'])
  beta_init <- aux$beta
  av_Q <- av_Q + aux$av_Q
  av_grad <- av_grad + aux$av_grad
  n_it <- n_it + aux$n_it
  tempo <- tempo + tempo_i
  
  res_gd_2_1_sem <- rbind(res_gd_2_1_sem,
                          data.frame(tol = tol,
                                     dif_otimo = sqrt(sum((gabarito - beta_init*c(1000,1000,400,40,0.7,0.3,0.03))^2)), 
                                     dif = aux$dif, 
                                     norma2_grad = sqrt(sum(aux$grad^2)),
                                     av_Q = av_Q, 
                                     n_it = n_it,
                                     av_grad = av_grad,
                                     tempo = tempo,
                                     dif_q_ot = Q2_beta(beta_init) -  Q2_beta(gabarito/(c(1000,1000,400,40,0.7,0.3,0.03))/0.01),
                                     ls = "Sem (Gamma = 10^-6)"))
}

gd21 <- rbind(res_gd_2_1_Grid,
              res_gd_2_1_gd,
              res_gd_2_1_sem)%>% 
  drop_na() %>% 
  transmute("Tolerância" = tol,
            "Log10 do nº de avaliações do gradiente" = log(av_grad, 10),
            "Log10 do nº de iterações" = log(n_it, 10),
            "Log10 do nº de avaliações de Q2 + 1" =log(av_Q + 1, 10),
            
            "Log10 da distância ao mínimo global" = log(dif_otimo, 10),
            "Log10 da norma 2 do gradiente" = log(norma2_grad, 10),
            "Log10 da diferença para a perda ótima" = log(dif_q_ot, 10),
            "Tempo (s) até a convergência" = tempo,
            "Método de line-search" = ls,
  ) %>% 
  pivot_longer(cols = 2:8, names_to = "kpi", values_to = "valor") %>% 
  ggplot(aes(x = `Tolerância`, y = valor, color = `Método de line-search`,
             shape = `Método de line-search`))+
  geom_point()+
  geom_line()+
  facet_wrap(~kpi, scales = "free")+
  scale_x_reverse(breaks = -4:2, labels = paste0("10^", -4:2))+
  theme_bw()+
  ggtitle('Ponto inicial: (1300,1500,500,75,1,0.4,0.05)')+
  guides(color = guide_legend(position = "top", 
                              title = "Método de line-search",
                              title.position="top"))+
  theme(plot.title = element_text(size=18, hjust = 0.5),
        axis.title = element_text(size=14),
        legend.title = element_text(size=14))

# Ponto inicial 2
Thurber <- data.frame(y = c(80.574, 84.248, 87.264, 87.195, 89.076, 89.608,
                            89.868, 90.101, 92.405, 95.854, 100.696, 101.06,
                            401.672, 390.724, 567.534, 635.316, 733.054,
                            759.087, 894.206, 990.785, 1090.109, 1080.914,
                            1122.643, 1178.351, 1260.531, 1273.514, 1288.339,
                            1327.543, 1353.863, 1414.509, 1425.208, 1421.384,
                            1442.962, 1464.35, 1468.705, 1447.894, 1457.628),
                      x = c(-3.067,-2.981,-2.921,-2.912,-2.84,-2.797,-2.702,-2.699,-2.633,-2.481,-2.363,-2.322,-1.501,-1.46,-1.274,-1.212,-1.1,-1.046,-0.915,-0.714,-0.566,-0.545,-0.4,-0.309,-0.109,-0.103, 0.01, 0.119, 0.377, 0.79, 0.963,
                            1.006, 1.115, 1.572, 1.841, 2.047, 2.2))

Q2_beta <- function(beta){
  Q <- 0
  for(i in 1:nrow(Thurber)){
    Q <- Q + (Thurber$y[i] - f2(Thurber$x[i], beta))^2
  }
  Q
}

f2 <- function(x, beta){
  betas_1 <- matrix(beta[1:4], nrow = 4)
  X1 <- matrix(c(1*1300, x*1500, x^2*500, x^3*75), ncol = 4)/100
  
  betas_2 <- matrix(beta[5:7], nrow = 3)
  X2 <-  matrix(c(x*1, x^2*0.4, x^3*0.05), ncol = 3)/100
  
  f2 <- as.numeric((1 + X2 %*% betas_2)^(-1)*X1 %*% betas_1)
  if(is.finite(f2)){
    return(f2)
  } else{
    return(as.numeric((0.0001)^(-1)*X1 %*% betas_1))
  }
}

df2_dbeta <- function(x, beta){
  betas_1 <- matrix(beta[1:4], nrow = 4)
  X1 <- matrix(c(1*1300, x*1500, x^2*500, x^3*75), ncol = 4)/100
  
  betas_2 <- matrix(beta[5:7], nrow = 3)
  X2 <-  matrix(c(x*1, x^2*0.4, x^3*0.05), ncol = 3)/100
  
  der <- matrix(c(as.numeric((1 + X2 %*% betas_2)^(-1))*X1, 
                  -as.numeric((1 + X2 %*% betas_2))^(-2)*X1 %*% betas_1 %*% X2), ncol = 7)
  if(sum(!is.finite(der))){
    return(matrix(c(as.numeric((0.0001)^(-1))*X1, 
                    -as.numeric((0.0001))^(-2)*X1 %*% betas_1 %*% X2), ncol = 7))
  } else{
    return(der)
  }
  
}
Q2_beta <- function(beta){
  Q <- 0
  for(i in 1:nrow(Thurber)){
    Q <- Q + (Thurber$y[i] - f2(Thurber$x[i], beta))^2
  }
  Q
}

res_gd_2_2_Grid <- data.frame(tol = NA, dif_otimo = NA, dif = NA, 
                              norma2_grad = NA, av_Q = NA, n_it = NA,
                              av_grad = NA, tempo = NA, dif_q_ot = NA, 
                              ls = "Grid")
beta_init <- rep(100, 7)
gabarito <- c(1.2881396800E+03,
              1.4910792535E+03,
              5.8323836877E+02,
              7.5416644291E+01,
              9.6629502864E-01,
              3.9797285797E-01,
              4.9727297349E-02)
av_Q <- 0
av_grad <- 0
n_it <- 0
tempo <- 0
for(tol in 2:-3){
  print("Tolerância")
  print(tol)
  tempo_i <- as.numeric(system.time({
    aux <- gradient_descent(y = Thurber$y,
                            x = Thurber$x,
                            beta_init = beta_init,
                            f =  f2,
                            df_dbeta = df2_dbeta,
                            tol = 10^tol, tipo_line_search = "Grid",
                            Q = Q2_beta, exp_max = -10, exp_min = -20)
  })['elapsed'])
  beta_init <- aux$beta
  av_Q <- av_Q + aux$av_Q
  av_grad <- av_grad + aux$av_grad
  n_it <- n_it + aux$n_it
  tempo <- tempo + tempo_i
  
  res_gd_2_2_Grid <- rbind(res_gd_2_2_Grid,
                           data.frame(tol = tol,
                                      dif_otimo = sqrt(sum((gabarito - beta_init*c(1300,1500,500,75,1,0.4,0.05))^2)), 
                                      dif = aux$dif, 
                                      norma2_grad = sqrt(sum(aux$grad^2)),
                                      av_Q = av_Q, 
                                      n_it = n_it,
                                      av_grad = av_grad,
                                      tempo = tempo,
                                      dif_q_ot = Q2_beta(beta_init) -  Q2_beta(gabarito/(c(1300,1500,500,75,1,0.4,0.05))/0.01),
                                      ls = "Grid"))
}

res_gd_2_2_gd <- data.frame(tol = NA, dif_otimo = NA, dif = NA, 
                            norma2_grad = NA, av_Q = NA, n_it = NA,
                            av_grad = NA, tempo = NA, dif_q_ot = NA, 
                            ls = "Gradient Descent")
beta_init <- rep(100, 7)
gabarito <- c(1.2881396800E+03,
              1.4910792535E+03,
              5.8323836877E+02,
              7.5416644291E+01,
              9.6629502864E-01,
              3.9797285797E-01,
              4.9727297349E-02)
av_Q <- 0
av_grad <- 0
n_it <- 0
tempo <- 0
for(tol in 2:-3){
  print("Tolerância")
  print(tol)
  tempo_i <- as.numeric(system.time({
    aux <- gradient_descent(y = Thurber$y,
                            x = Thurber$x,
                            beta_init = beta_init,
                            f =  f2,
                            df_dbeta = df2_dbeta,
                            tol = 10^tol, tipo_line_search = 'gd',
                            gamma_init = 10^-17,
                            Q = Q2_beta)
  })['elapsed'])
  beta_init <- aux$beta
  av_Q <- av_Q + aux$av_Q
  av_grad <- av_grad + aux$av_grad
  n_it <- n_it + aux$n_it
  tempo <- tempo + tempo_i
  
  res_gd_2_2_gd <- rbind(res_gd_2_2_gd,
                         data.frame(tol = tol,
                                    dif_otimo = sqrt(sum((gabarito - beta_init*c(1300,1500,500,75,1,0.4,0.05))^2)), 
                                    dif = aux$dif, 
                                    norma2_grad = sqrt(sum(aux$grad^2)),
                                    av_Q = av_Q, 
                                    n_it = n_it,
                                    av_grad = av_grad,
                                    tempo = tempo,
                                    dif_q_ot = Q2_beta(beta_init) -  Q2_beta(gabarito/(c(1300,1500,500,75,1,0.4,0.05))/0.01),
                                    ls = "Gradient Descent"))
}

res_gd_2_2_sem <- data.frame(tol = NA, dif_otimo = NA, dif = NA, 
                             norma2_grad = NA, av_Q = NA, n_it = NA,
                             av_grad = NA, tempo = NA, dif_q_ot = NA, 
                             ls = "Sem (gamma = 10^-6)")
beta_init <- rep(100, 7)
gabarito <- c(1.2881396800E+03,
              1.4910792535E+03,
              5.8323836877E+02,
              7.5416644291E+01,
              9.6629502864E-01,
              3.9797285797E-01,
              4.9727297349E-02)
av_Q <- 0
av_grad <- 0
n_it <- 0
tempo <- 0
for(tol in 2:-3){
  print("Tolerância")
  print(tol)
  tempo_i <- as.numeric(system.time({
    aux <- gradient_descent(y = Thurber$y,
                            x = Thurber$x,
                            beta_init = beta_init,
                            f =  f2,
                            df_dbeta = df2_dbeta,
                            tol = 10^tol, tipo_line_search = F,
                            gamma_init = 1e-6,
                            Q = Q2_beta)
  })['elapsed'])
  beta_init <- aux$beta
  av_Q <- av_Q + aux$av_Q
  av_grad <- av_grad + aux$av_grad
  n_it <- n_it + aux$n_it
  tempo <- tempo + tempo_i
  
  res_gd_2_2_sem <- rbind(res_gd_2_2_sem,
                          data.frame(tol = tol,
                                     dif_otimo = sqrt(sum((gabarito - beta_init*c(1300,1500,500,75,1,0.4,0.05))^2)), 
                                     dif = aux$dif, 
                                     norma2_grad = sqrt(sum(aux$grad^2)),
                                     av_Q = av_Q, 
                                     n_it = n_it,
                                     av_grad = av_grad,
                                     tempo = tempo,
                                     dif_q_ot = Q2_beta(beta_init) -  Q2_beta(gabarito/(c(1300,1500,500,75,1,0.4,0.05))/0.01),
                                     ls = "Sem (Gamma = 10^-6)"))
}

gd22 <- rbind(res_gd_2_2_Grid,
              res_gd_2_2_gd,
              res_gd_2_2_sem)%>% 
  drop_na() %>% 
  transmute("Tolerância" = tol,
            "Log10 do nº de avaliações do gradiente" = log(av_grad, 10),
            "Log10 do nº de iterações" = log(n_it, 10),
            "Log10 do nº de avaliações de Q2 + 1" =log(av_Q + 1, 10),
            
            "Log10 da distância ao mínimo global" = log(dif_otimo, 10),
            "Log10 da norma 2 do gradiente" = log(norma2_grad, 10),
            "Log10 da diferença para a perda ótima" = log(dif_q_ot, 10),
            "Tempo (s) até a convergência" = tempo,
            "Método de line-search" = ls,
  ) %>% 
  pivot_longer(cols = 2:8, names_to = "kpi", values_to = "valor") %>% 
  ggplot(aes(x = `Tolerância`, y = valor, color = `Método de line-search`,
             shape = `Método de line-search`))+
  geom_point()+
  geom_line()+
  facet_wrap(~kpi, scales = "free")+
  scale_x_reverse(breaks = -4:2, labels = paste0("10^", -4:2))+
  theme_bw()+
  ggtitle('Ponto inicial: (1000, 1000, 400, 40, 0.7, 0.3, 0.03)')+
  guides(color = guide_legend(position = "top", 
                              title = "Método de line-search",
                              title.position="top"))+
  theme(plot.title = element_text(size=18, hjust = 0.5),
        axis.title = element_text(size=14),
        legend.title = element_text(size=14))


# gd11+gd12
# gd21+gd22

