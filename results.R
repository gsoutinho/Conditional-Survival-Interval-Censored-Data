load('breastCancer.RData')

breastCancer[, c("rectime", "censrec", "survtime", "censdead", "left_rec", 
                 "right_rec", "left", "right", "age", "menopause", "hormone", "size", "grade", "nodes", 
                 "prog_recp", "estrg_recp")]

dim(breastCancer)

bandwidth_silverman <- function(U){
  n <- length(U)
  s <- min(sd(U), IQR(U)/1.349)
  h <- 1.06 * s * n^(-1/5)
  h
}

h <- bandwidth_silverman(breastCancer$nodes)  

times <- seq(0, 3650, by = 30)

fit_nodes_2 <- cond_surv_ic(
  L = breastCancer$left,
  R = breastCancer$right,
  Z = breastCancer$nodes,
  z0 = 2,
  times = times, 
  method='Turnbull'
)


fit_nodes_5 <- cond_surv_ic(
  L = breastCancer$left,
  R = breastCancer$right,
  Z = breastCancer$nodes,
  z0 = 5,
  times = times, 
  method='Turnbull'
)

fit_nodes_8 <- cond_surv_ic(
  L = breastCancer$left,
  R = breastCancer$right,
  Z = breastCancer$nodes,
  z0 = 8,
  times = times, 
  method='Turnbull'
)

plot(
  fit_nodes_2$survival$time / 365,
  fit_nodes_2$survival$survival,
  type = "s",
  lwd = 2,
  ylim = c(0, 1),
  xlab = "Time since surgery (years)",
  ylab = "Recurrence-free survival probability",
  main = "", 
  col = 1
)


lines(
  fit_nodes_5$survival$time / 365,
  fit_nodes_5$survival$survival,
  type = "s",
  lwd = 2,
  lty = 2, 
  col=2
)

lines(
  fit_nodes_8$survival$time / 365,
  fit_nodes_8$survival$survival,
  type = "s",
  lwd = 2,
  lty = 3,
  col = 4
)

legend(
  "topright",
  legend = c(
    "nodes = 2",
    "nodes = 5",
    "nodes = 8"
  ),
  lty = c(1, 2, 3),
  lwd = 2,
  bty = "n",
  col = c(1,2,4),
  cex=.7
)


midpoint_beran_nodes_2 <- cond_surv_ic(
  L = breastCancer$left,
  R = breastCancer$right,
  Z = breastCancer$nodes,
  z0 = 2,
  times = times, 
  method='Beran',
  h=h
)


midpoint_beran_nodes_5 <- cond_surv_ic(
  L = breastCancer$left,
  R = breastCancer$right,
  Z = breastCancer$nodes,
  z0 = 5,
  times = times, 
  method='Beran',
  h=h
)

midpoint_beran_nodes_8 <- cond_surv_ic(
  L = breastCancer$left,
  R = breastCancer$right,
  Z = breastCancer$nodes,
  z0 = 8,
  times = times, 
  method='Beran',
  h=h
)

plot(
  times / 365,
  midpoint_beran_nodes_2,
  type = "s",
  lwd = 2,
  ylim = c(0, 1),
  xlab = "Time since surgery (years)",
  ylab = "Recurrence-free survival probability",
  main = "", 
  col = 1
)

lines(
  times / 365,
  midpoint_beran_nodes_5,
  type = "s",
  lwd = 2,
  lty = 2, 
  col=2
)

lines(
  times / 365,
  midpoint_beran_nodes_8,
  type = "s",
  lwd = 2,
  lty = 3,
  col = 4
)

legend(
  "topright",
  legend = c(
    "nodes = 2",
    "nodes = 5",
    "nodes = 8"
  ),
  lty = c(1, 2, 3),
  lwd = 2,
  bty = "n",
  col = c(1,2,4),
  cex=.7
)


#-------------------------------------------------------------------------------

times <- seq(0, 3650, by = 30)

fit_1_6 <-cond_aZb_surv_ic(
  L = breastCancer$left,
  R = breastCancer$right,
  T = breastCancer$survtime,
  Delta = breastCancer$censdead,
  a = 30,
  b = 180,
  times = times,
  method='Turnbull'
)


fit_12_24 <-cond_aZb_surv_ic(
  L = breastCancer$left,
  R = breastCancer$right,
  T = breastCancer$survtime,
  Delta = breastCancer$censdead,
  a = 360,
  b = 720,
  times = times,
  method='Turnbull'
)

fit_24_36 <-cond_aZb_surv_ic(
  L = breastCancer$left,
  R = breastCancer$right,
  T = breastCancer$survtime,
  Delta = breastCancer$censdead,
  a = 720,
  b = 1080,
  times = times,
  method='Turnbull'
)


plot(
  fit_1_6$survival$time / 365,
  fit_1_6$survival$survival,
  type = "s",
  lwd = 2,
  ylim = c(0, 1),
  xlab = "Time since surgery (years)",
  ylab = "Conditional survival probability"#,
  #main = "Conditional survival according to recurrence time"
)


lines(
  fit_12_24$survival$time / 365,
  fit_12_24$survival$survival,
  type = "s",
  lwd = 2,
  col=2,
  lty = 2
)

lines(
  fit_24_36$survival$time / 365,
  fit_24_36$survival$survival,
  type = "s",
  lwd = 2,
  col=4,
  lty = 3
)


legend(
  "topright",
  legend = c(
    "1 < Z <= 6 months",
    "12 < Z <= 24 months",
    "24 < Z <= 36 months"
  ),
  lty = c(1, 2, 3), col=c(1, 2, 4),
  lwd = 2,
  bty = "n",
  cex=0.8
)


Beran_1_6<-cond_aZb_surv_ic(
  L = breastCancer$left,
  R = breastCancer$right,
  T = breastCancer$survtime,
  Delta = breastCancer$censdead,
  a = 30,
  b = 180,
  times = times,
  method='Beran'
)

Beran_12_24<-cond_aZb_surv_ic(
  L = breastCancer$left,
  R = breastCancer$right,
  T = breastCancer$survtime,
  Delta = breastCancer$censdead,
  a = 360,
  b = 720,
  times = times,
  method='Beran'
)

Beran_24_36<-cond_aZb_surv_ic(
  L = breastCancer$left,
  R = breastCancer$right,
  T = breastCancer$survtime,
  Delta = breastCancer$censdead,
  a = 720,
  b = 1080,
  times = times,
  method='Beran'
)


plot(
  fit_1_6$survival$time / 365,
  Beran_1_6,
  type = "s",
  lwd = 2,
  ylim = c(0, 1),
  xlab = "Time since surgery (years)",
  ylab = "Conditional survival probability"#,
  #main = "Conditional survival according to recurrence time"
)


lines(
  fit_12_24$survival$time / 365,
  Beran_12_24,
  type = "s",
  lwd = 2,
  col=2,
  lty = 2
)

lines(
  fit_24_36$survival$time / 365,
  Beran_24_36,
  type = "s",
  lwd = 2,
  col=4,
  lty = 3
)


legend(
  "topright",
  legend = c(
    "1 < Z <= 6 months",
    "12 < Z <= 24 months",
    "24 < Z <= 36 months"
  ),
  lty = c(1, 2, 3), col=c(1, 2, 4),
  lwd = 2,
  bty = "n",
  cex=0.8
)
