library(tidyverse)
library(ggplot2)
library(gridExtra)

'here it is assumed your in your working directory'

#load data
res_CG <- readRDS("result_CG.RData")
res_CHH <- readRDS("result_CHH.RData")
res_CHG <- readRDS("result_CHG.RData")

#vizualise bootstrap distribution
cg_boot_um <- as.data.table(res_CG$boot_est_u_m)
chh_boot_um <- as.data.table(res_CHH$boot_est_u_m)
chg_boot_um <- as.data.table(res_CHG$boot_est_u_m)

cg_boot_mu <- as.data.table(res_CG$boot_est_m_u)
chh_boot_mu <- as.data.table(res_CHH$boot_est_m_u)
chg_boot_mu <- as.data.table(res_CHG$boot_est_m_u)

#set theme
t <- theme_bw()

plot_cg_01 <- ggplot(data = cg_boot_um, aes(V1)) + geom_histogram() + labs(title= "CG: u->m", y = "", x= "") + t
plot_chh_01 <- ggplot(data = chh_boot_um, aes(V1)) + geom_histogram()+ labs(title= "CHH: u->m", y = "", x= "") + t 
plot_chg_01 <- ggplot(data = chg_boot_um, aes(V1)) + geom_histogram() + labs(title= "CHG: u->m", y = "", x= "") + t

plot_cg_10 <- ggplot(data = cg_boot_mu, aes(V1)) + geom_histogram() + labs(title= "CG: m->u", y = "", x= "") + t
plot_chh_10 <- ggplot(data = chh_boot_mu, aes(V1)) + geom_histogram() + labs(title= "CHH: m->u", y = "", x= "") + t
plot_chg_10 <- ggplot(data = chg_boot_mu, aes(V1)) + geom_histogram() + labs(title= "CHG: m->u", y = "", x= "") + t


boot_dist <- grid.arrange(plot_cg_01,plot_chh_01,plot_chg_01,plot_cg_10,plot_chh_10,plot_chg_10 , 
             ncol=3, nrow=2, top = "Bootstrap estimates for the different contexts")

#boxplot for transition probability from u to m:

d_u_m <- data.frame(group = rep(c("CG", "CHH", "CHG"), each = 1000), mean = rep(c(res_CG$mean_u_m, 
                                                                                  res_CHH$mean_u_m, res_CHG$mean_u_m), each = 1000))
d_u_m$value[1:1000] <- res_CG$boot_est_u_m
d_u_m$value[1001:2000] <- res_CHH$boot_est_u_m
d_u_m$value[2001:3000] <- res_CHG$boot_est_u_m

g0_u_m <- ggplot(d_u_m, aes(x = group, y = value))
g_box_u_m <- g0_u_m + geom_boxplot(fill = "grey", colour = "black") + theme_bw() +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", colour = "red") + 
  stat_summary(fun.y = mean, geom = "point", colour = "red") +ggtitle("Transition from u to m") 

#boxplot for transition probability from m to u:
d_m_u <- data.frame(group = rep(c("CG", "CHH", "CHG"), each = 1000), mean = rep(c(res_CG$mean_m_u, 
                                                                                  res_CHH$mean_m_u, res_CHG$mean_u_m), each = 1000))
d_m_u$value[1:1000] <- res_CG$boot_est_m_u
d_m_u$value[1001:2000] <- res_CHH$boot_est_m_u
d_m_u$value[2001:3000] <- res_CHG$boot_est_m_u

g0_m_u <- ggplot(d_m_u, aes(x = group, y = value))
g_box_m_u <- g0_m_u + geom_boxplot(fill = "grey", colour = "black") + theme_bw() +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", colour = "red") + 
  stat_summary(fun.y = mean, geom = "point", colour = "red") +ggtitle("Transition from m to u")

#figure with both of the boxplots:
grid.arrange(g_box_u_m , g_box_m_u, ncol=2, nrow=1)
