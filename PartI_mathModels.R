################################################################################
### Time will tell: How Attention Span affects Success in social Diffusion
### Part I: Images of the mathematical models
### Author: Malte Grönemann
### R version 3.6.3 (2020-02-29) on x86_64-pc-linux-gnu (Ubuntu 18.04)
################################################################################

library(ggplot2) # 3.3.2
library(deSolve) # 1.28
library(dplyr) # 1.0.0

# a colour palette suited for colourblind people
# adopted from http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette
cbPalette <- c("#0072B2", "#D55E00", "#009E73", "#56B4E9", "#F0E442", "#E69F00")
lines <- c("solid", "twodash", "longdash")


# examples of the modified bass model
df <- expand.grid(c("a: 0.2, b: 0", "a: 0.2, b: 0.6", "a: 0, b: 0.6"), # parameters
                  c(1, 3, 8), # r
                  c(0, .3, .9)) # P_act0
colnames(df) <- c("parameters", "r", "P_act0")
df$a <- rep(c(.2, .2, 0), length(df)/3)
df$b <- rep(c(0, .6, .6), length(df)/3)

mod_bass <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dp <- -b * p^2 - (a - b + 1/r) * p + a
    list(dp)
  })}

t <- seq(0, 12, by = 0.1)

temp <- ode(y = c(p = df$P_act0[1]), # approximating the curve with the ode starting from P_act0
            times = t, 
            func = mod_bass, 
            parms = c(a = df$a[1], 
                      b = df$b[1], 
                      r = df$r[1]))
df_m <- data.frame(parameters = rep(df$parameters[1], length(t)),
                   P_act0 = rep(df$P_act0[1], length(t)),
                   r = rep(df$r[1], length(t)),
                   t = temp[ , 1],
                   p = temp[ , 2])

for(c in 2:length(df$parameters)){
  temp <- ode(y = c(p = df$P_act0[c]), 
              times = t, 
              func = mod_bass, 
              parms = c(a = df$a[c], b = df$b[c], r = df$r[c]))
  temp <- data.frame(parameters = rep(df$parameters[c], length(t)),
                     P_act0 = rep(df$P_act0[c], length(t)),
                     r = rep(df$r[c], length(t)),
                     t = temp[ , 1],
                     p = temp[ , 2])
  df_m <- rbind(df_m, temp)
}
  
ggplot(df_m) +
  aes(x = t,
      y = p,
      colour = paste(P_act0 * 100, "%", sep = ""),
      linetype = paste(P_act0 * 100, "%", sep = "")) +
  geom_line() +
  facet_grid(r ~ parameters, 
             labeller = labeller(r = label_both)) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = c(0, .3, .6, .9),
                     labels = scales::percent) +
  scale_x_continuous(breaks = c(0, 4, 8, 12)) +
  scale_linetype_manual(values = lines) +
  scale_colour_manual(values = cbPalette) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "Time",
       y = "Active Population",
       colour = "Active Population at Start",
       linetype = "Active Population at Start")
ggsave("images/examples_mod_bass.pdf", height = 10, width = 15, units = "cm")


# derivative of the modified bass model
df2 <- expand.grid(c(0, .5, .9), # a
                  c(0, .5, .9), # b
                  c(1, 3, 8),  # r
                  (1:100)/100) # p
colnames(df2) <- c("a", "b", "r", "p")

df2 <- df2 %>% mutate(derivative = -b * p^2 - (a - b + 1/r) * p + a)

ggplot(df2) +
  aes(x = p,
      y = derivative,
      colour = as.character(r),
      linetype = as.character(r)) +
  geom_hline(yintercept = 0) +
  geom_line() +
  facet_grid(b ~ a, labeller = label_both) +
  scale_x_continuous(breaks = c(0, .4, .8),
                     labels = scales::percent) +
  scale_linetype_manual(values = lines) +
  scale_colour_manual(values = cbPalette) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "Active Population",
       y = "Value of the Derivative",
       colour = "Time active r: ",
       linetype = "Time active r: ")
ggsave("images/equi_mod_bass.pdf", height = 10, width = 15, units = "cm")
