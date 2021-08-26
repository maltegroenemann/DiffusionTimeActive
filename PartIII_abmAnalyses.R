################################################################################
### Time will tell: How Attention Span affects Success in social Diffusion
### Part III: Analysis of the Network ABM
### Author: Malte Grönemann
### R version 3.6.3 (2020-02-29) on x86_64-pc-linux-gnu (Ubuntu 18.04)
### NetLogo version 6.1.1, Model: DiffusionTimeActive.nlogo
################################################################################

library(ggplot2) # 3.3.2
library(dplyr)   # 1.0.0
library(margins) # 0.3.23
cbPalette <- c("#0072B2", "#D55E00", "#009E73", "#56B4E9", "#F0E442", "#E69F00")


# Data Manipulation
#####
DTA_long <- read.table("DiffusionTimeActive.csv", # change to location of data for replication
                      header = T,
                      sep = ",",
                      skip = 6) # first 6 rows only contain simulation information
colnames(DTA_long) <- c("run", "time_active", "probability", "threshold", 
                       "time", "n_links", "ave_path_length", "active")
DTA_long <- DTA_long %>% 
  arrange(run, time) %>%
  mutate(equil = 1 - 1 / time_active)


temp <- DTA_long %>%
  group_by(run) %>%
  summarise(slope = median(equil) / ifelse(min(time) == 0, NA, min(time))) # median just for summary, no variation in equil by run anyway
DTA_long <- merge(DTA_long, temp, by = "run")
remove(temp)

DTA_short <- DTA_long %>%
  group_by(run, time_active, threshold, n_links, ave_path_length, equil, slope) %>%
  summarise(died = min(active) == 0,
            reached = max(active) >= median(equil),
            slope = ifelse(min(active) == 0, NA, ))


# examples 
set.seed(83255)
DTA_long %>%
  filter(run %in% sample(max(DTA_long$run), 16)) %>%
  ggplot() +
  geom_line(aes(x = time,
                y = active,
                colour = "Observed Data",
                linetype = "Observed Data")) +
  geom_line(aes(x = time,
                y = equil,
                colour = "Analytical Equilibrium",
                linetype = "Analytical Equilibrium")) +
  facet_wrap(~run) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = c(0, .4, .8),
                     labels = scales::percent) +
  scale_x_continuous(breaks = c(0, 100, 200)) +
  scale_colour_manual(name = "", 
                      values = c("Analytical Equilibrium" = "grey", "Observed Data" = "black"), 
                      guide = "legend") +
  scale_linetype_manual(name = "", 
                        values = c("Analytical Equilibrium" = "dashed", "Observed Data" = "solid"), 
                        guide = "legend") +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.text.x = element_blank()) +
  labs(x = "Time", 
       y = "Proportion of Active Nodes")
ggsave("images/example_runs.pdf", height = 10, width = 15, units = "cm")


# survival probability
died_simple <- glm(data = subset(DTA_short, threshold == 1),
                   died ~ time_active + ave_path_length + n_links, 
                   family = binomial(link = "probit"))
died_complex <- glm(data = subset(DTA_short, threshold == 2),
                    died ~ time_active + ave_path_length + n_links, 
                    family = binomial(link = "probit"))
# some interactions were significant but resulted in worse model fit and not easily interpretable, therefore skipped

rbind(data.frame(cplot(died_simple, "time_active", draw = FALSE), 
                 threshold = "Threshold: 1"),
      data.frame(cplot(died_complex, "time_active", draw = FALSE), 
                 threshold = "Threshold: 2")) %>%
  ggplot() +
  aes(x = xvals,
      y = yvals) +
  geom_line(aes(colour = threshold,
                linetype = threshold)) +
  geom_ribbon(aes(ymin = lower,
                  ymax = upper,
                  fill = threshold),
              alpha = .2) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, by = .2)) +
  scale_x_continuous(breaks = seq(1, 10, by = 2)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_color_manual(values = cbPalette) +
  scale_fill_manual(values = cbPalette) +
  labs(x = "Average Time active",
       y = "Predicted Probability of Dying",
       fill = "", colour = "", linetype = "")
ggsave("images/died.pdf")


# assesing spreading speed
speed_simple <- summary(lm(data = DTA_simple,
                           slope ~ ave_path_length + n_links + attention_span))
speed_complex <- summary(lm(data = DTA_complex,
                            slope ~ ave_path_length + n_links + attention_span))
speed <- rbind(data.frame(var = rownames(speed_simple$coef), 
                          coef = speed_simple$coef[, 1], 
                          se = speed_simple$coef[, 2], 
                          model = "Threshold: 1"),
               data.frame(var = rownames(speed_complex$coef), 
                          coef = speed_complex$coef[, 1], 
                          se = speed_complex$coef[, 2], 
                          model = "Threshold: 2"))
speed <- speed %>% filter(var != "(Intercept)")
speed$var <- as.character(speed$var)
speed$var[speed$var == "attention_span"] <- "Attention Span"
speed$var[speed$var == "ave_path_length"] <- "Average Path Length"
speed$var[speed$var == "n_links"] <- "Number of Links"

ggplot(speed) +
  aes(x = var,
      y = coef,
      colour = model) +
  geom_point(position = position_dodge(width = 1/2)) +
  geom_linerange(aes(x = var, 
                     ymin = coef - se * (-qnorm((1-0.95)/2)), 
                     ymax = coef + se * (-qnorm((1-0.95)/2))),
                 lwd = 1, 
                 position = position_dodge(width = 1/2)) +
  geom_hline(yintercept = 0, 
             linetype = "dashed", 
             colour = "darkgrey") +
  coord_flip() +
  theme(legend.position = "bottom") +
  scale_color_manual(values = cbPalette) +
  labs(x = "", y = "", colour = "")
ggsave("images/slopes_coef.pdf")

ggplot(DTA_short) + 
  aes(x = attention_span, 
      y = slope, 
      colour = ave_path_length) + 
  geom_jitter(size = .2) +
  facet_grid(~threshold, labeller = label_both) +
  scale_x_continuous(breaks = 1:10) +
  theme(legend.position = "bottom") +
  labs(x = "Attention Span",
       y = "Average Slope",
       colour = "Average Path Length")
ggsave("images/slopes_scatter.pdf")
