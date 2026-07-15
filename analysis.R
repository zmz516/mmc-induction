library(tidyverse)
library(readxl)
library(dplyr)
library(ggplot2)


# import and tidy ---------------------------------------------------------

df <- read.csv("data-raw/MMC_induction.CSV")

# column headers give the well row
# the first row gives the well column
# the first two columns do not contain data

row_letter <- sub("\\..*", "", names(df)[-(1:2)]) 
col_number <- as.integer(df[1, -(1:2)])
well <- paste0(row_letter, col_number)

time <- df[-(1:2), 2]

readings <- as.data.frame(lapply(df[-(1:2),-(1:2)], as.numeric))
colnames(readings) <- well

readings <- readings |>
  mutate(time = time) |>
  relocate(time, .before = 1)

# plate map ---------------------------------------------------------------
plate_map <- tibble(row = rep(LETTERS[1:6], each = 12),
                    col = rep(1:12, times = 6),
                    mmc = rep(c(0, 0.1, 0.5, 1.0, 1.5, 3.0), each = 12),
                    strain = rep(rep(c("lys", "host"), each = 6), times = 6)) |> 
  mutate(well = paste0(row, col),
         rep = rep(rep(1:6, times = 2), times = 6)) # 6 reps per strain/concentration


# pivot longer ------------------------------------------------------------
longer <- readings |> 
  pivot_longer(cols = -time,       # which columns to collapse (everything except time)
               names_to = "well",   # name of the new header column
               values_to = "OD") # name of the new data value column


# convert time ------------------------------------------------------------
hours <- as.numeric(str_extract(longer$time,  "\\d+(?=\\s*h)"))
minutes <- as.numeric(str_extract(longer$time, "\\d+(?=\\s*min)"))
hours[is.na(hours)] <- 0
minutes[is.na(minutes)] <- 0
longer$time_mins <- (hours * 60) + minutes

# blank subtract ----------------------------------------------------------
  # wells in row G contained the LB control
  # so use these to subtract blank from OD values
blank <- longer |> 
  filter(str_starts(well, "G")) |> 
  group_by(time_mins) |> 
  summarise(blank_OD = mean(OD, na.rm = T), .groups = "drop")

# link the metadata -------------------------------------------------------
data <- longer |> inner_join(plate_map |> select(well, strain, mmc, rep), by = "well") |> 
  left_join(blank, by = "time_mins") |> 
  mutate(OD_corr = OD - blank_OD,
         mmc_fact = factor(mmc, levels = c(0, 0.1, 0.5, 1.0, 1.5, 3.0),
                           labels = paste0(c(0, 0.1, 0.5, 1.0, 1.5, 3.0), " ug/mL")))
        
# summarizing the data ----------------------------------------------------
summary_data <- data |> 
  group_by(strain, mmc_fact, time_mins) |> 
  summarise(mean_OD = mean(OD, na.rm = T), .groups = "drop")

# visualizing the data ----------------------------------------------------
custom_titles <- c("host" = "Host (009)",
                   "lys"  = "Lysogen")

avg_plot <- ggplot(summary_data, aes(x = time_mins, y = mean_OD, colour = mmc_fact)) +
  geom_line()+
  facet_wrap(~ strain, labeller = as_labeller(custom_titles))+
  labs(x = "Time (mins)",
       y = "OD 600nm",
       colour = "[MMC]") +
  theme_bw()
avg_plot
ggsave("plots/avg_plot.png", plot = avg_plot, width = 7, height = 5, dpi = 300)

# plotting the raw data for each concentration --------------------------------
growth_by_mmc <- ggplot(data, aes(x = time_mins, y = OD_corr, colour = strain, group = well))+
  geom_line()+
  theme_bw() +
  facet_wrap(~ mmc_fact,
             ncol = 3, nrow = 2) +
  labs(title = "Growth by MMC Concentration",
       x = "Time (mins)",
       y = "OD 600nm",
       colour = "Strain")+
  scale_colour_discrete(labels = c("host" = "Host (009)", 
                                 "lys" = "Lysogen"))
growth_by_mmc
ggsave("plots/growth_by_mmc.png", plot = growth_by_mmc, width = 7, height = 5, dpi = 300)


#  can see that one of the wells containing 0.0 MMC is an outlier
#  plot to identify which well

well_plot <- data |> filter(strain == "lys") |> filter(str_starts(well, "A")) |> 
  ggplot(aes(x = time_mins, y = OD_corr, colour = well, group = well))+ geom_line()+
  theme_bw()+
  labs(title = "Well A5 is causing problems!",
       x = "Time (mins)",
       y = "OD 600nm",
       colour = "Well")
well_plot
ggsave("plots/well_plot.png", plot = well_plot, width = 7, height = 5, dpi = 300)

# area under the curve ----------------------------------------------------

# remove well A5 from the dataset
# caluclate area under the curve to compare growth

remove_a5 <- data |> filter(well !="A5")

library(MESS)
auc_data <- remove_a5 |> 
  distinct(strain, mmc, mmc_fact, well, rep, time_mins, .keep_all = TRUE) |>  # De-duplicate identical time points within the same well/rep combo
  group_by(strain, mmc, mmc_fact, well, rep) |> 
  summarise(AUC = auc(time_mins, OD_corr, type = "spline"), 
    .groups = "drop")

auc_plot <- ggplot(data = auc_data, aes(x = mmc_fact, y = AUC, fill = strain))+
  geom_boxplot()+
  theme_bw()+
  labs(title = "AUC by MMC (ug/mL)",
       fill = "Strain")+
  scale_fill_discrete(labels = c("host" = "Host (009)", 
                                   "lys" = "Lysogen"))
auc_plot
ggsave("plots/auc_plot.png", plot = auc_plot, width = 7, height = 5, dpi = 300)

# t-test ------------------------------------------------------------------

# test the significance of any apparent differences in growth

library(dplyr)
library(rstatix)

t_test <- auc_data |> 
  group_by(mmc_fact) |> 
  t_test(AUC~strain) |> 
  adjust_pvalue(method = "BH") |>  # FDR adjustment for multiple testing
  add_significance() |> 
  add_x_position(x = "mmc_fact") # IMPORTANT: This tells ggplot to place the bracket over the x-axis groups


# add t-test results to plot ----------------------------------------------
library(ggpubr)

auc_sig_plot <- ggplot(data = auc_data, aes(x = mmc_fact, y = AUC, fill = strain))+
  geom_boxplot()+
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), size = 1, shape = 1)+ 
  theme_bw()+
  labs(title = "AUC by MMC (ug/mL)",
       fill = "Strain")+
  scale_fill_discrete(labels = c("host" = "Host (009)", 
                                 "lys" = "Lysogen")) +
  stat_pvalue_manual(t_test, 
                     label = "p.adj.signif", # displays *, **, or ns
                     y.position = max(auc_data$AUC) * 1.05) + # bracket sits just above the highest data point
  labs(
    title = "AUC by MMC (ug/mL)",
    x = "MMC Concentration",
    y = "AUC (OD 600nm)")
auc_sig_plot
ggsave("plots/auc_sig_plot.png", plot = auc_sig_plot, width = 7, height = 5, dpi = 300)

