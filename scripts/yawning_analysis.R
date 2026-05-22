library(lubridate)
library(tidyr)
library(ggplot2)
library(gsheet)
library(tictoc)
library(data.table)
library(dplyr)
library(here)
library(stringr)
library(gridExtra)
library(ggpmisc)
library(ggpubr)
library(ggh4x)
library(ggnewscale)
library(patchwork)
library(viridis)
library(spatstat)
library(broom)

source(here("scripts/yawning_functions.R"))
source(here("scripts/rest_functions.R"))

yawn_meta <- read.csv(here("yawn_meta.csv"))

url <- 'https://docs.google.com/spreadsheets/d/191cjK--bc0y1X6aWAW72KJKJ5LRaW1TpnLhkbqPb6fk/edit?usp=sharing'
yawns <- read.csv(text=gsheet2text(url, format='csv'), stringsAsFactors=FALSE)
yawns <- dplyr::select(yawns, 1:6, 10)
yawns <- na.omit(yawns)

yawns$sample_id <- paste0("FISH", yawns$sample_id)

samples <- unique(yawns$sample_id)

# load in rest data
als_files <- lapply(samples, function(x) list.files(path = "../als_hmm_files/", recursive = TRUE, pattern = x))

als_data <- lapply(als_files, function(x) read.csv(paste0("../als_hmm_files/", x)))

bout_data <- lapply(als_data, function(x) boutStructure(x, min_bout = 0))
bout_data <- Reduce(rbind, bout_data)
bout_data <- mutate(bout_data, start = as.POSIXct(bout_data$start, origin = "1970-01-01 00:00:00", format = "%Y-%m-%d %H:%M:%S"),
                    end = as.POSIXct(bout_data$end, origin = "1970-01-01 00:00:00", format = "%Y-%m-%d %H:%M:%S"))

als_data <- Reduce(rbind, als_data)

sleep_data <- dplyr::select(als_data, c(1:3, 7:16))
sleep_data$datetime <- as.POSIXct(sleep_data$datetime, origin = "1970-01-01 00:00:00", format = "%Y-%m-%d %H:%M:%S")

sleep_data <- filter(sleep_data, between(datetime, as.POSIXlt.character("1970-01-02 00:00:00", format = "%Y-%m-%d %H:%M:%S"), as.POSIXlt.character("1970-01-05 00:00:00", format = "%Y-%m-%d %H:%M:%S") ))
sleep_data <- full_join(sleep_data, yawn_meta)
sleep_data$X <- 1:1814205

sleep_data$phase <- ifelse(sleep_data$phase == "night", paste0("N", substr(sleep_data$phase, 2, nchar(sleep_data$phase))),
                           paste0("D", substr(sleep_data$phase, 2, nchar(sleep_data$phase))))

sleep_data$phase <- factor(sleep_data$phase, levels = c("Night","Dawn","Day", "Dusk")) 

neighbours <- list(
  "FISH20200923_c1_r0" = c("FISH20200923_c2_r0"),
  "FISH20200923_c2_r0" = c("FISH20200923_c1_r0"),
  "FISH20200923_c3_r0" = c("FISH20200923_c3_r1"),
  "FISH20200923_c3_r1" = c("FISH20200923_c3_r0", "FISH20200923_c4_r0"),
  "FISH20200923_c4_r0" = c("FISH20200923_c3_r1"),
  "FISH20200923_c5_r0" = c("FISH20200923_c6_r0"),
  "FISH20200923_c6_r0" = c("FISH20200923_c5_r0")
)

#### PREPARE DATA ####
merged_yawns <- prep_data(yawns, sleep_data)

merged_yawns <- filter(merged_yawns, state == "active")
merged_yawns$sample_id_short <- substr(merged_yawns$sample_id, 14, 18)

merged_yawns$neighbour_yawn_time <- as.POSIXct(merged_yawns$neighbour_yawn_time, origin = "1970-01-01 00:00:00", format = "%Y-%m-%d %H:%M:%S") 

# GENERAL STATS
sd(table(merged_yawns$sample_id))
var(table(merged_yawns$sample_id))
mean(table(merged_yawns$sample_id))

# set colours for plotting
sex_colours <- c("#4261eb", "#d28aff") 
sex_facets <- c("#7189ef", "#dda7ff")
sex_labels <- c("#203075", "#6d4685")
phase_colours <- c("#000093","#FFAB00", "#0080FF", "#9A2C6D")
phase_facets <- c("#7f7fc9", "#ffd57f", "#7fbfff", "#cc95b6")

# FIGURE 2: Yawning overview + sex differences

yawn_mf <- merged_yawns %>% group_by(sample_id, length_mm, sex) %>% summarise(total_yawns = n(), bout_length = mean(bout_length))
yawn_mf_phase <- merged_yawns %>% group_by(sample_id, length_mm, sex, phase, state) %>% summarise(n_yawns = n(), bout_length = mean(bout_length))

yawn_mf_durations <- filter(merged_yawns, duration > 0) %>% group_by(sample_id, length_mm, sex) %>% summarise(duration = mean(duration))
yawn_mf_durations_phase <- filter(merged_yawns, duration > 0) %>% group_by(sample_id, length_mm, sex, phase) %>% summarise(duration = mean(duration))

yawn_mf_phase$rate <- ifelse(yawn_mf_phase$phase == "Dawn" | yawn_mf_phase$phase == "Dusk", yawn_mf_phase$n_yawns/4, yawn_mf_phase$n_yawns/20)

day <- as.POSIXct('1970-01-02 08:00:00', format = "%Y-%m-%d %H:%M:%S", origin = "1970-01-01 00:00:00", tz = "EST")
night <- as.POSIXct('1970-01-02 20:00:00', format = "%Y-%m-%d %H:%M:%S", origin = "1970-01-01 00:00:00", tz = "EST")

dawn_start <-  as.POSIXct('1970-01-02 06:00:00', format = "%Y-%m-%d %H:%M:%S", origin = "1970-01-01 00:00:00", tz = "EST")
dawn_end <- as.POSIXct('1970-01-02 08:00:00', format = "%Y-%m-%d %H:%M:%S", origin = "1970-01-01 00:00:00", tz = "EST")

dusk_start <- as.POSIXct('1970-01-02 18:00:00', format = "%Y-%m-%d %H:%M:%S", origin = "1970-01-01 00:00:00", tz = "EST")
dusk_end <- format(night + hours(1), format="%H:%M:%S")

fig2a <- ggplot() +
  annotate("rect", xmin = night - hours(2), xmax = dawn_start + hours(2), ymin =-Inf, ymax=0.5, fill = "#0080FF", alpha = 0.75) +
  annotate("rect", xmin = night + days(1) - hours(2), xmax = dawn_start + days(1) + hours(2), ymin =-Inf, ymax=0.5, fill = "#0080FF", alpha = 0.75) +
  annotate("rect", xmin = -Inf, xmax = dawn_start, ymin =-Inf, ymax=0.5, fill = "#000093", alpha = 0.75) +
  annotate("rect", xmin = night + days(1), xmax =Inf, ymin =-Inf, ymax=0.5, fill = "#000093", alpha = 0.75) +
  annotate("rect", xmin = dusk_start + hours(2), xmax =dawn_start + days(1), ymin =-Inf, ymax=0.5, fill = "#000093", alpha = 0.75) +
  annotate("rect", xmin = dawn_start, xmax = day, ymin =-Inf, ymax=0.5, fill = "#FFAB00", alpha = 0.75) +
  annotate("rect", xmin = dawn_start + days(1), xmax = day + days(1), ymin =-Inf, ymax=0.5, fill = "#FFAB00", alpha = 0.75) +
  annotate("rect", xmin = dusk_start, xmax = night, ymin =-Inf, ymax=0.5, fill = "#9A2C6D", alpha = 0.75) +
  annotate("rect", xmin = dusk_start + days(1), xmax = night + days(1), ymin =-Inf, ymax=0.5, fill = "#9A2C6D", alpha = 0.75) +
  geom_tile(data = merged_yawns, aes(x = datetime, y = sample_id_short, colour = sex), width = 0.1, position = position_dodge(0.8)) +
  scale_x_datetime(expand = c(0, 0)) + 
  scale_y_discrete(expand = expansion(add = c(0.75, 0))) + 
  scale_colour_manual(values = c("#4261eb", "#d28aff")) +
  guides(colour = "none") +
  labs(x = "Time", y = "Sample ID", title = expression(paste(italic("B. microlepis")," yawning across time"))) +
  theme_pubr() +
  theme(text = element_text(size = 16))

fig2b <- basic_box(yawn_mf, sex, total_yawns, sex, strip = NA, scales = "fixed",
                   fill_colours = sex_colours, label_colours = sex_labels, 
                   facet_phase = F, facet_sex = F,  
                   "Yawns by sex", "Sex", "Number of yawns per individual", include_points = T)

fig2c <- basic_box(yawn_mf_durations, sex, duration, sex, strip = NA, scales = "fixed",
                   fill_colours = sex_colours, label_colours = sex_labels,
                   facet_phase = F, facet_sex = F,  
                   "Yawn duration by sex", "Sex", "Duration (sec)", T)


fig2d <- basic_box(yawn_mf_phase, phase, rate, phase, strip = sex_facets, 
                   fill_colours = phase_colours, label_colours = NA,
                   facet_phase = F, facet_sex = T, scales = "fixed",
                   "Rate of yawning across time", "Phase", "Rate (yawns per hour)", T)
design <- "AAB
CDD"

fig2 <- fig2a + fig2b + fig2c + fig2d  + plot_layout(axis_titles = "collect", design = design) 

#ggsave(filename = "fig2.pdf", path = "plots/yawning/", plot = fig2, width = 8.5*1.5, height = 6*1.5, units = "in", device='pdf', dpi=700)

# FIGURE 3 & S1-2: Yawning and state changes

facets3a <- c("#7189ef", "#dda7ff","#7f7fc9", "#7f7fc9")
fig3a <- plot_mc_analysis(filter(merged_yawns, phase == "Night"), "Pre-sleep", facets3a, 30, T, F, T)

facetsS1 <- c("#7189ef", "#7189ef", "#7189ef","#7189ef","#dda7ff", "#dda7ff", "#dda7ff", "#dda7ff",
              "#7f7fc9","#ffd57f", "#7fbfff", "#cc95b6", "#7f7fc9", "#ffd57f", "#7fbfff", "#cc95b6")
figS1 <- plot_mc_analysis(filter(merged_yawns), "Pre-sleep", facetsS1, 30, T, F, T)

strip <- strip_themed(background_x = elem_list_rect(fill = phase_facets))
fig3b <- ggplot(filter(merged_yawns, phase == "Night"), aes(x = sex, y = next_bout)) +
  geom_violin(aes(fill = sex),width = 1) +
  scale_fill_manual(values = sex_colours) +
  facet_wrap2(~ phase, scales = "free", strip = strip) +
  stat_summary( aes(colour = sex, label = round(after_stat(y), 0)), fun=mean, geom="text", position = position_dodge(0.25)) +
  scale_x_discrete(limits = c("Male", "Female"))+
  guides(fill = "none", colour = "none") +
  scale_colour_manual(values = sex_labels) +
  stat_pwc(method = "t_test", label = "p = {p}, t = {round(statistic, 2)}, df = {round(df, 2)}", remove.bracket = TRUE) +
  labs(x = "Sex", y = "Time until next rest bout (sec)", title = paste0("Time between yawns and state changes")) +
  theme_bw() +
  theme(panel.grid = element_blank(), text = element_text(size = 16))

# summarize sleep info
sleep_stats_phase <- sleep_data %>%
  group_by(sample_id, sex, length_mm, rest, phase) %>% 
  summarise(count = n())

sleep_stats_phase$hours <- sleep_stats_phase$count/3600
sleep_stats_phase <- filter(sleep_stats_phase, rest == T)
sleep_stats_phase$sex <- ifelse(sleep_stats_phase$sex == "f", "Female", "Male")
sleep_stats_phase$sex <- factor(sleep_stats_phase$sex, levels = c("Male","Female")) 

fig3c <- basic_box(sleep_stats_phase, sex, hours, sex,  strip = phase_facets, 
                   fill_colours = sex_colours, label_colours = sex_labels,
                   facet_sex = F, facet_phase = T, scales = "free", 
                   title = "Phase-specific sex differences in sleep", "Sex", "Total sleep (hours)", T)

# correlations
yawn_mf_phase <- merge(yawn_mf_phase, sleep_stats_phase, by = c("sample_id", "phase", "sex", "length_mm"))

corr <- get_correlations(yawn_mf_phase, c("length_mm", "n_yawns", "hours"), group_sex = TRUE)
dur_corr <- get_correlations(filter(merged_yawns, duration > 0), c("duration", "next_bout", "mean_y_nt"), TRUE)

fullgrid_strip <- c(rep("#7189ef", 4), rep("#dda7ff", 4), rep(phase_facets,2))

fig3d <- corr_plot(filter(corr, phase == "Night"), comparison = "n", facets3a)

figS2a <- corr_plot(corr, comparison = "n", fullgrid_strip, nrow = 2)
 
facet_colours <- c("#7189ef", "#dda7ff","#ffd57f", "#ffd57f")
fig3e <- corr_plot(filter(dur_corr, phase == "Dawn"), comparison = "duration", facet_colours)

figS2b <- corr_plot(dur_corr, comparison = "duration", fullgrid_strip, nrow = 2)
 
design <- "AB
CC"

fig3upper <- fig3a + fig3b + fig3c + plot_layout(axis_titles = "collect", design = design)
fig3lower <- fig3d + fig3e + plot_layout(axis_titles = "collect")
figS2 <- figS2a + figS2b + plot_layout(nrow = 2)

#ggsave(filename = "fig3upper.pdf", path = "plots/yawning/", plot = fig3upper, width = 8.5*1.5, height = 6*1.5, units = "in", device='pdf', dpi=700)
#ggsave(filename = "fig3lower.pdf", path = "plots/yawning/", plot = fig3lower, width = 8.5*1.5, height = 3*1.5, units = "in", device='pdf', dpi=700)
#ggsave(filename = "figS2.pdf", path = "plots/yawning/", plot = figS2, width = 8.5*1.5, height = 11*1.5, units = "in", device='pdf', dpi=700)
#ggsave(filename = "figS1.pdf", path = "plots/yawning/", plot = figS1, width = 8.5*1.5, height = 6*1.5, units = "in", device='pdf', dpi=700)

# FIGURE 4 & S3-4: Social yawning

y_pos <- group_by(merged_yawns, sample_id, phase) %>%
  summarise(count = n(), y_pos = mean(mean_y_nt), x_pos = mean(mean_x_nt), next_bout = mean(next_bout), prev_bout = mean(prev_bout))

fig5a <- basic_box(y_pos, phase, y_pos, phase, NA, phase_colours, NA,
                  facet_sex = F, facet_phase = F, scales = "free",
                  "y-position of yawns", "Phase", "y-position", T)

left_neighbours <- c("FISH20200923_c1_r0", "FISH20200923_c3_r0", "FISH20200923_c3_r1", "FISH20200923_c5_r0")
fig5b <- plot_densities(merged_yawns, T, to_plot = "Yawns")

fig5c <- plot_mc_analysis(filter(merged_yawns, phase == "Day"), "Social", c("#7fbfff"), 30, F, F, T)

fig5d <- basic_box(filter(merged_yawns, phase == "Day",), sex, neighbour_yawn_timediff, sex,  c("#7fbfff"),
                   sex_colours, sex_labels, facet_sex = F, facet_phase = T, scales = "free",
                   "Yawn response times", "Sex", "Time from neighbour yawn (sec)")

merged_yawns$distance <- ifelse(merged_yawns$neighbour_id %in% left_neighbours, 
                                sqrt(((merged_yawns$mean_x_nt + 1) - merged_yawns$neighbour_xpos)^2 + (merged_yawns$mean_y_nt - merged_yawns$neighbour_ypos)^2), 
                                sqrt((merged_yawns$mean_x_nt - (merged_yawns$neighbour_xpos + 1))^2 + (merged_yawns$mean_y_nt - merged_yawns$neighbour_ypos)^2))

fig5e <- basic_box(filter(merged_yawns, phase == "Day"), sex, distance, sex,  c("#7fbfff"), 
                   sex_colours, sex_labels, F, T, "free",
                   "Neighbour proximity", "Sex", "Distance from neighbour")


design <- 
"ABB
CDE"

fig5 <- fig5a + fig5b + fig5c + fig5d + fig5e + plot_layout(guides = 'collect', axis_titles = 'collect', design = design) 

#ggsave(filename = "locations.png", path = "plots/yawning/", plot = fig5b, width = 8.5*1.5, height = 3*1.5, units = "in", device='png', dpi=700)
#ggsave(filename = "fig5.pdf", path = "plots/yawning/", plot = fig5, width = 8.5*1.5, height = 6*1.5, units = "in", device='pdf', dpi=700)

sleeppos <- plot_densities(sleep_data, F, "Sleep")
yawnpos <- plot_densities(merged_yawns, F, "Yawns")

figS3 <- sleeppos + yawnpos + plot_layout(axis_titles = "collect", guides = "collect")

#ggsave(filename = "figS3.png", path = "plots/yawning/", plot = figS3, width = 8.5*1.5, height = 6*1.5, units = "in", device='png', dpi=700)

figS4a <- plot_mc_analysis(filter(merged_yawns), "Social", phase_facets, 30, F, F, T)

merged_yawns$trigger <- ifelse(merged_yawns$datetime == merged_yawns$neighbour_yawn_time, "both", ifelse(merged_yawns$datetime <= merged_yawns$neighbour_yawn_time, merged_yawns$sample_id, merged_yawns$neighbour_id))
merged_yawns$trigger <- as.character(merged_yawns$trigger)

merged_yawns$socialyawn_type <- ifelse(merged_yawns$sample_id == merged_yawns$trigger, "Origin", "Response") 
yawn_triggers <- group_by(merged_yawns, sample_id, sex, socialyawn_type, phase) %>% summarise(count = n())

figS4b <- ggplot(filter(yawn_triggers, phase == "Day"), aes(x = sex, y = count, fill = sex)) + 
  geom_boxplot()+
  geom_point(size = 1.5) +
  facet_wrap2(~socialyawn_type, strip = strip_themed(background_x = elem_list_rect(fill = c("#7fbfff")))) +
  theme_pubr() +
  stat_summary(aes(colour = sex, label = round(after_stat(y), 2)), fun=mean, geom="text") +
  scale_fill_manual(values = sex_colours) +
  scale_colour_manual(values = sex_labels) +
  guides(fill = "none", colour = "none") +
  labs(x = "Sex", y = "Number of yawns", title = "Yawn contagion by sex") +
  stat_pwc(method = "t_test", label = "p = {p}, t = {round(statistic, 2)}, df = {round(df, 2)}", remove.bracket = TRUE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 16))

merged_yawns$mc_distance <- ifelse(merged_yawns$neighbour_id %in% left_neighbours, 
                                   sqrt(((merged_yawns$matched_xpos + 1) - merged_yawns$mc_neighbour_xpos)^2 + (merged_yawns$matched_ypos - merged_yawns$mc_neighbour_ypos)^2), 
                                   sqrt((merged_yawns$matched_xpos - (merged_yawns$mc_neighbour_xpos + 1))^2 + (merged_yawns$matched_ypos - merged_yawns$mc_neighbour_ypos)^2))

figS4c <- basic_box(filter(merged_yawns, phase == "Day"), sex, mc_distance, sex, c("#7fbfff","#7fbfff"),
                    sex_colours, sex_labels, facet_sex = F, facet_phase = T, scales = "free",
                    "MC neighbour proximity", "Sex", "Distance from neighbour")
design = "AAAA
BBBC"

figS4 <- figS4a + figS4b + figS4c + plot_layout(axis_titles = "collect", design = design)

#ggsave(filename = "figS4.pdf", path = "plots/yawning/", plot = figS4, width = 8.5*1.5, height = 6*1.5, units = "in", device='pdf', dpi=700)
