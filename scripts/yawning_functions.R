prep_data <- function(yawns, sleep_data) {
  yawns$day <- ifelse(as.numeric(substr(yawns$sample_id, 5, 12)) - as.numeric(str_replace_all(yawns$Date.of.video, "/", "")) == 0,
                      1, as.numeric(str_replace_all(yawns$Date.of.video, "/", "")) - as.numeric(substr(yawns$sample_id, 5, 12)) + 1)
  
  yawns$Video.start.time <- paste0("1970-01-0", yawns$day, " ", yawns$Video.start.time)
  
  yawns$Video.start.time <- as.POSIXct(yawns$Video.start.time, origin = "1970-01-01 00:00:00", format = "%Y-%m-%d %H:%M:%S")
  
  yawns$Actual.time. <- paste0("1970-01-0", yawns$day, " ", yawns$Actual.time.)
  yawns$Actual.time. <- as.POSIXct(yawns$Actual.time., origin = "1970-01-01 00:00:00", format = "%Y-%m-%d %H:%M:%S")
  
  yawns$day <- ifelse(hour(yawns$Video.start.time) == 23 & hour(yawns$Actual.time.) == 0, yawns$day + day(1), yawns$day)
  
  yawns$datetime <- paste0("1970-01-0", yawns$day, substr(yawns$Actual.time., 11,19))
  yawns$datetime <- as.POSIXct(yawns$datetime, origin = "1970-01-01 00:00:00", format = "%Y-%m-%d %H:%M:%S")
  
  yawns <- filter(yawns, between(datetime, as.POSIXlt.character("1970-01-02 00:00:00", format = "%Y-%m-%d %H:%M:%S"), as.POSIXlt.character("1970-01-04 00:00:00", format = "%Y-%m-%d %H:%M:%S") ))
  
  yawns <- filter(yawns, Timestamp != "N/A")
  merged_yawns <- full_join(yawns, sleep_data, by = c("sample_id", "datetime"))
  merged_yawns <- na.omit(merged_yawns)
  
  merged_yawns <- filter(merged_yawns, between(datetime, as.POSIXlt.character("1970-01-02 00:00:00", format = "%Y-%m-%d %H:%M:%S"), as.POSIXlt.character("1970-01-04 00:00:00", format = "%Y-%m-%d %H:%M:%S") ))
  
  merged_yawns$sex <- ifelse(merged_yawns$sex == "f", "Female", "Male")
  merged_yawns$sex <- factor(merged_yawns$sex, levels = c("Male","Female")) 
  
  merged_yawns$hour <- hour(merged_yawns$datetime)
  
  # SLEEP-WAKE
  index <- apply(merged_yawns, 1, function(x) which(bout_data$sample_id == x["sample_id"] & bout_data$start <= x["datetime"] & bout_data$end >= x["datetime"]))
  
  merged_yawns$index <- unlist(lapply(index, function(x) x[1]))
  
  merged_yawns <- mutate(merged_yawns, state = bout_data$state[index], bout_start = bout_data$start[index], bout_end = bout_data$end[index], 
                         bout_length = bout_data$length[index], next_bout_length = bout_data$length[index + 1])
  
  merged_yawns$prev_bout <- as.numeric(abs(difftime(merged_yawns$bout_start, merged_yawns$datetime)))
  merged_yawns$next_bout <- as.numeric(abs(difftime(merged_yawns$bout_end, merged_yawns$datetime)))
  
  merged_yawns$matched_time <- merged_yawns$datetime + days(1)
  
  matched_index <- apply(merged_yawns, 1, function(x) which(bout_data$sample_id == x["sample_id"] & bout_data$start <= x["matched_time"] & bout_data$end >= x["matched_time"]))
  pos_index <- apply(merged_yawns, 1, function(x) which(sleep_data$sample_id == x["sample_id"] & sleep_data$datetime == x["matched_time"]))
  
  merged_yawns$matched_index <- unlist(lapply(matched_index, function(x) x[1]))
  merged_yawns$pos_index <- unlist(lapply(pos_index, function(x) x[1]))
  
  merged_yawns <- mutate(merged_yawns, matched_bout_start = bout_data$start[matched_index], matched_bout_end = bout_data$end[matched_index])
  merged_yawns <- mutate(merged_yawns, matched_xpos = sleep_data$mean_x_nt[pos_index], matched_ypos = sleep_data$mean_y_nt[pos_index])
  
  merged_yawns$matched_prev_bout <- abs(difftime(merged_yawns$matched_bout_start, merged_yawns$matched_time))
  merged_yawns$matched_next_bout <- abs(difftime(merged_yawns$matched_bout_end, merged_yawns$matched_time))
  
  # SOCIAL YAWNING
  for (i in 1:nrow(merged_yawns)) {
    sample <- merged_yawns$sample_id[i]
    yawn <- merged_yawns$datetime[i]
    if (merged_yawns$day[i] == 2) { matched_yawn <- merged_yawns$datetime[i] + days(1) }
    if (merged_yawns$day[i] == 3) { matched_yawn <- merged_yawns$datetime[i] - days(1) }
    
    neighbour <- neighbours[[sample]]
    nonneighbour <- sample(samples[!(samples %in% c(neighbours[[sample]], sample))], 1)
    
    neighbour_yawns <- merged_yawns[merged_yawns$sample_id %in% neighbour, ]
    nonneighbour_yawns <- merged_yawns[merged_yawns$sample_id %in% nonneighbour, ]
    
    diff <- abs(difftime(neighbour_yawns$datetime, yawn, units = "sec"))
    nn_diff <- abs(difftime(nonneighbour_yawns$datetime, yawn, units = "sec"))
    matched_diff <- abs(difftime(neighbour_yawns$datetime, matched_yawn, units = "sec"))
    
    merged_yawns$neighbour_yawn_timediff[i] <- min(diff)
    merged_yawns$nn_yawn_timediff[i] <- min(nn_diff)
    merged_yawns$mc_neighbour_yawn_timediff[i] <- min(matched_diff)
    merged_yawns$neighbour_id[i] <- neighbour_yawns$sample_id[which(diff == min(diff))]
    merged_yawns$neighbour_yawn_time[i] <- as.POSIXct(neighbour_yawns$datetime[which(diff == min(diff))], origin = "1970-01-01 00:00:00", format = "%Y-%m-%d %H:%M:%S") 
    merged_yawns$neighbour_xpos[i] <- neighbour_yawns$mean_x_nt[which(diff == min(diff))]
    merged_yawns$neighbour_ypos[i] <- neighbour_yawns$mean_y_nt[which(diff == min(diff))]
    merged_yawns$neighbour_state[i] <- neighbour_yawns$rest[which(diff == min(diff))]
    merged_yawns$neighbour_state[i] <- ifelse(merged_yawns$neighbour_state[i] == FALSE, "active", "rest")
    
    merged_yawns$mc_neighbour_xpos[i] <- neighbour_yawns$mean_x_nt[which(matched_diff == min(matched_diff))]
    merged_yawns$mc_neighbour_ypos[i] <- neighbour_yawns$mean_y_nt[which(matched_diff == min(matched_diff))]
    
  }
  return(merged_yawns)
}

get_correlations <- function(data, variables, group_sex) {
  grid <- expand.grid(variables, variables)
  
  if (group_sex) {
    data <- data %>%
      group_by(sex, phase) 
    n <- 2
  }
  else {
    data <- data %>%
      group_by(phase) 
    n <- 1
  }
  
  output.grid <- apply(grid, 1, function(x) {
    if (x[1] == x[2]) {
      return (list("estimate" = rep(NA, 4*n), "p.value" = rep(1, 4*n)))
    } 
    else {
      
      corr <- data %>%
        reframe(tidy(cor.test(eval(as.name(paste(x[1]))), eval(as.name(paste(x[2]))), method = "pearson")))
      
      return(corr)
    }
    
  })
  
  if (group_sex) { full_grid <- expand.grid(c("Night", "Dawn", "Day", "Dusk"), c("Male", "Female"), variables, variables) }
  else { full_grid <- expand.grid(c("Night", "Dawn", "Day", "Dusk"), variables, variables) }
  
  rs <- (lapply(output.grid, function(x) x$estimate))
  ps <- (lapply(output.grid, function(x) x$p.value))
  
  rs <- unlist(rs)
  ps <- unlist(ps)
  
  out <- data.frame(namex = full_grid$Var3,
                    namey = full_grid$Var4,
                    sex = full_grid$Var2,
                    phase = full_grid$Var1,
                    r = rs,
                    p.value = ps)
  
  out$label <- ifelse(out$p.value < 0.0001, "****", ifelse(out$p.value < 0.001, "***", ifelse(out$p.value < 0.01, "**", ifelse(out$p.value < 0.05, "*", "")) ) )
  return(out)
}

### PLOTTING ###
plot_mc_analysis <- function(data, comparison = c("Pre-sleep", "Social"), strip, interval = 30, filt_sex = FALSE, sex_diff = FALSE, phase = FALSE) {
  if (comparison == "Pre-sleep") {
    data$yawn <- data$next_bout < interval
    data$nonyawn <- data$matched_next_bout < interval
  } 
  if (comparison == "Social") {
    data$yawn <- data$neighbour_yawn_timediff < interval
    data$nonyawn <- data$mc_neighbour_yawn_timediff < interval
  }
  
  values = c("grey","darkgrey")
  
  if (phase) {
    yawn_summary <- data %>% group_by(sample_id, phase, sex) %>% 
      summarise(count = n(), yawn = sum(yawn, na.rm = TRUE), nonyawn = sum(nonyawn, na.rm = TRUE))
  }
  else {
    yawn_summary <- data %>% group_by(sample_id, sex) %>% 
      summarise(count = n(), yawn = sum(yawn, na.rm = TRUE), nonyawn = sum(nonyawn, na.rm = TRUE))
  }
  
  yawn_summary <- mutate(yawn_summary, yawn_prop = yawn/count, nonyawn_prop = nonyawn/count)
  
  yawn_props <- pivot_longer(yawn_summary, cols = 4:length(yawn_summary), names_to = "condition")
  yawn_props <- filter(yawn_props, grepl("prop", condition))
  yawn_props <- mutate(yawn_props, condition = ifelse(condition == "yawn_prop", "Yawns", "Non-yawns"))

  stat_method <- "t.test"
  plot <- basic_box(yawn_props, condition, value, condition, strip, 
                    fill_colours = c("darkgrey", "grey"), label_colours = c("black","black"),
                    facet_sex = filt_sex, facet_phase = phase, scales = "free_x", 
                    title = paste0(comparison, " yawning proportions"), xlab = "Condition", ylab = "Proportion", 
                    include_points = T, paired = T)
  return(plot)
}

basic_box <- function(data, x_var, y_var, fill_var, strip, fill_colours, label_colours, facet_sex = FALSE, facet_phase = FALSE, scales = "fixed", title, xlab, ylab, include_points = FALSE, paired = F) {
  x_vals <- ungroup(data) %>% dplyr::select({{ x_var }})
  x_vals <- x_vals[[1]]
  n_x <- length(unique(x_vals))
  
  y_vals <- ungroup(data) %>% dplyr::select({{ y_var }})
  y_vals <- y_vals[[1]]
  
  strip <- strip_themed(background_x = elem_list_rect(fill = strip))
  
  stat_method <- ifelse(n_x > 2, "anova", "t.test")
  
  plot <- ggplot(data, aes(x = {{x_var}}, y = {{y_var}})) + 
    geom_boxplot(outlier.size = 0.5, aes(fill = {{fill_var}})) + 
    {if (include_points) geom_point(size = 1.5)} +
    {if (facet_phase) facet_wrap2(~ phase, strip = strip, scales = scales, ncol = 4)}+
    {if (facet_sex) facet_wrap2(~ sex, strip = strip, scales = scales, ncol = 2)}+
    {if (facet_phase & facet_sex) facet_wrap2(~ sex + phase, strip = strip, scales = scales, ncol = 4)}+
    {if (stat_method == "t.test") stat_pwc(method = "t_test", label = "p = {p}, t = {round(statistic, 2)}, df = {round(df, 2)}", remove.bracket = TRUE, method.args = list(paired = paired))} +
    {if (stat_method == "anova") stat_anova_test(label.x.npc = "left", label = "as_detailed_expression", size = 3)} +
    scale_fill_manual(values = fill_colours) +
    {if (!is.na(label_colours[1])) stat_summary(aes(colour = {{fill_var}}, label = round(after_stat(y), 2)), fun=mean, geom="text", position = position_dodge(width = 0.75))} +
    guides(fill = "none", colour = "none", alpha = "none") +
    scale_colour_manual(values = label_colours) +
    ggtitle(title) +
    labs(x = xlab, y = ylab, title = title) +
    guides(fill = "none") +
    theme_pubr() +
    theme(panel.grid = element_blank(), text = element_text(size = 16))
  
  return(plot)
}

plot_densities <- function(data, phase = T, to_plot =  c("Yawns", "Sleep")) {
  if (to_plot == "Yawns") {
    colour <- "magma"
  }
  if (to_plot == "Sleep") {
    data <- filter(data, rest == TRUE)
    colour <- "mako"
  }
  
  data$mean_x_nt <- ifelse(data$sample_id %in% left_neighbours, data$mean_x_nt, 1 - data$mean_x_nt)
  data <- filter(data, sample_id != "FISH20200923_c3_r1")
  
  if (phase) {
    night <- filter(data, phase == "Night")
    day <- filter(data, phase == "Day")
    dawn <- filter(data, phase == "Dawn")
    dusk <- filter(data, phase == "Dusk")
    
    night_x <- night$mean_x_nt
    night_y <- night$mean_y_nt
    
    day_x <- day$mean_x_nt
    day_y <- day$mean_y_nt
    
    dawn_x <- dawn$mean_x_nt
    dawn_y <- dawn$mean_y_nt
    
    dusk_x <- dusk$mean_x_nt
    dusk_y <- dusk$mean_y_nt
    
    night_pp  <- ppp(night_x, night_y, window=owin(xrange=range(0,1), yrange=range(c(0,1))))
    day_pp  <- ppp(day_x, day_y, window=owin(xrange=range(0,1), yrange=range(c(0,1))))
    dawn_pp  <- ppp(dawn_x, dawn_y, window=owin(xrange=range(0,1), yrange=range(c(0,1))))
    dusk_pp  <- ppp(dusk_x, dusk_y, window=owin(xrange=range(0,1), yrange=range(c(0,1))))
    
    night_d  <- density(night_pp, sigma=bw.diggle(night_pp))
    day_d  <- density(day_pp, sigma=bw.diggle(day_pp))
    dawn_d  <- density(dawn_pp, sigma=bw.diggle(dawn_pp))
    dusk_d  <- density(dusk_pp, sigma=bw.diggle(dusk_pp))
    
    night_df <- as.data.frame(night_d, xy=TRUE)
    colnames(night_df) <- c("x","y","density")
    night_df$norm_density <- night_df$density / max(night_df$density)
    night_df$phase <- "Night"
    
    day_df <- as.data.frame(day_d, xy=TRUE)
    colnames(day_df) <- c("x","y","density")
    day_df$norm_density <- day_df$density / max(day_df$density)
    day_df$phase <- "Day"
    
    dawn_df <- as.data.frame(dawn_d, xy=TRUE)
    colnames(dawn_df) <- c("x","y","density")
    dawn_df$norm_density <- dawn_df$density / max(dawn_df$density)
    dawn_df$phase <- "Dawn"
    
    dusk_df <- as.data.frame(dusk_d, xy=TRUE)
    colnames(dusk_df) <- c("x","y","density")
    dusk_df$norm_density <- dusk_df$density / max(dusk_df$density)
    dusk_df$phase <- "Dusk"
    
    dn_df <- rbind(night_df, day_df)
    dd_df <- rbind(dawn_df, dusk_df)
    df <- rbind(dn_df, dd_df)
    
    df$phase <- factor(df$phase, levels = c("Night","Dawn","Day", "Dusk")) 
    
    strip <- strip_themed(background_x = elem_list_rect(fill = phase_facets))
  }
  else {
    x <- data$mean_x_nt
    y <- data$mean_y_nt
    
    pp  <- ppp(x, y, window=owin(xrange=range(0,1), yrange=range(c(0,1))))
    d  <- density(pp, sigma=bw.diggle(pp))
    
    df <- as.data.frame(d, xy=TRUE)
    colnames(df) <- c("x","y","density")
    
    df$norm_density <- df$density / max(df$density)
  }
  
  plot <- ggplot(df, aes(x = x, y = y, fill = norm_density)) +
    theme_bw() +
    geom_tile(colour = NA) +
    #guides(fill = "none") +
    scale_fill_viridis(name="Density", option=colour) +
    {if (phase) facet_wrap2(~ phase, ncol = 4, strip = strip)} +
    labs(title = to_plot, x = "x coordinate", y = "y coordinate") +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0)) +
    theme(line = element_blank(), panel.border = element_blank(), rect = element_blank(), 
          axis.text = element_blank(), text = element_text(size = 18))
  
  plot
  return(plot)
}

corr_plot <- function(data, comparison = c("n", "duration"), strip, nrow = 1) {
  strip <- strip_themed(background_x = elem_list_rect(fill = strip))
  
  title <- ifelse(comparison == "n", "Number of yawns correlations", "Yawn duration correlations")
  
  plot <- ggplot(data, aes(x = namex, y = namey, fill = r)) + 
    geom_tile() + scale_fill_viridis(direction = 1,option="E") + 
    facet_wrap2(~ sex + phase, strip = strip, nrow = nrow) +
    geom_text(aes(label = paste(round((r), 2), "\n", label), colour = r)) +
    scale_colour_viridis(direction = -1,option="E") + 
    theme_bw() +
    guides(colour = "none") +
    {if (comparison == "n") scale_y_discrete(expand = c(0,0), labels = c("length_mm" = "Body length\n (mm)", "n_yawns" = "Number\n of yawns", "hours" = "Total sleep\n (hours)"))} +
    {if (comparison == "n") scale_x_discrete(expand = c(0,0), labels = c("length_mm" = "Body length\n (mm)", "n_yawns" = "Number\n of yawns", "hours" = "Total sleep\n (hours)"))} +
    {if (comparison == "duration") scale_y_discrete(expand = c(0,0), labels = c("duration" = "Duration\n (sec)", "next_bout" = "Time until\n next bout", "mean_y_nt" = "y-position"))} +
    {if (comparison == "duration") scale_x_discrete(expand = c(0,0), labels = c("duration" = "Duration\n (sec)", "next_bout" = "Time until\n next bout", "mean_y_nt" = "y-position"))} +
    labs(title = title, fill = "R", x = "Variable 1", y = "Variable 2") +
    theme(panel.grid = element_blank(), text = element_text(size = 16), axis.text.y = element_text(angle = 90, hjust = 0.5))
  
  return(plot)
}

