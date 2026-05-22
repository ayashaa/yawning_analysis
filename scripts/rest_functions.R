fillNA <- function(als_data = als_data, columns) {
  full_time <- data.frame(datetime = seq(from = als_data$datetime[1], to = als_data$datetime[nrow(als_data)], by = "1 sec"))
  df_filled <- left_join(full_time, als_data, by = "datetime")
  #df_filled$sample_id <- df_filled$sample_id[1]
  
  df_filled[columns] <- lapply(df_filled[columns], function(col) {
    if (is.numeric(col)) { return(na.approx(col, na.rm = FALSE)) } 
    else { return(col) }
  })

  return(df_filled)
}

# Inserts rest column based on activity in rolling 60 sec windows (default = < 5% movement above 15 mm/s)
getRest <- function(als_data = als_data, threshold = 15, pct = 0.05, n = 600) {
  #require("runner")
  tic("rest calculated")
  
  n <- n 
  
  if (n == 600) { als_data <- als_data[!is.na(als_data$speed_mm),] }
  if (n == 60) { als_data <- als_data[!is.na(als_data$mean_speed_mm),] }
  als_data <- tail(als_data, -1)
  
  if (n == 600) { als_data$movement <- als_data$speed_mm > threshold }
  if (n == 60) { als_data$movement <- als_data$mean_speed_mm > threshold }
  
  als_data$rest <- frollmean(als_data$movement, n = n, align = "center") < pct
  
  als_data <- tail(als_data, -(n/2))
  als_data <- head(als_data, -(n/2))
  toc()
  
  return(als_data)
}

# Inserts phase (dawn/day/dusk/night) for each observation, day/night times can be adjusted but defaulted to 7am and 7pm
getPhase <- function(als_data = als_data, day = "7:00:00", night = "19:00:00") {
  tic("phase determined")
  
  day <- as.POSIXct(paste('1970-01-01 ', day), format = "%Y-%m-%d %H:%M:%S", origin = "1970-01-01 00:00:00", tz = "GMT")
  night <- as.POSIXct(paste('1970-01-01 ', night), format = "%Y-%m-%d %H:%M:%S", origin = "1970-01-01 00:00:00", tz = "GMT")
  
  dawn_start <- format(day - hours(1), format="%H:%M:%S")
  dawn_end <- format(day + hours(1), format="%H:%M:%S")
  
  dusk_start <- format(night - hours(1), format="%H:%M:%S")
  dusk_end <- format(night + hours(1), format="%H:%M:%S")
  
  als_data <- mutate(als_data, phase = case_when(between(format(as.POSIXct(datetime, origin = "1970-01-01 00:00:00"), format="%H:%M:%S", origin = "1970-01-01 00:00:00"), dusk_start, dusk_end) ~ 'dusk',
                                                 between(format(as.POSIXct(datetime, origin = "1970-01-01 00:00:00"), format="%H:%M:%S", origin = "1970-01-01 00:00:00"), dawn_start, dawn_end) ~ 'dawn',
                                                 between(format(as.POSIXct(datetime, origin = "1970-01-01 00:00:00"), format="%H:%M:%S", origin = "1970-01-01 00:00:00"), dawn_end, dusk_start) ~ 'day',
                                                 .default = 'night'))
  
  toc()
  
  return(als_data)
}

findBoutPhase <- function(start_time, start_phase, end_time, end_phase, day = "07:00:00", night = "19:00:00") {
  day <- as.POSIXct(paste(format(start_time, "%Y-%m-%d"), day), format = "%Y-%m-%d %H:%M:%S", origin = "1970-01-01 00:00:00", tz = "GMT")
  night <- as.POSIXct(paste(format(start_time, "%Y-%m-%d"), night), format = "%Y-%m-%d %H:%M:%S", origin = "1970-01-01 00:00:00", tz = "GMT")
  
  start_time <- as.POSIXct(start_time, origin = "1970-01-01 00:00:00", tz = "GMT")
  end_time <- as.POSIXct(end_time,  origin = "1970-01-01 00:00:00", tz = "GMT")
  
  dawn_start <- day - hours(1)
  dawn_end <- day + hours(1)
  
  dusk_start <- night - hours(1)
  dusk_end <- night + hours(1)
  
  #if (difftime(start_time, end_time, units = "hours") > 24) { return(NA) }
  
  if (start_phase == end_phase & day(start_time) == day(end_time)) { return(start_phase) }
  else if (start_phase == "dawn" & (end_phase == "dusk" | end_phase == "night")) { return("day") }
  else if (start_phase == "dusk" & (end_phase == "dawn" | end_phase == "day")) { return("night") }
  else if (start_phase == "night" & end_phase == "dusk") { return("day") }
  else if (start_phase == "day" & end_phase == "dawn") { return("night") }
  else {
    next_phase_time <- case_when(end_phase == "night" ~ dusk_end, end_phase == "dawn" ~ dawn_start,
              			 end_phase == "day" ~ dawn_end, end_phase == "dusk" ~ dusk_start)
    
    if (day(start_time) < day(end_time)) { next_phase_time <- next_phase_time + day(1) }
    
    if (start_phase == "night" & end_phase == "day") { next_phase_time <- day }
    else if (start_phase == "day" & end_phase == "night") { next_phase_time <- night }
    
    phase <- ifelse(abs(as.numeric(difftime(start_time, next_phase_time, units = "secs"))) > 
                    abs(as.numeric(difftime(end_time, next_phase_time, units = "secs"))), 
                    start_phase, end_phase)
  }
  
  return(phase)
}

findDay <- function(start_time, end_time) {
  day_end = as.POSIXct(paste(format(end_time, "%Y-%m-%d"), "00:00:00"), format = "%Y-%m-%d %H:%M:%S", origin = "", tz = "GMT")
  
  day <- ifelse(abs(as.numeric(difftime(start_time, day_end, units = "secs"))) > 
                abs(as.numeric(difftime(end_time, day_end, units = "secs"))), 
                day(start_time), day(end_time))
}

# Creates new dataframe with the bout structure (length, start and end times, and phases)
boutStructure <- function(als_data = als_data, min_bout = 1) {
  tic("bout structure determined")
  als_data$rest <- ifelse(als_data$rest == TRUE, "rest", "active")
  
  seq <- rle(als_data$rest)
  bouts <- as.data.frame(cbind(seq$values, seq$lengths))
  colnames(bouts) <-  c("state", "length")
  bouts$length <- as.numeric(bouts$length)
  #bouts <- as.data.frame(bouts)
  
  bouts$day <- NA
  bouts$start <- NA
  bouts$start_phase <- NA
  bouts$end <- NA
  bouts$end_phase <- NA
  bouts$sample_id <- NA
  bouts$species_six <- NA
  bouts$tribe <- NA
  bouts$diel <- NA
  bouts$mean_y_pos <- NA

  time_elapsed = 0
  
  # Assign start and end times and phases
  for (i in 1:nrow(bouts)) {
    time_elapsed <- time_elapsed + bouts$length[i]
    
    bouts$end[i] <- als_data$datetime[time_elapsed]
    bouts$end_phase[i] <- als_data$phase[time_elapsed]
    
    if (i != nrow(bouts)) {
      bouts$start[i + 1] <- als_data$datetime[time_elapsed + 1]
      bouts$start_phase[i + 1] <- als_data$phase[time_elapsed + 1]
    }
    
    else { 
      last_bout_index <- nrow(als_data) - bouts$length[i] + 1
      bouts$start[i] <- als_data$datetime[last_bout_index]
      bouts$start_phase[i] <- als_data$phase[last_bout_index]
    }

    bouts$mean_y_pos[i] <- mean(als_data$mean_y_nt[time_elapsed + 1 - bouts$length[i]]:als_data$mean_y_nt[time_elapsed])
  }
  bouts$start[1] <- als_data$datetime[1]
  bouts$start_phase[1] <- als_data$phase[1]
  
  bouts$start <- as.POSIXct(bouts$start, origin = "1970-01-01 00:00:00", format = "%Y-%m-%d %H:%M:%S")
  bouts$end <- as.POSIXct(bouts$end, origin = "1970-01-01 00:00:00", format = "%Y-%m-%d %H:%M:%S")
  
  bouts <- mutate(bouts, length = as.numeric(difftime(bouts$end, bouts$start, unit = "secs")))
  #bouts <- mutate(bouts, day = ifelse(day(bouts$start) == day(bouts$end), day(bouts$start), 
  #                                    findDay(bouts$start, bouts$end)))
  
 
  
  bouts$sample_id <- unique(als_data$sample_id)
  bouts$species_six <- unique(als_data$species_six)
  bouts$tribe <- unique(als_data$tribe)
  bouts$diel <- als_data$diel[1]
  
  #bouts$overall_phase <- mapply(findBoutPhase, bouts$start, bouts$start_phase, bouts$end, bouts$end_phase)
  
  #total <- sum(bouts$length)
  #bouts$proportion <- bouts$length/total
  bouts <- filterBouts(bouts, min_bout) 

  toc()
  
  return(bouts)
}


filterBouts <- function(bout_data = bout_data, threshold = 0) {
  # Remove all bouts below threshold
  bouts <- filter(bout_data, bout_data$length > threshold)
  
  state_seq <- bouts$state
  state_length <- (rle.state_seq = rle(state_seq)$lengths)
  
  # Adjust bout lengths after removing short bouts
  for (i in 1:length(state_length)) {
    if (state_length[i] > 1) { 
      bouts$end[i] <- bouts$end[i - 1 + state_length[i]]
      #bouts$mean_y_pos[i] <- mean(bouts$mean_y_pos[i]:bouts$mean_y_pos[i - 1 + state_length[i]])
      bouts$length[i] <- difftime(bouts$end[i], bouts$start[i], unit = "secs")
      ifelse(state_length[i] == 2, bouts <- bouts[-(i + 1),], bouts <- bouts[-((i + 1):(i - 1 + state_length[i])),]) 
      #bouts$mean_y_pos <- mean(bouts$mean_y_pos[state_length[i]]:bouts$mean_y_pos[i - 1 + state_length[i]])
    }
  }
  
  return(bouts)
}

# Returns table with daily summaries of total rest, bout frequency and phase ratio 
boutSummary <- function(bout_data = bout_data, summarise_by = c("day", "week"), include_crepuscular = FALSE) {
  tic("bouts summarised")
  phases <- unique(bout_data$overall_phase)
  
  if (include_crepuscular == FALSE) {
    bout_data <- filter(bout_data, overall_phase == "day" | overall_phase =="night")
  }
  
  daily <- group_by(bout_data, day, tribe, species_six, sample_id, state, overall_phase) 
  
  summary <- summarise(daily, total_sec = sum(length), total_hour = total_sec/3600, freq = length(state), 
                       mean_bout_length = mean(length), longest_bout = max(length), sfi = freq/total_sec,
                       pct_cons = mean(length)/sum(length)*100,
                       max_cons = max(length)/total_sec*100, mean_y = mean(mean_y_pos))
  
  if (summarise_by == "week") { 
    weekly <- summary

    for (days in 1:8) {
    	if (sum(weekly[which(weekly$day == days), 'total_hour']) < 18) { weekly <- filter(weekly, day != days) }
    }
    #if (sum(summary[which(summary$day == 8), 'total_hour']) < 18) { weekly <- filter(weekly, day != 8) }
    #if (sum(summary[which(summary$day == 1), 'total_hour']) < 18) { weekly <- filter(weekly, day != 1) }
    
    weekly <- group_by(weekly, tribe, species_six, sample_id, state, overall_phase) 
    
    week_summary <- summarise(weekly, mean_total_sec = mean(total_sec), absolute_total_sec = sum(total_sec), mean_total = mean(total_hour), mean_freq = mean(freq), 
                              mean_bout_length = mean(mean_bout_length), mean_longest_bout = mean(longest_bout), abs_longest = max(longest_bout), mean_sfi = mean(sfi),
                              pct_cons = mean_bout_length/mean(mean_total_sec)*100, max_cons = mean_longest_bout/mean_total_sec*100, abs_max_cons = abs_longest/mean_total_sec*100,
                              mean_max_cons = mean(max_cons), mean_y = mean(mean_y))
  }
  
  toc()
  if (summarise_by == "day") { return(summary) }
  if (summarise_by == "week") { return(week_summary) }
}

getmode <- function(v) {
  uniqv <- unique(v)
  
  return(uniqv[which.max(tabulate(match(v, uniqv)))])
}

N50consolidation <- function(bout_lengths = bout_lengths) {
  
  bout_lengths <- sort(bout_lengths, TRUE)
  
  bout_sum <- cumsum(bout_lengths)
  bout_ind <- which(bout_sum > sum(bout_lengths)*0.5)
  
  L <- min(bout_ind)
  
  return(L)
}

L50consolidation <- function(bout_lengths = bout_lengths) {
  bout_lengths <- sort(bout_lengths, TRUE)
  return(bout_lengths[N50consolidation(bout_lengths)])
}


# Write CSVs for bout data
restData <- function(als_data = als_data, speed_threshold = 15, pct_movement = 0.05, day = "07:00:00", night = "19:00:00", min_bout = 4) {
  require("tictoc")
  
  als_data <- getRest(als_data, threshold = speed_threshold, pct = pct_movement)
  als_data <- getPhase(als_data, day, night)
  
  rest_data <- vector(mode='list', length=3)
  
  rest_data$bout_data <- boutStructure(als_data, min_bout)
  rest_data$daily_summary <- boutSummary(rest_data$bout_data, "day")
  rest_data$weekly_summary <- boutSummary(rest_data$bout_data, "week")
  
  write.csv(rest_data$bout_data, paste("/home/ayasha/projects/def-mshafer/ayasha/cichlid_bout_data/filtered5/y_pos/", als_data$sample_id[1], "_bout_structure.csv", sep = ""))
  write.csv(rest_data$daily_summary, paste("/home/ayasha/projects/def-mshafer/ayasha/cichlid_bout_data/filtered5/y_pos/", als_data$sample_id[1], "_daily_bout_summary.csv", sep = ""))
  write.csv(rest_data$weekly_summary, paste("/home/ayasha/projects/def-mshafer/ayasha/cichlid_bout_data/filtered5/y_pos/", als_data$sample_id[1], "_weekly_bout_summary.csv", sep = ""))
  
  return(rest_data)
}
