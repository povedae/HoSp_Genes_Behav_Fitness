
# Code from CK Hui. THANK YOU!
library(dplyr)
library(lubridate)
library(ggplot2)
library(data.table)

## DETERMINE ONSET AND OFFSET TIMES ##

#load data
alltags <- read.csv("~/HoSp Data/prepared datasets/alltags.csv")

# Combine the date and time columns into a single datetime column
alltags <- alltags %>% mutate(datetime = as.POSIXct(paste(date, time), format="%Y-%m-%d %H:%M:%S"))

# Arrange and group data, then summarize by the minute
df_clean <- alltags %>%  arrange(ID) %>%                                  
  mutate(datetime = floor_date(datetime, "minute"),  # Round `datetime` to minute precision
    datetime_str = format(datetime, "%Y-%m-%d %H:%M")) %>%  # Format to remove seconds
  group_by(date, ID, datetime) %>%  # Group by formatted datetime
  summarise(activity = sum(activity), .groups = "drop")  # Summarize activity over full minute

df_clean$datetime <- as.POSIXct(df_clean$datetime_str, format = "%Y-%m-%d %H:%M")
df_clean$datetime_str <- NULL


#Split birds into their own data frame
unique_ids <- unique(df_clean$ID)

for (id in unique_ids) {
  df_name <- paste("ID", id, sep = "_")
  assign(df_name, subset(df_clean, ID == id))
}

######
# Now all the birds are in different data frames
######

############################################
# How to determine activity on- and offset #
############################################
crossings<-function(ID_num){
  ID_num <- ID_num %>% arrange(datetime)
  
  #rolling means for 20 minutes
  roll_ID<- ID_num %>% 
    mutate(rolling_avg = data.table::frollmean(activity,21),index=seq(1,nrow(ID_num),1))
  roll_ID$rolling_avg[is.na(roll_ID$rolling_avg)] <- 0
  
  # daily means
  daily_mean_ID <- ID_num %>%
    group_by(date) %>%
    summarise(mean_activity = mean(activity, na.rm = TRUE))
  
  # Set threshold for on- and offset marker
  # mean across all days used because means are expected to be similar
  threshold <- mean(daily_mean_ID$mean_activity)
  
  # Placeholder for crossing data
  crossings_data <- data.frame(datetime = as.POSIXct(character()), value = numeric())
  
  # Loop to calculate interpolated crossing points
  for (i in seq_along(roll_ID$rolling_avg)[-1]) {
    if ((roll_ID$rolling_avg[i-1] < threshold && roll_ID$rolling_avg[i] > threshold) ||
        (roll_ID$rolling_avg[i-1] > threshold && roll_ID$rolling_avg[i] < threshold)) {
      # Linear interpolation for more precise timestamp
      time_diff <- diff(as.numeric(roll_ID$datetime[(i-1):i]))
      value_diff <- diff(roll_ID$rolling_avg[(i-1):i])
      proportion <- abs(threshold - roll_ID$rolling_avg[i-1]) / value_diff
      exact_time <- as.POSIXct(roll_ID$datetime[i-1] + time_diff * proportion, origin = "1970-01-01")
      
      # Store the crossing point
      crossings_data <- rbind(crossings_data, data.frame(datetime = exact_time, value = threshold))
    }
  }
  
  # Ensure datetime format for plotting
  crossings_data$datetime <- as.POSIXct(crossings_data$datetime, origin = "1970-01-01")
  
  #Plot to see if point are accurate and which ones to ignore
  plot_ID<-ggplot(roll_ID, aes(x = datetime, y = rolling_avg)) +
    geom_line(color = "black") +  # Plot the rolling average
    geom_hline(yintercept = threshold, color = "red", linetype = "dashed") +  # Plot the threshold line
    geom_point(data = crossings_data, aes(x = datetime, y = value), color = "blue", size = 3, shape = 19) +  # Plot crossing points
    labs(title = "Rolling Average and Interpolated Crossings",
         x = "Time",
         y = "Rolling Average") +
    theme_minimal()
  
  #Here is your data (select only one around day activity)
  print(plot_ID)
  print(crossings_data)
  
}

##########################

## Run through each ID object with:
crossings(ID_27)

  # use your best judgement to choose which intercept is most feasible when there are multiple
  # record onset and offset times in excel sheet also with Julian day, site, sex, ID, and smi.
    # find difference between onset/offset and sunrise/sunset for relative onset/offset.

