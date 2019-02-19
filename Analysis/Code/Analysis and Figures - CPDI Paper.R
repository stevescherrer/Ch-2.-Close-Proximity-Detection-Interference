  ###### Script for CPDI Doughnut Effect Paper - Producing analysis results and figures for publication
  
    #### Written by: Stephen R. Scherrer 
    #### 26 January 2017
  
  ####### Cleaning workspace and setting directories ----
    ### Clearing Workspace
    rm(list=ls()) # Clear workspace
  
    ### Starting Script Timer
    script_timer <- proc.time()
    
    ### Linking to Project Directories
      setwd('/Users/stephenscherrer/Google Drive/Weng Lab/Manuscripts/Scherrer - CPDI Model Paper/CPDI Submission Repository')
      project_dir = getwd()
      data_dir = file.path(project_dir, 'Data')
      results_dir = file.path(project_dir, 'Results')
      figure_dir = file.path(project_dir, 'Figures')
      source_dir = file.path(project_dir, 'Code')
    
    ### Setting Project Directory 
    setwd(project_dir)
    
    ### Establishing History 
    savehistory(file= file.path(results_dir, "Rhistory"))
  
  ####### Importing principle dependencies ----
    # install.packages('geosphere')
    library('geosphere') # distGeo() Note: wrapped in old lldist function
    # install.packages('reshape')
    library('reshape') # melt()
    # install.packages('MuMIn')
    library('MuMIn') # AICc()
    # install.packages('dplyr')
    library('dplyr') # filter()
    #install.packages('doParallel')
    library('doParallel') # do_parallel()
    #install.package('lubridate')
    library('lubridate') # floor_date()
    ## install.packages('beepr')
    library('beepr') # beep()
    ## install.packages('notifyR')
    library('notifyR') # send_push()
    ## install.packages('ggplot2')
    library('ggplot2') # geom_bar()
    ## install.packages('data.table') 
    library('data.table') # uniqueN()
    ## install.packages('marmap')
    library('marmap') # readBathy(), subsetBathy(), getNOAAbathy(), plot.bathy(), scaleBathy()
    ## install.packages('data.table')
    library('data.table') # uniqueN()  
    # install.packages('mgcv')
    library('mgcv') # gam()
  
  ### Sourcing Mechanistic CPDI Model
    source(file.path(source_dir, 'Mechanistic Model Implemented - R/Mechanistic CPDI Model.R')) # model_receiver_interferenece()
  
  ####### Setting Up Parallel Environment ----
    ### Setting up parallel processing for later analysis using pdredge() and foreach()
    ## Setting number of cores
    n_cores = detectCores() # Default to the number of cores available on the computer
    ## Creating a cluster
    clust = makeCluster(n_cores)
    clusterEvalQ(clust, library('MuMIn')) # loading MuMIn package to clusters for using pdredge function
    clusterEvalQ(clust, library('mgcv')) # loading mgcv package to clusters for using GAM function
    
  ####### Utility Functions ----
  datenum2posix = function(x, timez = "HST") {
    ## Function to convert matlab datenum to R posix date
    # Author: Luke Miller   Feb 20, 2011
      # Convert a numeric  MATLAB datenum (days since 0000-1-1 00:00) to seconds in 
      # the Unix epoch (seconds since 1970-1-1 00:00). Specify a time zone if the 
      # input datenum is anything other than the GMT/UTC time zone. 
    days = x - 719529 	# 719529 = days from 1-1-0000 to 1-1-1970
    secs = days * 86400 # 86400 seconds in a day
      # This next string of functions is a complete disaster, but it works.
      # It tries to outsmart R by converting the secs value to a POSIXct value
      # in the UTC time zone, then converts that to a time/date string that 
      # should lose the time zone, and then it performs a second as.POSIXct()
      # conversion on the time/date string to get a POSIXct value in the user's 
      # specified timezone. Time zones are a goddamned nightmare.
    return(as.POSIXct(strftime(as.POSIXct(x = secs, units = 'secs', origin = '1970-1-1', 
                                          tz = 'UTC'), format = '%Y-%m-%d %H:%M:%S', 
                               tz = 'UTC', usetz = FALSE), tz = timez))
  }
  
  convert_tz = function(datetime, new.tz = 'HST'){
    ## Function to convert GMT/UTC times to HST time
    datetime.new.tz = strptime(datetime, format = '%Y-%m-%d %H:%M:%S', tz = new.tz)
    dateoffset = datetime-datetime.new.tz
    datetime.new.tz = datetime.new.tz + dateoffset
    return(datetime.new.tz)
  }
  
  convert_lat_lon = function(ll_deg, ll_min = FALSE){
    ## Converts latitude and longitude between ll minutes and ll decimal degrees
    # 2 usages:
      # Convert decimal degrees to degree minutes
        # 1 argument
          # ll_pref is a single argument of latitude or longitude in decimal degrees
          # Returns a prefix and decimal for that argument
      # Convert degree minutes to decimal degrees
        # 2 arguments
          # ll_pref is the latitude or longitude's degree
          # ll_min is the degree minutes
          # returns a single float of ll in decimal degrees
    if (ll_min[1] == FALSE){ #then we are going from one number to two
      ll_deg = as.numeric(as.character(ll_deg))
      ll_bin = matrix(0, length(ll_deg), 2)
      for (r in 1:length(ll_deg)){
        if (isTRUE(ll_deg[r] >= 0)){
          ll_dec = ll_deg[r] - floor(ll_deg[r])
          ll_bin[r, ] = c(floor(ll_deg[r]), (ll_dec)*60)
        } else {
          ll_dec = (ll_deg[r] - ceiling(ll_deg[r]))*-1
          ll_bin[r, ] = c(ceiling(ll_deg[r]), (ll_dec)*60)
        }
      }
    }else{ #if we are converting from two numbers to one
      ll_deg = as.numeric(as.character(ll_deg))
      ll_min = as.numeric(as.character(ll_min))
      ll_bin = matrix(0, length(ll_deg), 1)
      for (r in 1:length(ll_deg)){
        ll_dec_deg = abs(ll_deg[r]) + (abs(ll_min[r])/60)
        if (isTRUE(ll_deg[r] < 0)){
          ll_dec_deg = ll_dec_deg*(-1)
        }
        ll_bin[r] = ll_dec_deg
      }
    }
    return (ll_bin)
  }
  
  lldist = function(point1, point2){
    ## Wrapper function to calculate distance between two lat lon points on a geodesic sphere.
      # Requires distGeo() from the geosphere package
    distance = distGeo(p1 = point1, p2 = point2) / 1000
    return(distance)
  }
  
  clean_receiver = function(receiver){
    ## Returns receiver serial number as a factor, remvoing the 'VR2W-' Prefix
    cleaned_receiver = as.factor(substring(receiver, 6))
    return (cleaned_receiver)
  }
  
  load_receiver = function(filename, filepath = FALSE){
    ## Function to load in receiver data .csv file
    proj_dir = getwd()
    if (filepath != FALSE){
      setwd(filepath)
    }
    receiver_dates = receiver_col_names(read.csv(filename))
    receiver_dates$deployment_date = strptime(receiver_dates$deployment_date, format = '%m/%d/%y %H:%M', tz = 'HST')
    receiver_dates$recovery_date = strptime(receiver_dates$recovery_date, format = '%m/%d/%y %H:%M', tz = 'HST')
    receiver_dates$lat = convert_lat_lon(receiver_dates$lat_deg, receiver_dates$lat_min)
    receiver_dates$lon = convert_lat_lon(receiver_dates$lon_deg, receiver_dates$lon_min)
    setwd(proj_dir)
    return (receiver_dates)
  }
  
  receiver_col_names = function(receiver_file_raw){
    ## Function that renames columns for an imported receiver data file
    receiver_file = receiver_file_raw
    colnames(receiver_file)[1] = 'serviced'
    colnames(receiver_file)[2] = 'station_name'
    colnames(receiver_file)[3] = 'consecutive_deployment_number'
    colnames(receiver_file)[4] = 'deployment_date'
    colnames(receiver_file)[5] = 'recovery_date'
    colnames(receiver_file)[6] = 'recovered'
    colnames(receiver_file)[7] = 'in_data_set'
    colnames(receiver_file)[8] = 'lat_deg'
    colnames(receiver_file)[9] = 'lat_min'
    colnames(receiver_file)[10] = 'lon_deg'
    colnames(receiver_file)[11] = 'lon_min'
    colnames(receiver_file)[12] = 'depth'
    colnames(receiver_file)[13] = 'vr2w_serial'
    colnames(receiver_file)[14] = 'acoustic_release_serial'
    colnames(receiver_file)[15] = 'acoustic_release_battery_life'
    colnames(receiver_file)[16] = 'acoustic_release_voltage_at_deployment'
    colnames(receiver_file)[17] = 'acoustic_release_serial_code'
    colnames(receiver_file)[18] = 'temperature_logger_serial'
    colnames(receiver_file)[19] = 'location_code'
    colnames(receiver_file)[20] = 'deployed_by'
    colnames(receiver_file)[21] = 'recovered_by'
    colnames(receiver_file)[22] = 'comments_deployment'
    colnames(receiver_file)[23] = 'comments_recovery'
    return (receiver_file)
  }
  
  clean_vue_lat_lon = function(vue_data_df, receiver_data_df){
    ## Function that replaces the lat lon positions a receiver may have been initialized with, to field records of coordinates for the actual deployment location.
    station = rep(NA, times = length(vue_data_df$datetime))
    for (i in 1:length(receiver_data_df$station_name)){
      receiver_subset_index = which(vue_data_df$receiver == receiver_data_df$vr2w_serial[i])
      deploy_subset_index = which(vue_data_df$datetime >= receiver_data_df$deployment_date[i])
      recover_subset_index = which(vue_data_df$datetime < na.omit(receiver_data_df$recovery_date[i]))
      ind = Reduce(intersect, list(receiver_subset_index, deploy_subset_index, recover_subset_index))
      vue_data_df$lat[ind] = receiver_data_df$lat[i]
      vue_data_df$lon[ind] = receiver_data_df$lon[i]
      station[ind] = as.character(receiver_data_df$station_name[i])
    }
    vue_data_df$station = as.factor(station)
    return(vue_data_df)
  }
  
  clean_tag_id = function(tag_id){
    ## Function that returns tag ID number as a factor, removing the 'A69-####-' prefix
    cleaned_id = as.factor(substring(tag_id, 10))
    return (cleaned_id)
  }
  
  vue_col_names = function(vue_data_raw){
    # Function that renames columns for an imported vemco data set
    colnames(vue_data_raw)[1]  <- 'datetime'
    colnames(vue_data_raw)[2]  <- 'receiver'
    colnames(vue_data_raw)[3]  <- 'tag_id'
    colnames(vue_data_raw)[4]  <- 'name'
    colnames(vue_data_raw)[5]  <- 'tag_serial'
    colnames(vue_data_raw)[6]  <- 'sensor_value'
    colnames(vue_data_raw)[7]  <- 'sensor.unit'
    colnames(vue_data_raw)[8]  <- 'station'
    colnames(vue_data_raw)[9]  <- 'lat'
    colnames(vue_data_raw)[10] <- 'lon'
    return (vue_data_raw)
  }
  
  load_vemco = function(filename, filepath = FALSE, format = '%Y-%m-%d %H:%M:%S'){
    ## Function that loads in Vemco Database from a VUE Export file 
    proj_dir = getwd()
    if (isTRUE(filepath != FALSE)) {setwd(filepath)}
    vue_data_raw = read.csv(filename)
    vue_data_cleaned = vue_col_names(vue_data_raw)
    vue_data_cleaned$datetime = strptime(vue_data_cleaned$datetime, 
                                         format = format,
                                         tz = "GMT")
    vue_data_cleaned$datetime = convert_tz(vue_data_cleaned$datetime, new.tz = 'HST')
    vue_data_cleaned$tag_id = clean_tag_id(vue_data_cleaned$tag_id)
    vue_data_cleaned$receiver = clean_receiver(vue_data_cleaned$receiver)
    setwd(proj_dir)
    return (vue_data_cleaned)
  }
  
  compute_standard_error = function(numeric_vector){
    stderror = sd(numeric_vector) / sqrt(length(numeric_vector))
    return(stderror)
  }
  
  as.number = function(x){
    if(class(x) == 'factor'){
      x = levels(x)[x]
    }
    as.numeric(x)
  }
  
  ####### General Parameters ----
  expected_detections = data.frame(
    ## From: https://vemco.com/collision-calculator/
    'n_tags' = c(1:18),
    'hourly_detections_total' = c(57, 102, 138, 165, 186, 200, 210, 216, 218, 218, 215, 211, 205, 198, 191, 183, 175, 166),
    stringsAsFactors = FALSE
  )
  expected_detections$hourly_detections_per_tag = expected_detections$hourly_detections_total / expected_detections$n_tags
  
  ####### Experiment 1: Deep Water Range Test ----
  
    #### Loading In and Cleaning Data Files ----
    setwd(data_dir)
  
    ### Loading Receiver Data
    receiver_data = load_receiver(filename = 'DEPLOYMENT_RECOVERY_LOG.csv')
    # dim(receiver_data)
      # 217  28
    
    ### Loading VUE data 
    vue_data = load_vemco(filename = 'Range_Test_June_2014_All_Receivers.csv', format = '%m/%d/%y %H:%M')
    #  dim(vue_data)
      # 130284     10
    
    ### Importing Windspeed Data from from NOAA Honolulu station Station ID 1612340
      # http://www.ndbc.noaa.gov/view_text_file.php?filename=oouh1h2014.txt.gz&dir=data/historical/stdmet/
      # parent: http://www.ndbc.noaa.gov/station_history.php?station=oouh1
      # Data dates range between dates 2014-01-01 00:00:00 and 2014-12-31 23:54:00
        # with the following format: 
        # YY  MM DD hh mm WDIR WSPD GST  WVHT   DPD   APD MWD   PRES  ATMP  WTMP  DEWP  VIS  TIDE
      # with the following units:   
        # yr  mo dy hr mn degT m/s  m/s     m   sec   sec degT   hPa  degC  degC  degC   mi    ft
      wind_data = as.data.frame(read.table("oouh1h2014.txt"), header = TRUE, sep = " ")
      # dim(wind_data)
        # 85864    18
      colnames(wind_data) = c("YY",  "MM", "DD", "hh", "mm", "WDIR", "WSPD", "GST",  "WVHT",   "DPD",   "APD", "MWD",   "PRES",  "ATMP",  "WTMP",  "DEWP",  "VIS",  "TIDE")
      ### Convert dates and times to POSIXct format. (Station data time is in local (UTC) time by default)
        wind_data$datetime = as.POSIXct(paste(wind_data$YY, "-", wind_data$MM, "-", wind_data$DD, " ", wind_data$hh, ":", wind_data$mm, sep = ""), format = "%Y-%m-%d %H:%M", tz = "GMT")
      ### Convert datetime to HST
        wind_data$datetime = strftime(wind_data$datetime, format = "%Y-%m-%d %H:%M")
        wind_data$hourly = strftime(wind_data$datetime, format = "%Y-%m-%d %H")
        
    ### Grouping windspeed data into hourly means (mean_wspd)
      wind_data_grouped = group_by(wind_data, hourly)
      wspd_by_hour = summarize(wind_data_grouped, mean_wspd = mean(WSPD))
      wspd_by_hour$hourly = as.POSIXct(wspd_by_hour$hourly, format = "%Y-%m-%d %H")
      
    ### Grouping windgust data into hourly means
      gst_by_hour = summarize(wind_data_grouped, mean_gst = mean(GST))
      gst_by_hour$hourly = as.POSIXct(gst_by_hour$hourly, format = "%Y-%m-%d %H")
    
    ### Importing tide data from: NOAA Honolulu station, accessed 28 OCT 2015
      # http://tidesandcurrents.noaa.gov/waterlevels.html?id=1612340&units=standard&bdate=20110930&edate=20100930&timezone=LST&datum=MLLW&interval=h&action=## Meta data available here: 
      # Data accessed for dates ranging between 2013-09-30 00:00:00 HST and 2014-09-30 13:00:00 HST
      # Data is already in HST
        tide_data = read.csv('noaa tide data 2013-2014.csv')
        # dim(tide_data)
          # 8774    5
        colnames(tide_data) = c("hour_bin", "water_level", 'sigma', 'I', "L")
        tide_data$hour_bin = as.POSIXct(tide_data$hour_bin, format = "%Y-%m-%d %H:%M")
      ### Determining if tide was going in or out.
      ## First assume all data is headed in
        tide_data$direction = "In"
      ## Then loop through data. If current water level is less than previous water level, tide is going out.
        for(i in 2:length(tide_data$direction)){
          if(tide_data$water_level[i-1] > tide_data$water_level[i]){
          tide_data$direction[i] = "Out"
          }
        }
        tide_data$direction = as.factor(tide_data$direction)
        
      ### Importing Tag Deployment Meta Data File
      tag_meta_data = read.csv('Range Test June 2014 - Tag Meta Data Locations and Depth.csv')
      # dim(tag_meta_data)
        # 12  9
      ### Converting tag_meta_data Lat Lons from degree minutes to decimal degrees
        tag_meta_data$lat = convert_lat_lon(tag_meta_data$lat_deg, tag_meta_data$lat_min)
        tag_meta_data$lon = convert_lat_lon(tag_meta_data$lon_deg, tag_meta_data$lon_min)
        
      ### Determining distance of tag from receiver string using lon and lat waypoints from gps
        for (i in 1:length(tag_meta_data$distance)){
          tag_meta_data$distance[i] = round(1000*(lldist(point1 = c(tag_meta_data$lon[i], 
                                                              tag_meta_data$lat[i]), 
                                                   point2 = c(receiver_data$lon[receiver_data$station_name == 'Range Test - June 2014 - Diamond Head 1m'][1], 
                                                              receiver_data$lat[receiver_data$station_name == 'Range Test - June 2014 - Diamond Head 1m'][1]))))
        }
        unique(tag_meta_data$distance)
        # 0 199 399 578 766 959
    
      ### Removing any other receivers in the database
        receiver_1m = 123736 # serial number of receiver in 1 m depth
        receiver_30m = 123732 # serial number of receiver in 30 m depth
        vue_data = vue_data[which(vue_data$receiver %in% c(receiver_1m, receiver_30m)), ]
      
      ### Removing detections prior to experiment (occurred durring experiment set up on boat) and determining an end to the experiment.
        start_date = as.POSIXct('2014-06-07 8:00:00', tz = "HST")
        # "2014-06-07 08:00:00 HST"
        end_date = as.POSIXct("2014-06-26 05:00:00", format = "%Y-%m-%d %H:%M:%S")
        # "2014-06-26 05:00:00 HST"
        vue_data = vue_data[which(vue_data$datetime >= start_date & vue_data$datetime < end_date), ]
    
      ### Rounding off of final partial hour to account for transmissions received while gear was recovered
        vue_data = vue_data[which(vue_data$datetime >= ceiling_date(min(vue_data$datetime), unit = 'hour') & vue_data$datetime < floor_date(max(vue_data$datetime), unit = 'hour')), ]
        range(vue_data$datetime)
          # "2014-06-07 09:00:00 HST" "2014-06-26 03:59:00 HST"
    
        ## How many total detections occurred durring the course of the experiment? (includes tags not part of range test)
          detections_all = dim(vue_data)[1]
            # 129622
    
      ### Removing tags not in experiment. This could be a tagged fish that swam by or a false detection
        tags_0m    = c(18236, 18237)
        tags_200m  = c(18238, 18239)
        tags_400m  = c(18240, 18241)
        tags_600m  = c(18242, 18243)
        tags_800m  = c(18244, 18245)
        tags_1000m = c(18246, 18247)
        tag_ids = c(tags_0m, tags_200m, tags_400m, tags_600m, tags_800m, tags_1000m) 
        vue_data = vue_data[as.numeric(levels(vue_data$tag_id))[vue_data$tag_id]
                            %in% tag_ids, ] 
      
        ## How many detections were from tags that were part of the experiment?
          detections_experiment = dim(vue_data)[1]
            # 128178
          detections_experiment / detections_all * 100
            # 98.88599 %
      
        ## How many detections were from tags that were not part of the experiment?
          detections_other = detections_all - detections_experiment
            # 1444
          detections_other / detections_all * 100
            # 1.114008
          
      ### Assigning Lon Lat positions to vue_data
          vue_data = clean_vue_lat_lon(vue_data_df = vue_data, 
                                       receiver_data_df = receiver_data)
      ### Getting estimated distances between receivers and tags based on deployment coordinates
      for (i in 1:length(tag_meta_data$lon)){
        tag_meta_data$distance[i] = lldist(point1 = c(unique(vue_data$lon)[1], 
                                                      unique(vue_data$lat)[1]),
                                           point2 = c((tag_meta_data$lon)[i],
                                                      (tag_meta_data$lat)[i])) *
          1000 # m/km
      }
        round(unique(tag_meta_data$distance))
          # 0 199 399 578 766 959
    
    #### Hourly Detection Data Analysis ----
      ### Separating data by tag, receiver, and hour bin
      ## First assign each datum to an date and hour bin
        vue_data$hour_bin = floor_date(vue_data$datetime, unit = "hour")
      ## Then loop through combinations of date/hour, tag id, and receiver to get the number of detections from a tag at a receiver each hour
        detection_df_1 = data.frame()
        registerDoParallel(cores = n_cores)
        detection_df_1 = foreach(b = c(1:length(unique(vue_data$hour_bin))), .combine = rbind) %:% 
          foreach(i = c(1:length(unique(vue_data$tag_id))), .combine = rbind) %:% 
          foreach(r = c(1:length(unique(vue_data$receiver))), .combine = rbind) %dopar%{
            filtered = filter(vue_data, hour_bin == unique(vue_data$hour_bin)[b], 
                              tag_id == unique(vue_data$tag_id)[i],
                              receiver == unique(vue_data$receiver)[r])
            write_line = cbind(as.character(unique(vue_data$hour_bin))[b],
                               as.character(unique(vue_data$tag_id))[i], 
                               as.character(unique(vue_data$receiver))[r], 
                               dim(filtered)[1])
            return(write_line)
          }
        ## Cleaning detection_df_1
        detection_df_1 = as.data.frame(detection_df_1)
        colnames(detection_df_1) = c('hour_bin', 'tag_id', 'receiver', 'detections')
        detection_df_1$detections = as.numeric(as.character(detection_df_1$detections))
          # dim(detection_df_1)
            # 10824     4
        detections_per_hour = aggregate(detection_df_1$tag_id[detection_df_1$detections != 0], by = list(detection_df_1$hour_bin[detection_df_1$detections != 0]), FUN = uniqueN)
        
      ### Combining detection dataframe with distance, bottom depth, height condition, lat, and lon tag_meta_data
        detection_df_1 = merge(detection_df_1, tag_meta_data, 
                            by.x = colnames(detection_df_1) == 'tag_id', 
                            by.y = colnames(tag_meta_data) == 'tag_id')  
        detection_df_1$distance = detection_df_1$distance
        detection_df_1$height_off_bottom = as.factor(detection_df_1$height_off_bottom)
        detection_df_1$recovery_rate = detection_df_1$detections / 60
        detection_df_1$day = as.factor(strftime(detection_df_1$hour_bin, format = "%Y-%m-%d"))
        detection_df_1$time = as.factor(strftime(detection_df_1$hour_bin, format = "%H:%M:%S"))
        colnames(detection_df_1)[colnames(detection_df_1) == 'height_off_bottom'] = "height_of_tag_off_bottom"
        detection_df_1$hour_bin = as.POSIXct(detection_df_1$hour_bin, format = "%Y-%m-%d %H:%M:%S")
  
        ## Calculating the number of unique tags detected during each hour bin
        tags_per_hour = aggregate(detection_df_1$tag_id[detection_df_1$detections != 0], by = list(detection_df_1$hour_bin[detection_df_1$detections != 0]), FUN = uniqueN)
        colnames(tags_per_hour) = c("hour_bin", "tags_per_hour")
        detection_df_1 = merge(detection_df_1, tags_per_hour, by.x = colnames(detection_df_1) == "hour_bin", by.y = colnames(tags_per_hour) == 'hour_bin')
        
      ### Combining detection data with wind speed, and wind gust data
          detection_df_1 = merge(detection_df_1, wspd_by_hour,
                               by.x = colnames(detection_df_1) == 'hour_bin',
                               by.y = colnames(wspd_by_hour) == 'hourly')
          detection_df_1 = merge(detection_df_1, gst_by_hour,
                               by.x = colnames(detection_df_1) == 'hour_bin',
                               by.y = colnames(gst_by_hour) == 'hourly')
          detection_df_1 = merge(detection_df_1, tide_data,
                               by.x = colnames(detection_df_1) == 'hour_bin',
                               by.y = colnames(tide_data) == 'hour_bin')
    
          detection_df_1$hour_bin = as.factor(detection_df_1$hour_bin)
    
      ### Assigning diurnal period to detections based around 12 hour time bins. Before 6am or after 6pm considered night, between 6am and 6pm considered day
        detection_df_1$day_night = 'Night'
        for(i in 1:length(detection_df_1$hour_bin)){
          if(hour(detection_df_1$hour_bin[i]) < 18 & hour(detection_df_1$hour_bin[i]) > 6){
            detection_df_1$day_night[i] = 'Day'
          }
        }
        detection_df_1$day_night = as.factor(detection_df_1$day_night)
        
      ### Assigning receiver height
        detection_df_1$height_of_receiver_off_bottom = NULL
        detection_df_1$height_of_receiver_off_bottom[detection_df_1$receiver == receiver_1m] = '1'
        detection_df_1$height_of_receiver_off_bottom[detection_df_1$receiver == receiver_30m] = '30'
        detection_df_1$height_of_receiver_off_bottom = as.factor(detection_df_1$height_of_receiver_off_bottom)
      
        ### Determining recovery rate based on the number of detections logged and an average of 60 transmissions per tag per hour
        detection_df_1$recovery_rate = NA
        for(i in 1:length(detection_df_1$tags_per_hour)){
          detection_df_1$recovery_rate[i] = detection_df_1$detections[i] / expected_detections$hourly_detections_per_tag[which(expected_detections$n_tags %in% detection_df_1$tags_per_hour[i])]
        }
                          
        # dim(detection_df_1)
          # 10824    20
    #### Statistical Analysis - Modeling Detection Probability ----
      
    ### Setting Model Parameters
    alpha = 0.05
    transmissions_per_hour = 60
    detection_threshold = (alpha * transmissions_per_hour)
      # 3
    
    
    ### Fitting a GAM model 
      global.gam_1 = gam( formula = detections ~ s(distance, k = 6) + s(day_night, bs = 're') + height_of_tag_off_bottom + height_of_receiver_off_bottom + s(direction, bs = "re") + s(mean_wspd, bs = "re") + s(mean_gst, bs = "re") + s(water_level, bs = "re"),
                    data = detection_df_1, keepData = TRUE, family = poisson(link = 'log'),  na.action = "na.fail")
      global_gam.summary_1 = summary(global.gam_1)
      
      ## R sqr value from GAM
        global_gam.summary_1$r.sq
          # r.sq = 0.6465411
       
      ## Checking GAM smoother terms are appropriate
        gam.check(global.gam_1)
  
    ### Using Dredge function to compare this GAM to all possible combinations of variables
        ## exporting mod.gam and data to cluster for using pdredge function
        clusterExport(clust, c('global.gam_1', 'detection_df_1'))
        mod.dredged_1 = pdredge(global.gam_1, cluster = clust)
        
        ### Comparing predictions for various canadate models from dredge that are within delta 2 of the lowest AICc
        candidate_models_1 = get.models(mod.dredged_1, subset = delta <= 2)
          # length(candidate_models_1)
            # 8
        
        #### Making Figures and Tables ----
        ## Setting up repository for results
        table_1 = data.frame(NULL)
        candidate_model_summary_1 = as.data.frame(NULL)
        candidate_model_predictions_1 = list(NULL)
        
        ## Looping through different combinations of predictor variables
        for(time_of_day in unique(detection_df_1$day_night)){
          for(receiver_height in as.numeric(levels(unique(detection_df_1$height_of_receiver_off_bottom))[unique(detection_df_1$height_of_receiver_off_bottom)])){
            for(tag_height in as.numeric(levels(unique(detection_df_1$height_of_tag_off_bottom))[unique(detection_df_1$height_of_tag_off_bottom)])){
              for(tide_direction in unique(detection_df_1$direction)){
  
        
        ## A Master list of all predictor variables from global model to be used for predicting with candidate models
                  ## Updated with values from loops
          var_predictors = list(
            day_night = time_of_day,
            height_of_tag_off_bottom = tag_height,
            height_of_receiver_off_bottom = receiver_height,
            mean_wspd = median(detection_df_1$mean_wspd),
            mean_gst = median(detection_df_1$mean_gst),
            direction = tide_direction,
            water_level = median(detection_df_1$water_level)
          )
          
          min_median_minus_se = 10000 # Something absurdly high
          max_median_plus_se = 0 # Something absurdly low
          
          min_cpdi_estimation = 1000 ## Something higher than we thing possible
          max_cpdi_estimation = 0 ## Something equal to or lower than we think possible
  
        ## Looping through candidate models
        for(i in 1:length(candidate_models_1)){
          candidate_model_id = paste('model', names(candidate_models_1)[i], sep = '_')
          candidate_model = candidate_models_1[[i]]
          
          ## Subsetting terms from master list of predictors from global model if they appear as model terms in the candidate model 
          candidate_model_terms = attr(candidate_model$terms, which = 'term.labels')
          new_mod_data = var_predictors[names(var_predictors) %in% candidate_model_terms]
          
          ### Predicting detection rates using the candidate model under evaluation
          predicted_rates = c()
          predicted_error_fits.fit = c()
          predicted_error_fits.se = c()
          for(dist in 0:1000){
            new_mod_data$distance = dist
            predicted_rates = c(predicted_rates, predict(candidate_model, newdata = new_mod_data, type = "response")[[1]])
            predicted_error_fits.fit = c(predicted_error_fits.fit, predict(candidate_model, newdata = new_mod_data, type = "response", se.fit = TRUE)$fit[[1]])
            predicted_error_fits.se = c(predicted_error_fits.se, predict(candidate_model, newdata = new_mod_data, type = "response", se.fit = TRUE)$se.fit[[1]])
          }
  
          ## Saving model predictions
          model_predictions = list()
          model_predictions$predicted_rates = predicted_rates
          model_predictions$predicted_error_fits.fit = predicted_error_fits.fit
          model_predictions$predicted_error_fits.se = predicted_error_fits.se
          candidate_model_predictions_1[[i]] = model_predictions
          
          predicted_lower_ci = predicted_error_fits.fit - predicted_error_fits.se
          predicted_upper_ci = predicted_error_fits.fit + predicted_error_fits.se
          
          ### Summary Stats
          ## Determining the average maximum detection radius
          ave_max_distance = 500+(which(abs(predicted_rates[500:1000] - detection_threshold) == min(abs(predicted_rates[500:1000] - detection_threshold))))-1 #subtract 1 because predictions start at 0 but indexing starts at 1 1
          upper_ci_max_distance = 500+(which(abs(predicted_upper_ci[500:1000] - detection_threshold) == min(abs(predicted_upper_ci[500:1000] - detection_threshold))))-1 #subtract 1 because predictions start at 0 but indexing starts at 1 1
          lower_ci_max_distance = 500+(which(abs(predicted_lower_ci[500:1000] - detection_threshold) == min(abs(predicted_lower_ci[500:1000] - detection_threshold))))-1 #subtract 1 because predictions start at 0 but indexing starts at 1 1
          
          if(upper_ci_max_distance > max_median_plus_se){
            max_median_plus_se = upper_ci_max_distance
          }
          if(lower_ci_max_distance < min_median_minus_se){
            min_median_minus_se = lower_ci_max_distance
          }
          
          ## Determining CPDI Extent
          #cpdi_extent = which.max(round(predicted_rates, digits = 0))-1 # subtract 1 because predictions start at 0, while indexing starts at 1
          ## Changing our CPDI Criterion
          max_rate_index = which.max(predicted_rates)
          max_rate_minus_se = predicted_rates[max_rate_index] - predicted_error_fits.se[max_rate_index]
          all_rates_plus_se = predicted_rates[0:max_rate_index] + predicted_error_fits.se[0:max_rate_index]
          cpdi_extent = which.min(abs(max_rate_minus_se - all_rates_plus_se)) - 1 #  subtract 1 because predictions start at 0, while indexing starts at 1
          
          
          ### Updating our min/max cpdi extents
          if(min_cpdi_estimation > cpdi_extent){
            min_cpdi_estimation = cpdi_extent
          }
          if(max_cpdi_estimation < cpdi_extent){
            max_cpdi_estimation  = cpdi_extent
          }
          
          
          ## Determining the model's R^2 value
          r.sq = summary(candidate_model)$r.sq
          
          ## Determinning deviance explained by the model
          dev.exp = summary(candidate_model)$dev.expl
          
          ## Saving summary stats to to model summary dataframe
           candidate_model_summary_1 = rbind(candidate_model_summary_1, data.frame('model_id' = candidate_model_id, 'model_terms' = paste(candidate_model_terms, sep = ' ', collapse = ' '), 'r.sq' = r.sq, 'dev.expl' = dev.exp, 'aic' = candidate_model$aic, 'ave_max_distance' = ave_max_distance,  'cpdi_extent' = cpdi_extent, 'min_with_se' = min_median_minus_se, 'max_with_se' = max_median_plus_se, stringsAsFactors = FALSE))
        }
         colnames(candidate_model_summary_1) = c('model_id', 'model_terms', 'r.sq', 'dev.expl', 'aic', 'ave_max_distance', 'cpdi_extent', 'min_with_se', 'max_with_se')
         # names(candidate_model_predictions_1) = as.character(candidate_model_summary_1$model_id)
        
      ## R sq range from candidate models 
        rsq = fivenum(as.numeric(as.character(candidate_model_summary_1$r.sq)))
          # 0.6465015 0.6465026 0.6465209 0.6465408 0.6465411
        
      ## Deviance explained by candidate models
        dev_explained = fivenum(as.numeric(as.character(candidate_model_summary_1$dev.expl)))
        
      ## Average maximum detection radius # min - SE, median, max + se
        amdr = c(min(as.numeric(as.character(candidate_model_summary_1$min_with_se))), median(as.numeric(as.character(candidate_model_summary_1$ave_max_distance))), max(as.numeric(as.character(candidate_model_summary_1$max_with_se))))
        amdr_for_table = paste(amdr[2], ' (', amdr[1], '-', amdr[3], ')', sep = "", collapse = "")
      ## CPDI extent from candidate models
        cpdi_extent = c(min_cpdi_estimation, median(candidate_model_summary_1$cpdi_extent), max_cpdi_estimation)
        cpdi_extent_for_table = paste(cpdi_extent[2], ' (', cpdi_extent[1], '-', cpdi_extent[3], ')', sep = "", collapse = "")
    ### Running Mechanistic CPDI model for receiver position
         cpdi_model = predict_cpdi_interference(bottom_depth = 300, 
                                             ave_max_detection_radius = amdr[2], 
                                             receiver_depth = (300 - receiver_height), 
                                             speed_of_sound = 1530, 
                                             max_horizontal_dist = 1000, 
                                             evaluation_interval = 1,
                                             blanking_interval = .260,
                                             plot = FALSE)
         
         ## Getting CPDI prediction from model
         cpdi_prediction = min(which(cpdi_model[(300 - tag_height), ] == 1)) - 1 # because distances are indexed starting at 1 and not 0
         ## Writing output of all this to table
         table_1 = rbind(table_1,  data.frame('Receiver Height (m)' = receiver_height, 'Tag Height (m)' = tag_height, 'Tidal Phase' = tide_direction, "Diurnal Period" = time_of_day, "GAM Estimated AMDR (m)" = amdr_for_table, "GAM Estimated CPDI (m)" = cpdi_extent_for_table, "Model Predicted CPDI (m)" = cpdi_prediction))
          
         ##### The following plot is now at the end of the script. DWRT and SWRT have been combined in one image.
         # #### Plotting Mechanistic Model Prediction and GAM Fit prediction
         # ## Setting directory to figure directory
         # setwd(figure_dir)
         # ## Subsetting data 
         # subset_prediction_data = detection_df_1[which(detection_df_1$height_of_receiver_off_bottom == receiver_height & detection_df_1$height_of_tag_off_bottom == tag_height & detection_df_1$direction == tide_direction & detection_df_1$day_night == time_of_day), ]
         # 
         # ### Creating Figure
         # image_name = paste('DWRT-RH', receiver_height, '-TH', tag_height, '-TD', tide_direction, '-TOD', time_of_day, '.png', sep = "")
         # png(filename = image_name, width = 1000, height = 500)
         # par(mfrow = c(1,2), oma = c(0, 0, 3, 0))
         # plot_heat_map(cpdi_model, save = FALSE, save_description = parameter_description_for_figure_name, main = "Mechanistic Model Prediction", compress_x_labels = 100, compress_y_labels = 50)
         # 
         # ## Creating confidence intervals from predicted error fits
         # upper_ci = predicted_rates + (predicted_error_fits.se)
         # lower_ci = predicted_rates - (predicted_error_fits.se)
         # 
         # ## Converting detections to recovery rate
         # recovery_rates = predicted_rates / 60
         # 
         # ## Plotting Empty Model
         # plot(x = subset_prediction_data$distance, y = subset_prediction_data$detections,  
         #      ylim=c(0,1), 
         #      main = 'GAM Criteron Estimation',
         #      xlab = 'Distance Between Tag and Receiver (m)', ylab = '# Transmissions Detected/Tag', 
         #      col = 'black',
         #      pch = 19,
         #      cex = .25
         #     )
         # 
         # ## Looping through model fits
         # for(r in 1:length(candidate_model_predictions_1)){
         #   predicted_recovery_rate = (candidate_model_predictions_1[[r]]$predicted_rates)/60
         #   lines(x = 0:1000 ,y = predicted_recovery_rate)
         # }
         # 
         # ## Plotting individual data points
         # points(x = subset_prediction_data$distance, 
         #        y = subset_prediction_data$recovery_rate,
         #        pch = 19, cex = .25, col = 'black')
         # 
         # ## Adding average maximum distance text
         # max_distance = 500+(which(abs(round(predicted_rates[500:1000]) - detection_threshold) == min(abs(round(predicted_rates[500:1000]) - detection_threshold))))-1 #subtract 1 because predictions start at 0, not 1
         # max_distance_lower_ci = 500+(which(abs(lower_ci[500:1000] - detection_threshold) == min(abs(lower_ci[500:1000] - detection_threshold)))) # subtract 1 because predictions start at 0, not 1
         # max_distance_upper_ci = 500+(which(abs(upper_ci[500:1000] - detection_threshold) == min(abs(upper_ci[500:1000] - detection_threshold)))) # subtract 1 because predictions start at 0, not 1
         # 
         # ## Adding AMDR to plot
         # points(x = amdr[2], y = .05, pch = 8)
         # text(x = max_distance +50, y = .15, labels = paste('AMDR \n', amdr[2]-1, 'm'), col = 'black', cex = .75)
         # 
         # ## Adding CPDI to plot
         # if(cpdi_extent[2] != 0){
         #   points(x = cpdi_extent[2], y = predicted_rates[cpdi_extent[2]], pch = 4)
         #   text(x = cpdi_extent[2], y = max(recovery_rates) + .1, labels = paste('CPDI Extent \n ', cpdi_extent[2], 'm'), cex = .75)
         # }
         # 
         # mtext('CPDI Predictions for Deep Water Ranging Experiment \n ', font = 2,outer = TRUE)
         # mtext(paste(' \n Receiver Height =', receiver_height, 'm, Tag Height =', tag_height, 'm, Tide Direction =', tide_direction, 'Diurnal Period =', time_of_day),outer = TRUE)
         # dev.off()
         # beep(1)
              }
            }
          }
        }
        
        ## Reorganizing Table 1
        table_1 = table_1[order(table_1$Receiver.Height..m., table_1$Tag.Height..m., table_1$Diurnal.Period, table_1$Tidal.Phase), ]
        
        ## Saving table_1
        setwd(results_dir)
        write.csv(table_1, 'CPDI_TABLE_1.csv', row.names = FALSE)
   
    #### Analysis Cleanup -----
      ### Saving workspace in results folder
      setwd(results_dir)
      save.image(file = 'Deep Water Range Test Results.R')
  
  ###### Experiment 2: Shallow Water Range Test ----
  
    #### Loading In and Cleaning Data Files ----
    setwd(data_dir)
  
    ### Loading Receiver Data
    receiver_data = load_receiver(filename = 'DEPLOYMENT_RECOVERY_LOG.csv')
      # dim(receiver_data)
        # 217  28
    
    ### Loading VUE data
    vue_data = load_vemco('sand_island_range_test_nov_dec_2014.csv')
      # dim(vue_data)
        # 426181     10
    
    ### Importing Windspeed Data from from NOAA Honolulu station Station ID 1612340
    # http://www.ndbc.noaa.gov/view_text_file.php?filename=oouh1h2015.txt.gz&dir=data/historical/stdmet/
    # parent: http://www.ndbc.noaa.gov/station_history.php?station=oouh1
    # Data dates range between dates 2014-01-01 00:00:00 and 2014-12-31 23:54:00
    # with the following format: 
    # YY  MM DD hh mm WDIR WSPD GST  WVHT   DPD   APD MWD   PRES  ATMP  WTMP  DEWP  VIS  TIDE
    # with the following units:   
    # yr  mo dy hr mn degT m/s  m/s     m   sec   sec degT   hPa  degC  degC  degC   mi    ft
      wind_data = as.data.frame(read.table("oouh1h2014.txt"), header = TRUE, sep = " ")
      #  dim(wind_data)
        # 85864    18
      colnames(wind_data) = c("YY",  "MM", "DD", "hh", "mm", "WDIR", "WSPD", "GST",  "WVHT",   "DPD",   "APD", "MWD",   "PRES",  "ATMP",  "WTMP",  "DEWP",  "VIS",  "TIDE")
    ### Convert dates and times to POSIXct format. (Station data time is in local (UTC) time by default)
      wind_data$datetime = as.POSIXct(paste(wind_data$YY, "-", wind_data$MM, "-", wind_data$DD, " ", wind_data$hh, ":", wind_data$mm, sep = ""), format = "%Y-%m-%d %H:%M", tz = "GMT")
    ### Convert datetime to HST
      wind_data$datetime = strftime(wind_data$datetime, format = "%Y-%m-%d %H:%M")
      wind_data$hourly = strftime(wind_data$datetime, format = "%Y-%m-%d %H")
  
    ### Grouping windspeed data into hourly means (mean_wspd)
      wind_data_grouped = group_by(wind_data, hourly)
      wspd_by_hour = summarize(wind_data_grouped, mean_wspd = mean(WSPD))
      wspd_by_hour$hourly = as.POSIXct(wspd_by_hour$hourly, format = "%Y-%m-%d %H")
  
    ### Grouping windgust data into hourly means
      gst_by_hour = summarize(wind_data_grouped, mean_gst = mean(GST))
      gst_by_hour$hourly = as.POSIXct(gst_by_hour$hourly, format = "%Y-%m-%d %H")
  
  
    ### Importing tide data from: NOAA Honolulu station, accessed 27 Jan 2017 
    # http://tidesandcurrents.noaa.gov/waterlevels.html?id=1612340&units=standard&bdate=20110930&edate=20100930&timezone=LST&datum=MLLW&interval=h&action=## Meta data available here: 
    # Data accessed for dates ranging between 2014-09-30 00:00:00 HST and 2015-09-30 13:00:00 HST
    # Data is already in HST
      tide_data = read.csv('noaa tide data 214-2015.csv')
        # dim(tide_data)
          #  8784    5
      colnames(tide_data) = c("hour_bin", "water_level", 'sigma', 'I', "L")
      tide_data$hour_bin = as.POSIXct(tide_data$hour_bin, format = "%Y-%m-%d %H:%M")
    
    ### Determining if tide was going in or out.
    ## First assume all data is headed in
      tide_data$direction = "in"
    ## Then loop through data. If current water level is less than previous water level, tide is going out.
      for(i in 2:length(tide_data$direction)){
        if(tide_data$water_level[i-1] > tide_data$water_level[i]){
          tide_data$direction[i] = "out"
        }
      }
      tide_data$direction = as.factor(tide_data$direction)
      
    ### Importing Tag Deployment Meta Data File
    tag_meta_data = read.csv('Range Test Nov 2014 - Tag Meta Data Locations and Depth.csv')
      # 18  9
    ### Converting Lat Lons from degree minutes to decimal degrees
      tag_meta_data$lat = convert_lat_lon(tag_meta_data$lat_deg, tag_meta_data$lat_min)
      tag_meta_data$lon = convert_lat_lon(tag_meta_data$lon_deg, tag_meta_data$lon_min)
  
    ### Determining distance of tag from receiver string using lon and lat waypoints from gps
      for (i in 1:length(tag_meta_data$distance)){
        tag_meta_data$distance[i] = round(1000*(lldist(point1 = c(tag_meta_data$lon[i], 
                                                            tag_meta_data$lat[i]), 
                                                 point2 = c(receiver_data$lon[receiver_data$station_name == 'Range Test - Nov 2014 - Sand Island 1m'][1], 
                                                            receiver_data$lat[receiver_data$station_name == 'Range Test - Nov 2014 - Sand Island 1m'][1]))))
      }
     unique(tag_meta_data$distance)
        # 0   75  150  300  600 1200
  
    ### Removing detections prior to experiment (occurred durring experiment set up on boat) and determining an end to the experiment.
      start_date = as.POSIXct("2014-11-21 10:15:00 HST")
      end_date = as.POSIXct("2014-12-02 07:22:00 HST")
      vue_data = vue_data[which(vue_data$datetime >= start_date & vue_data$datetime < end_date), ]
     
    ### Rounding off of final partial hour to account for transmissions received while gear was recovered
      vue_data = vue_data[which(vue_data$datetime >= ceiling_date(min(vue_data$datetime), unit = 'hour') & vue_data$datetime < floor_date(max(vue_data$datetime), unit = 'hour')), ]
      range(vue_data$datetime)
        # "2015-03-17 19:00:23 HST" "2015-03-25 10:58:53 HST"
        # dim(vue_data)
          # 423384     10
  
  
    ### Grouping receivers by height condition (m of sea floor)
      rec_1m  = c(103911, 110298, 110308)
      rec_7.5m = c(102200, 123736, 123732) 
      rec_15m = c(123733, 104543, 123734)
  
    ## How many total detections occurred durring the course of the experiment?
      detections_all = dim(vue_data)[1]
        # 423384
  
    ### Removing tags not in experiment. This could be a tagged fish that swam by or a false detection
        tags_0m    = c(18255, 18257, 18256)
        tags_75m  = c(18260, 18273, 18270)
        tags_150m  = c(18265, 18259, 18269)
        tags_300m  = c(18262, 18258, 18271)
        tags_600m  = c(18261, 18275, 18274)
        tags_1200m = c(18272, 18268, 18263)
        tag_ids = c(tags_0m, tags_75m, tags_150m, tags_300m, tags_600m, tags_1200m) 
        vue_data = vue_data[as.numeric(levels(vue_data$tag_id))[vue_data$tag_id]
                            %in% tag_ids, ] 
    
      ## How many detections were from tags that were part of the experiment?
        detections_experiment = dim(vue_data)[1]
          # Total detections of experimental tags = 421342
        detections_experiment.percent = detections_experiment / detections_all * 100
          # 99.5177 %
    
      ## How many detections were from tags that were not part of the experiment?
        detections_other = detections_all - detections_experiment
          # Total detections from other tags = 2042
        detections_other.percent = detections_other / detections_all * 100
          # 0.4823045 %
    
    ## Assigning Lon Lat positions to detections
    vue_data = clean_vue_lat_lon(vue_data_df = vue_data, 
                                 receiver_data_df = receiver_data)
    
    ## Getting estimated distances between receivers and tags based on deployment coordinates
    for (i in 1:length(tag_meta_data$lon)){
      tag_meta_data$distance[i] = lldist(point1 = c(unique(vue_data$lon)[1], 
                                                    unique(vue_data$lat)[1]),
                                         point2 = c((tag_meta_data$lon)[i],
                                                    (tag_meta_data$lat)[i])) *
        1000 # m/km
    }
    round(unique(tag_meta_data$distance))
    # Distances are: 0   75  150  300  600 1200
  
    #### Hourly Detection Data Analysis ----
    ### Binning data by tag, receiver, and hour bin
    ## First assign each datum to an date and hour bin
      vue_data$hour_bin = floor_date(vue_data$datetime, unit = "hour")
  
    ## Then loop through combinations of date/hour, tag id, and receiver to get the number of detections from a tag at a receiver each hour
      detection_df_2 = data.frame()
      registerDoParallel(cores = n_cores)
      detection_df_2 = foreach(b = c(1:length(unique(vue_data$hour_bin))), .combine = rbind) %:% 
        foreach(i = c(1:length(unique(vue_data$tag_id))), .combine = rbind) %:% 
        foreach(r = c(1:length(unique(vue_data$receiver))), .combine = rbind) %dopar%{
          filtered = filter(vue_data, hour_bin == unique(vue_data$hour_bin)[b], 
                            tag_id == unique(vue_data$tag_id)[i],
                            receiver == unique(vue_data$receiver)[r])
          write_line = cbind(as.character(unique(vue_data$hour_bin))[b],
                             as.character(unique(vue_data$tag_id))[i], 
                             as.character(unique(vue_data$receiver))[r], 
                             dim(filtered)[1])
          return(write_line)
        }
      ## Cleaning detection_df_2
      detection_df_2 = as.data.frame(detection_df_2)
      colnames(detection_df_2) = c('hour_bin', 'tag_id', 'receiver', 'detections')
        # dim(detection_df_2)
          # 32760     4
      detection_df_2$detections = as.numeric(as.character(detection_df_2$detections))
      detection_df_2$day = as.factor(strftime(detection_df_2$hour_bin, format = "%Y-%m-%d"))
      detection_df_2$time = as.factor(strftime(detection_df_2$hour_bin, format = "%H:%M:%S"))
      detection_df_2$height_of_tag_off_bottom = 7.5
      detection_df_2$height_of_tag_off_bottom = as.factor(detection_df_2$height_of_tag_off_bottom)
      detection_df_2$hour_bin = as.POSIXct(detection_df_2$hour_bin, format = "%Y-%m-%d %H:%M:%S")
    
    
    ### Combining detection dataframe with distance, bottom depth, height condition, lat, and lon tag_meta_data
      detection_df_2 = merge(detection_df_2, tag_meta_data, 
                          by.x = colnames(detection_df_2) == 'tag_id', 
                          by.y = colnames(tag_meta_data) == 'tag_id')  
  
    ## Combining detection data with wind speed, and wind gust data
      detection_df_2 = merge(detection_df_2, wspd_by_hour,
                           by.x = colnames(detection_df_2) == 'hour_bin',
                           by.y = colnames(wspd_by_hour) == 'hourly')
      detection_df_2 = merge(detection_df_2, gst_by_hour,
                           by.x = colnames(detection_df_2) == 'hour_bin',
                           by.y = colnames(gst_by_hour) == 'hourly')
      detection_df_2 = merge(detection_df_2, tide_data,
                           by.x = colnames(detection_df_2) == 'hour_bin',
                           by.y = colnames(tide_data) == 'hour_bin')
      
      detection_df_2$hour_bin = as.factor(detection_df_2$hour_bin)
      
    ### Assigning Receivers to a height condition (height of sea floor)
      detection_df_2$height_of_receiver_off_bottom = NA
      detection_df_2$height_of_receiver_off_bottom[which(detection_df_2$receiver %in% rec_1m)] = 1
      detection_df_2$height_of_receiver_off_bottom[which(detection_df_2$receiver %in% rec_7.5m)] = 7.5
      detection_df_2$height_of_receiver_off_bottom[which(detection_df_2$receiver %in% rec_15m)] = 15
      detection_df_2$height_of_receiver_off_bottom = as.factor(detection_df_2$height_of_receiver_off_bottom)
  
    ### Assigning diurnal period to detections based around 12 hour time bins. Before 6am or after 6pm considered night, between 6am and 6pm considered day
      detection_df_2$day_night = 'night'
      for(i in 1:length(detection_df_2$hour_bin)){
        if(hour(detection_df_2$hour_bin[i]) < 18 & hour(detection_df_2$hour_bin[i]) > 6){
          detection_df_2$day_night[i] = 'day'
        }
      }
      detection_df_2$day_night = as.factor(detection_df_2$day_night)
      
    ### Determining recovery rate based on the number of detections logged and an average of 60 transmissions per tag per hour
      detection_df_2$recovery_rate = detection_df_2$detections / 60
      # dim(detection_df_2)
        # 32760    27
      
      range(aggregate(detection_df_2$tag_id[detection_df_2$detections != 0], by = list(detection_df_2$hour_bin[detection_df_2$detections != 0]), FUN = uniqueN)$'x')
        # 10 - 13
      
    #### Statistical Analysis - Modeling Detection Probability ----
  
    ### Setting Model Parameters
      alpha = 0.05
      transmissions_per_hour = 60
      detection_threshold = (alpha * transmissions_per_hour)
        # 3
      
    ### Fitting a GAM model 
      global.gam_2 = gam(formula = detections ~ s(distance, k = 5) + height_of_receiver_off_bottom + s(direction, bs = "re") + s(mean_wspd, bs = "re") + s(mean_gst, bs = "re") + s(water_level, bs = "re") + s(day_night, bs = "re"),  # + s(receiver, bs = "re"),
                    data = detection_df_2, family = poisson(link = 'log'),  na.action = "na.fail")
      gam.summary_2 = summary(global.gam_2)
      
      ## R sqr value from gam
        gam.summary_2$r.sq
          # r.sq = 0.6870176
  
      ## Checking GAM smoother terms are appropriate
        gam.check(global.gam_2)
  
      ### Using Dredge function to compare this GAM to all possible combinations of variables
        ## Exporting mod.gam and data to cluster for using pdredge function
          clusterExport(clust, c('global.gam_2', 'detection_df_2'))
          mod.dredged_2 = pdredge(global.gam_2, cluster = clust)
          
          ### Comparing predictions for various canadate models from dredge that are within delta 2 of the lowest AICc
          candidate_models_2 = get.models(mod.dredged_2, subset = delta <= 2)
          table_2 = data.frame(NULL)
          
          ## Looping through different combinations of predictor variables
          for(time_of_day in unique(detection_df_2$day_night)){
            for(receiver_height in as.numeric(levels(unique(detection_df_2$height_of_receiver_off_bottom))[unique(detection_df_2$height_of_receiver_off_bottom)])){
                for(tide_direction in unique(detection_df_2$direction)){
          
          ## A Master list of all predictor variables from global model to be used for predicting with candidate models
          var_predictors = list(
            day_night = time_of_day,
            height_of_receiver_off_bottom = receiver_height,
            mean_wspd = median(detection_df_2$mean_wspd),
            mean_gst = median(detection_df_2$mean_gst),
            direction = tide_direction,
            water_level = median(detection_df_2$water_level)
          )
        
          ## Setting up output dataframes for summary statistics and predictions from candidate models
          candidate_model_summary_2 = as.data.frame(NULL)
          candidate_model_predictions_2 = list()
          
          ## For calculating confidence intervals
          min_median_minus_se = 10000 # Something absurdly high
          max_median_plus_se = 0 # Something absurdly low
          
          min_cpdi_estimation = 120000 # Something absurdly high
          max_cpdi_estimation = 0 # Something equal to or lower than our lowest value
          
          ## Looping through candidate models
          for(i in 1:length(candidate_models_2)){
            candidate_model_id = paste('model', names(candidate_models_2)[i], sep = '_')
            candidate_model = candidate_models_2[[i]]
            
            ## Subsetting terms from master list of predictors from global model if they appear as model terms in the candidate model 
            candidate_model_terms = attr(candidate_model$terms, which = 'term.labels')
            new_mod_data = var_predictors[names(var_predictors) %in% candidate_model_terms]
            
            ### Predicting detection rates using the candidate model under evaluation
            predicted_rates = c()
            predicted_error_fits.fit = c()
            predicted_error_fits.se = c()
            for(dist in 0:1200){
              new_mod_data$distance = dist
              predicted_rates = c(predicted_rates, predict(candidate_model, newdata = new_mod_data, type = "response")[[1]])
              predicted_error_fits.fit = c(predicted_error_fits.fit, predict(candidate_model, newdata = new_mod_data, type = "response", se.fit = TRUE)$fit[[1]])
              predicted_error_fits.se = c(predicted_error_fits.se, predict(candidate_model, newdata = new_mod_data, type = "response", se.fit = TRUE)$se.fit[[1]])
            }
            
            ## Creating confidence intervals from predicted error fits
            upper_ci = predicted_rates + (predicted_error_fits.se)
            lower_ci = predicted_rates - (predicted_error_fits.se)
            
            ## Saving model predictions
            model_predictions = list()
            model_predictions$predicted_rates = predicted_rates
            model_predictions$predicted_error_fits.fit = predicted_error_fits.fit
            model_predictions$predicted_error_fits.se = predicted_error_fits.se
            candidate_model_predictions_2[[i]] = model_predictions
            
            ### Summary Stats
            ## Determining the average maximum detection radius
            ave_max_distance = which(abs(predicted_rates - detection_threshold) == min(abs(predicted_rates - detection_threshold)))-1 #subtract 1 because predictions start at 0, not 1
            upper_ci_max_distance = which(abs(upper_ci - detection_threshold) == min(abs(upper_ci - detection_threshold)))-1 #subtract 1 because predictions start at 0 but indexing starts at 1 1
            lower_ci_max_distance = which(abs(lower_ci - detection_threshold) == min(abs(lower_ci - detection_threshold)))-1 #subtract 1 because predictions start at 0 but indexing starts at 1 1
            
            if(upper_ci_max_distance > max_median_plus_se){
              max_median_plus_se = upper_ci_max_distance
            }
            if(lower_ci_max_distance < min_median_minus_se){
              min_median_minus_se = lower_ci_max_distance
            }
            
            ## Determining CPDI Extent
            #cpdi_extent = which.max(round(predicted_rates, digits = 0))-1 # subtract 1 because predictions start at 0, while indexing starts at 1
            ## Changing our CPDI Criterion
            max_rate_index = which.max(predicted_rates)
            max_rate_minus_se = predicted_rates[max_rate_index] - predicted_error_fits.se[max_rate_index]
            all_rates_plus_se = predicted_rates[0:max_rate_index] + predicted_error_fits.se[0:max_rate_index]
            cpdi_extent = which.min(abs(max_rate_minus_se - all_rates_plus_se)) - 1 #  subtract 1 because predictions start at 0, while indexing starts at 1
            
            ## Updating CPDI bounds
            if(min_cpdi_estimation > cpdi_extent){
              min_cpdi_estimation = cpdi_extent
            }
            if(max_cpdi_estimation < cpdi_extent){
              max_cpdi_estimation = cpdi_extent
            }
            
            ## Determining the model's R^2 value
            r.sq = summary(candidate_model)$r.sq
            
            ## Determining the deviance explained for the model
            dev.explained = summary(candidate_model)$dev.expl
            
            ## Saving summary stats to to model summary dataframe
            candidate_model_summary_2 = rbind(candidate_model_summary_2, cbind(as.character(candidate_model_id), paste(candidate_model_terms, sep = ' ', collapse = ' '), r.sq, dev.explained, candidate_model$aic, ave_max_distance, cpdi_extent))
          }
          colnames(candidate_model_summary_2) = c('model_id', 'model_terms', 'r.sq', 'dev.expl', 'aic', 'ave_max_distance', 'cpdi_extent')
          names(candidate_model_predictions_2) = as.character(candidate_model_summary_2$model_id)
  
          
      ### R sq range from candidate models 
        rsq = fivenum(as.numeric(as.character(candidate_model_summary_2$r.sq)))
        #  0.6883606 0.6883628 0.6883662 0.6883737 0.6883817
      
      ### Deviance explained from candidate models
        dev_explained = fivenum(as.numeric(as.character(candidate_model_summary_2$dev.expl)))
        # 0.7266961 0.7266961 0.7266962 0.7266962 0.7266963
        
      ### Average maximum detection radius
        amdr = c(lower_ci_max_distance, median(as.numeric(as.character(candidate_model_summary_2$ave_max_distance))), upper_ci_max_distance)
        amdr_for_table = paste(amdr[2], ' (', amdr[1], '-', amdr[3], ')', sep = "",collapse = "")
      ### CPDI extent from candidate models
        cpdi_extent = c(min_cpdi_estimation, median(as.numeric(as.character(candidate_model_summary_2$cpdi_extent))), max_cpdi_estimation)
        cpdi_extent_for_table = paste(cpdi_extent[2], ' (', cpdi_extent[1], '-', cpdi_extent[3], ')', sep = "", collapse = "")
      ### Running Mechanistic CPDI model for receiver position
        cpdi_model = predict_cpdi_interference(bottom_depth = 25, 
                                                 ave_max_detection_radius = amdr[3], 
                                                 receiver_depth = (25 - receiver_height), 
                                                 speed_of_sound = 1530, 
                                                 max_horizontal_dist = 1200, 
                                                 blanking_interval = .260,
                                                 evaluation_interval = 1)
        
        ## Getting CPDI prediction from model
        cpdi_prediction = min(which(cpdi_model[(25 - receiver_height), ] == 1))-1
        
        ## Writing output of all this to table
        table_2 = rbind(table_2, data.frame('Receiver Height' = receiver_height, 'Tidal Phase' = tide_direction, 'Diurnal Period' = time_of_day, 'GAM Estimated AMDR (m)' = amdr_for_table, 'GAM Estimated CPDI (m)' = cpdi_extent_for_table, 'Model Predicted CPDI (m)' = cpdi_prediction))
  
        ##### The following code has been depreciated. Figure is now made at end of script and part of figure comparing DWRT to SWRT
        # ### Plotting Mechanistic Model Prediction and GAM Fit prediction Figures
        # ## Setting directory to output a figure of model
        # setwd(figure_dir)
        # 
        # ## Subsetting data to be only data used in fit
        # subset_prediction_data = detection_df_2[which(detection_df_2$height_off_bottom == receiver_height & detection_df_2$direction == tide_direction & detection_df_2$day_night == time_of_day), ]
        # 
        # ## Creating Figure
        # image_name = paste('SWRT-RH', receiver_height, tide_direction, '-TOD', time_of_day, '.png', sep = "")
        # png(filename = image_name, width = 1000, height = 500)
        # par(mfrow = c(1,2), oma = c(0, 0, 3, 0))
        # plot_heat_map(cpdi_model, save = FALSE, save_description = parameter_description_for_figure_name, main = "Mechanistic Model Prediction", compress_x_labels = 100, compress_y_labels = 5)
        # 
        # ### Creating confidence intervals from predicted error fits
        # upper_ci = predicted_rates + (2 * predicted_error_fits.se)
        # lower_ci = predicted_rates - (2 * predicted_error_fits.se)
        # 
        # ### Converting detections to recovery rate
        # recovery_rates = predicted_rates / 60
        # 
        # ## Plotting Model Fit
        # plot(x = subset_prediction_data$distance, y = subset_prediction_data$detections,  
        #      ylim=c(0,1), xlim = c(0, 1200),
        #      main = 'GAM Criteron Estimation',
        #      xlab = 'Distance Between Tag and Receiver (m)', ylab = '# Transmissions Detected/Tag', 
        #      col = 'black',
        #      pch = 19,
        #      cex = .25
        # )
        # 
        # ## Looping through model fits
        # for(r in 1:length(candidate_model_predictions_2)){
        #   predicted_recovery_rate = (candidate_model_predictions_2[[r]]$predicted_rates)/60
        #   lines(x = 0:1200 ,y = predicted_recovery_rate)
        # }
        # 
        # ## Plotting individual data points
        # points(x = subset_prediction_data$distance, 
        #        y = subset_prediction_data$recovery_rate,
        #        pch = 19, cex = .25, col = 'black')
        # 
        # ## Adding average maximum distance text
        # max_distance = (which(abs(predicted_rates - detection_threshold) == min(abs(predicted_rates - detection_threshold))))-1 #subtract 1 because predictions start at 0, not 1
        # max_distance_lower_ci = (which(abs(lower_ci - detection_threshold) == min(abs(lower_ci - detection_threshold))))-1 # subtract 1 because predictions start at 0, not 1
        # max_distance_upper_ci = (which(abs(upper_ci - detection_threshold) == min(abs(upper_ci - detection_threshold))))-1 # subtract 1 because predictions start at 0, not 1
        # 
        # ## Adding AMDR to plot
        # points(x = amdr[3], y = .05, pch = 8)
        # text(x = max_distance +50, y = .15, labels = paste('AMDR \n', max_distance-1, 'm'), col = 'black', cex = .75)
        # 
        # ## Adding CPDI to plot
        # if(cpdi_extent[3] != 0){
        #   points(x = cpdi_extent[3], predicted_rates[cpdi_extent[3]], pch = 4)
        #   text(x = cpdi_extent[3], y = max(recovery_rates) + .1, labels = paste('CPDI Extent \n ', cpdi_extent[3], 'm'), cex = .75)
        # }
        # 
        # mtext('CPDI Predictions for Shallow Water Ranging Experiment \n ', font = 2,outer = TRUE)
        # mtext(paste(' \n Receiver Height =', receiver_height, 'm, Tag Height =', tag_height, 'm, Tide Direction =', tide_direction, 'Diurnal Period =', time_of_day),outer = TRUE)
        # dev.off()
        # 
        # 
        # beep(1)
                }
            }
          }
    
          ## Reorganizing the results of table 2
          table_2 = table_2[order(table_2$Receiver.Height, table_2$Diurnal.Period, table_2$Tidal.Phase), ]
  
          
          
    #### Analysis Cleanup -----
      setwd(results_dir)
      ### Saving table_2
      write.csv(table_2, 'CPDI_TABLE_2.csv', row.names = FALSE)
      ### Saving workspace in results folder
      setwd(results_dir)
      save.image(file = 'Shallow Water Range Test Results.R')
      
      
  ###### Experiment 3: Depth Dependent Validation Experiment ----
      
    #### Predicting CPDI / No CPDI for two different depth scenarios
     rec_50m =  predict_cpdi_interference(bottom_depth = 50,
                                  ave_max_detection_radius = 843,
                                  receiver_depth = 49,
                                  speed_of_sound = 1530,
                                  max_horizontal_dist = 50,
                                  evaluation_interval = 1)
      # plot_heat_map(rec_50m, save = TRUE, save_description = 'Experiment_3_depth_dependent_50m.png')
      
      
     rec_215m =  predict_cpdi_interference(bottom_depth = 215,
                                  ave_max_detection_radius = 846,
                                  receiver_depth = 214,
                                  speed_of_sound = 1530,
                                  max_horizontal_dist = 50,
                                  evaluation_interval = 1)
     # plot_heat_map(rec_215m, save = TRUE, save_description = 'Experiment_3_depth_dependent_215m.png')
     
     
     ## Because receiver for 215 m missed the mark, making sure predictions are the same for 212 m
     rec_212m =  predict_cpdi_interference(bottom_depth = 212,
                                             ave_max_detection_radius = 846,
                                             receiver_depth = 211,
                                             speed_of_sound = 1530,
                                             max_horizontal_dist = 50,
                                             evaluation_interval = 1)
     # plot_heat_map(rec_212m, save = TRUE, save_description = 'Experiment_3_depth_dependent_212m.png')
     
      
   ##### Evaluating experiment results
    #### Loading In and Cleaning Data Files ----
      setwd(data_dir)
    ### Loading Receiver Data
      receiver_data = load_receiver(filename = 'DEPLOYMENT_RECOVERY_LOG.csv')
      ### Loading VUE data
        vue_data = load_vemco(filename = 'VUE_Export_CPDI_Depth_Validation_Experiment_March_2015.csv')
    ### Loading in Receiver Meta Data Files
      metadata_212m = read.csv("VUE_Export_Rec_102202_metadata.csv")
      metadata_50m = read.csv("VUE_Export_Rec_110317_metadata.csv")
  
    
    ### Assigning receiver SNs and tag IDs to respective conditions
      receiver_50m = 110317 # serial number of receiver in 50 m depth
      receiver_212m = 102202 # serial number of receiver in 212 m depth
      tag_50m = 18275
      tag_212m = 18274
      tag_ids = c(tag_50m, tag_212m)
  
    ### Removing detections prior to deployment (occurred durring set up on boat)
      start_date = as.POSIXct('2015-03-17 17:40:00', tz = "HST")
      vue_data = vue_data[vue_data$datetime >= start_date, ]
      
    ### Removing detections after end of experiment
      ## Determining an end to the experiment -
      ## Receiver that was found floating off airport could have broken free at anytime. 
      ## Because tag attached to that receiver (18275) was detected at the 212 m receiver,
      ## The last time this tag was detected is a stand in for when the receiver broke free.
        end_date = min(vue_data$datetime[which(vue_data$receiver == receiver_212m & vue_data$tag_id == tag_50m)])
        vue_data = vue_data[vue_data$datetime < end_date, ]
  
      ## Removing detections before start of first full hour of test and after last full hour of test
        vue_data = vue_data[vue_data$datetime >= ceiling_date(min(vue_data$datetime), 
                                                              unit = "hour"), ]
        vue_data = vue_data[vue_data$datetime <= floor_date(max(vue_data$datetime), 
                                                            unit = "hour"), ]
        
    ## How many total detections of any tag occurred during the course of the experiment?  
      detections_all = dim(vue_data)[1]
        # 12382
  
    ### Removing detections of tags not associated with current experiment
      vue_data = vue_data[vue_data$tag_id %in% tag_ids, ]
    
    ## How many detections of only experiment tags?
      detections_experiment = dim(vue_data)[1]
        # 12381
      detections_experiment.percent = detections_experiment / detections_all * 100
        # 99.99192 %
  
    ## How many detections were from non-experiment tags?
      detections_other = detections_all - detections_experiment
        # 1
      detections_other.percent = detections_other / detections_all * 100
        # 0.00807624
  
  
    #### Hourly Detection Data Analysis----
    ### Binning Data by Hour
    ## Separating data by tag, receiver, and hour bin
    # First assign each datum to an date and hour bin
      vue_data$hour_bin = floor_date(vue_data$datetime, unit = "hour")
    # Then loop through combinations of date/hour, tag id, and receiver to get the number of detections from a tag at a receiver each hour
      detection_df_3 = data.frame()
      registerDoParallel(cores = n_cores)
      detection_df_3 = foreach(b = c(1:length(unique(vue_data$hour_bin))), .combine = rbind) %:% 
        foreach(r = c(1:length(unique(vue_data$receiver))), .combine = rbind) %dopar%{
          filtered = filter(vue_data, hour_bin == unique(vue_data$hour_bin)[b], 
                            receiver == unique(vue_data$receiver)[r])
          write_line = cbind(as.character(unique(vue_data$hour_bin))[b],
                             as.character(unique(vue_data$receiver))[r], 
                             dim(filtered)[1])
          return(write_line)
        }
        detection_df_3 = as.data.frame(detection_df_3)
        colnames(detection_df_3) = c('hour_bin', 'receiver', 'detections')
        detection_df_3$detections = as.numeric(as.character(detection_df_3$detections))
    ## Assigning depth condition based on receiver sn
      detection_df_3$depth = NA
      detection_df_3$depth[detection_df_3$receiver == receiver_50m] = '50 m'
      detection_df_3$depth[detection_df_3$receiver == receiver_212m] = '212 m'
  
    #### Statistical Analysis - Comparing Hourly Detections ----
  
    #### Hourly Receiver Comparison
    ### ANOVA - Detection rates of receiver not lost for period before and after first receiver lost – p-value
    
    ### Comparing Number of Detections
      ## Shapiro Wilk test that distribution of each case are non-parametric (p-value)
        shapiro.test(detection_df_3$detections[detection_df_3$depth == '50 m']) 
          # p = 0.02473, p < 0.05 reject null hypothesis that distribution is normal
          # Conclusion: Distribution is not normal
        shapiro.test(detection_df_3$detections[detection_df_3$depth == '212 m']) 
          # p < 2.482e-13 * 10^-13, p < 0.05 reject null hypothesis that distribution is normal
          # Conclusion: Distribution is not normal  
  
        ## Wilcox test - Selected to compare two groups (50m and 500m) non-parametrically
        wilcox.test(x = detection_df_3$detections[detection_df_3$depth == '50 m'], 
                    y = detection_df_3$detections[detection_df_3$depth ==  '212 m'],
                    paired = TRUE)
          # p < 2.2 * 10^-16, p < 0.05 reject null hypothesis that difference between groups is zero
          # Conclusion Wilcoxon Sign rank test shows significant difference between two receiver conditions (p-value)
        
    # Ratio for Mean hourly detection rates for 50m receiver and 212 m receiver 
        mean_detections_50m  = mean(detection_df_3$detections[detection_df_3$depth == '50 m'])
          # 56.56989
        mean_detections_212m = mean(detection_df_3$detections[detection_df_3$depth == '212 m'])
          # 9.994624
        mean_detections_50m/mean_detections_212m
          # Receiver at 50 m has 5.660032 x as many detections as receiver in 212 m
        
    # Number of hours where deeper receiver detection rates exceeds shallow receiver 
      total_hour_bins = length(unique(detection_df_3$hour_bin))
      times_212m_exceeds_50m = 0
      times_212m_equal_to_50m = 0
        for(hour in unique(detection_df_3$hour_bin)){
          subset_detection_df_3 = filter(detection_df_3, hour_bin == hour)
          if(subset_detection_df_3$detections[subset_detection_df_3$depth == '212 m'] >
            subset_detection_df_3$detections[subset_detection_df_3$depth == '50 m']){
            times_212m_exceeds_50m = times_212m_exceeds_50m + 1
          }else if(subset_detection_df_3$detections[subset_detection_df_3$depth == '212 m'] ==
                   subset_detection_df_3$detections[subset_detection_df_3$depth == '50 m']){
            times_212m_equal_to_50m = times_212m_equal_to_50m + 1
          }
        }
      times_212m_exceeds_50m
        # 0
      times_212m_equal_to_50m
        # 0
      
      
    #### Daily Meta Data Analysis ----
      start_date = as.POSIXct('2015-3-19')
      
    #### Cleaning receiver meta data 
      ## Receiver SN 102202
        colnames(metadata_212m) = c('datetime', 'receiver', 'description', 'data', 'units')
        metadata_212m$datetime = as.POSIXct(metadata_212m$datetime)
        metadata_212m = metadata_212m[which(metadata_212m$datetime >= start_date & metadata_212m$datetime < floor_date(end_date, unit = "day")), ]
        metadata_212m$data = as.numeric(as.character(metadata_212m$data))
      ## Receiver SN receiver_50m
        colnames(metadata_50m) = c('datetime', 'receiver', 'description', 'data', 'units')
        metadata_50m$datetime = as.POSIXct(metadata_50m$datetime)
        metadata_50m = metadata_50m[which(metadata_50m$datetime >= start_date & metadata_50m$datetime < floor_date(end_date, unit = "day")), ]
        metadata_50m$data = as.numeric(as.character(metadata_50m$data))
      
    ### Creating daily dataframe with criteria: Pings, Syncs, Rejects, Detections
      metadata_df = data.frame()
      for(i in 1:length(unique(metadata_212m$datetime))){
        ## Subset metadata log for receiver receiver_212m and write to dataframe 
          subset_metadata = metadata_212m[metadata_212m$datetime == unique(metadata_212m$datetime)[i], ]
          metadata_df = rbind(metadata_df, cbind(as.character(subset_metadata$datetime[1]), 'receiver_212m', subset_metadata$data[subset_metadata$description == 'Daily Pings on 69 kHz'], subset_metadata$data[subset_metadata$description == 'Daily Syncs on 69 kHz'], subset_metadata$data[subset_metadata$description == 'Daily Detections on A69-1601'], subset_metadata$data[subset_metadata$description == 'Daily Rejects on 69 kHz']))
        ## Subset metadata log for receiver receiver_50m and write to dataframe
          subset_metadata = metadata_50m[metadata_50m$datetime == unique(metadata_50m$datetime)[i], ]
          metadata_df = rbind(metadata_df, cbind(as.character(subset_metadata$datetime[1]), 'receiver_50m', subset_metadata$data[subset_metadata$description == 'Daily Pings on 69 kHz'], subset_metadata$data[subset_metadata$description == 'Daily Syncs on 69 kHz'], subset_metadata$data[subset_metadata$description == 'Daily Detections on A69-1601'], subset_metadata$data[subset_metadata$description == 'Daily Rejects on 69 kHz']))
      }
      colnames(metadata_df) = c('date', 'receiver', 'pings', 'syncs', 'detections', 'rejections')
      
      ## Adding the following metrics based on Simpfendorfer et al 2008: Code Detection Efficiency (cde), rejection coefficient (rc), noise quotient (nq)
        metadata_df$cde = as.numeric(as.character(metadata_df$detections)) / as.numeric(as.character(metadata_df$syncs))
        metadata_df$rc  = as.numeric(as.character(metadata_df$rejections)) / as.numeric(as.character(metadata_df$syncs))
        metadata_df$nq  = as.numeric(as.character(metadata_df$pings)) - (as.numeric(as.character(metadata_df$syncs)) * 7)
        
    #### Statistical Analysis - Comparing Daily Detection Metrics ----  
  
      ### Comparing Code Detection Efficiency
      ## Shapiro Wilk Test of Code Detection Efficiency
        shapiro.test(metadata_df$cde[metadata_df$receiver == 'receiver_212m'])
          # p-value = 0.7799
        shapiro.test(metadata_df$cde[metadata_df$receiver == 'receiver_50m'])
          # p-value = 0.002853
      
      ## Paired Wilcoxon test for code detection efficiency metric for each receiver
        wilcox.test(x = metadata_df$cde[metadata_df$receiver == 'receiver_212m'],
                    y = metadata_df$cde[metadata_df$receiver == 'receiver_50m'], 
                    paired = TRUE)
            # p-value = 0.03125
        median(x = metadata_df$cde[metadata_df$receiver == 'receiver_212m'])
          # 0.0865077
        median(x = metadata_df$cde[metadata_df$receiver == 'receiver_50m'])
          # 1
      
      ### Comparing Rejection Coefficient  
      ## Shapiro Wilk Test of Rejection Coefficient
        ### NOTE: Commmented out as cannot shapiro.test when all values are the same. rejection coefficent is 0 for all hour bins, therefore throws error when sourced 
        # shapiro.test(metadata_df$rc[metadata_df$receiver == 'receiver_50m'])
          # Rejection Coefficient = 0. Cannot do a shapiro test.
        shapiro.test(metadata_df$rc[metadata_df$receiver == 'receiver_212m'])
          # p-value = 0.9427
        
      ## Paired Wilcoxon test for rejection coefficient metric for each receiver.
        wilcox.test(x = metadata_df$rc[metadata_df$receiver == 'receiver_212m'],
                    y = metadata_df$rc[metadata_df$receiver == 'receiver_50m'], 
                    paired = TRUE)
          # p-value = 0.03125
        median(metadata_df$rc[metadata_df$receiver == 'receiver_212m'])
          # 0.01375137
        median(metadata_df$rc[metadata_df$receiver == 'receiver_50m'])
          # 0.00
      
      ### Comparing Noise Quotient
      ## Shapiro Wilk Test of Noise Quotient
        shapiro.test(metadata_df$nq[metadata_df$receiver == 'receiver_50m'])
          # p-value = 0.002315
        shapiro.test(metadata_df$nq[metadata_df$receiver == 'receiver_212m'])
          # p-value = 0.25
      ## Paired wilcoxon test for noise quotient metric for each receiver 
        wilcox.test(x = metadata_df$nq[metadata_df$receiver == 'receiver_212m'],
                    y = metadata_df$nq[metadata_df$receiver == 'receiver_50m'], 
                    paired = TRUE)
            # p-value = 0.03125
        median(metadata_df$nq[metadata_df$receiver == 'receiver_212m'])
          # -13793.5
        median(metadata_df$nq[metadata_df$receiver == 'receiver_50m'])
          # 1356
      
    #### Analysis Cleanup -----
      ### Saving workspace in results folderb
        setwd(results_dir)
        save.image(file = 'Depth Validation Test Results.R')
  
        
  ###### Experiment 4: Depth and Distance Validation Experiment ----
        
    #### Loading In and Cleaning Data Files ############# 
      setwd(data_dir)
    
    ### Loading Receiver Data 
      receiver_data = load_receiver(filename = 'DEPLOYMENT_RECOVERY_LOG.csv')
      
    ### Loading VUE data 
      vue_data = load_vemco(filename = 'VUE_Export_Depth_and_Distance_Validation_Experiment_June_2015.csv')
      
    ### Loading receiver metadata logs
      metadata_receiver_60m = read.csv('VUE_Export_Rec_123739_metadata.csv')
      metadata_receiver_508m = read.csv('VUE_Export_Rec_102199_metadata.csv')
        
    ### Cleaning Vue Data
      ## Assigning receiver SNs and tag IDs to respective conditions
        receiver_60m = 123739 # serial number of receiver in 50 m depth
        receiver_508m = 102199 # serial number of receiver in 212 m depth
        tag_ids = c(46915, 46912, 46913)
        
      ## Removing detections prior to deployment (occurred durring set up on boat) and after equipment recovered
        start_date = as.POSIXct('2015-05-24 19:12:25', tz = "HST")
        end_date = as.POSIXct('2015-05-30 09:46:00', tz = 'HST')
        vue_data = vue_data[which(vue_data$datetime >= start_date & vue_data$datetime < end_date), ]
  
      ## Removing detections from partial hour bins, aka hour equipment was deployed and hour it was recovered before first full hour of test
        vue_data = vue_data[which(vue_data$datetime >= ceiling_date(min(vue_data$datetime), unit = "hour") & 
                                    vue_data$datetime <= floor_date(max(vue_data$datetime), unit = "hour")), ]
      
      ## How many total tag detections were there?
        detections_all = dim(vue_data)[1]
          # Total detections of all tags = 1699
      
      ## Removing detections of tags not associated with current experiment
        vue_data = vue_data[vue_data$tag_id %in% tag_ids, ]
        detections_experiment = dim(vue_data)[1]
          # Total detections of experimental tags = 1669
      
      ## What was the percentage of tag detections from experiment tags
        detections_experiment.percent = detections_experiment / detections_all * 100
          # 98.23426 %
      
      ## What was the number of logged detections from other tags
        detections_other = detections_all - detections_experiment
          # Total detections from other tags =  30
        detections_other.percent = detections_other / detections_all * 100
          # 1.765745
      
    #### Hourly Detection Data Analysis ----
    ### Creating a dataframe (detection_df_4), to store results binned by hour
  
      ## Separating data by tag, receiver, and hour bin
      ## First assign each datum to an date and hour bin
        vue_data$hour_bin = floor_date(vue_data$datetime, unit = "hour")
      
      ## Then loop through combinations of date/hour, tag id, and receiver to get the number of detections from a tag at a receiver each hour
        detection_df_4 = data.frame()
        registerDoParallel(cores = n_cores)
        detection_df_4 = foreach(b = c(1:length(unique(vue_data$hour_bin))), .combine = rbind) %:% 
          foreach(r = c(1:length(unique(vue_data$receiver))), .combine = rbind) %dopar%{
            filtered = filter(vue_data, hour_bin == unique(vue_data$hour_bin)[b], 
                              receiver == unique(vue_data$receiver)[r])
            write_line = cbind(as.character(unique(vue_data$hour_bin))[b],
                               as.character(unique(vue_data$receiver))[r], 
                               dim(filtered)[1])
            return(write_line)
          }
        detection_df_4 = as.data.frame(detection_df_4)
        colnames(detection_df_4) = c('hour_bin', 'receiver', 'detections')
        detection_df_4$detections = as.numeric(as.character(detection_df_4$detections))
      
      ## Assigning distance condition based on receiver sn
        detection_df_4$distance = NA
        detection_df_4$distance[detection_df_4$receiver == receiver_60m] = '60 m'
        detection_df_4$distance[detection_df_4$receiver == receiver_508m] = '512 m'
      
    #### Statistical Analysis - Comparing Hourly Detections ----
      
      ### Hourly Receiver Comparison
      ## Shapiro Wilk test that distribution of each case are non-parametric (p-value)
        shapiro.test(detection_df_4$detections[detection_df_4$distance == '60 m']) 
          # p = 0.00151, p < 0.05 reject null hypothesis that distribution is normal
          # Conclusion: Distribution is not normal
        shapiro.test(detection_df_4$detections[detection_df_4$distance == '512 m']) 
          # p < 4.738e-09, p < 0.05 reject null hypothesis that distribution is normal
          # Conclusion: Distribution is not normal  
  
      ## Wilcox test - Selected to compare two groups (60m and 500m) non-parametrically
        wilcox.test(x = detection_df_4$detections[detection_df_4$distance == '60 m'], 
                    y = detection_df_4$detections[detection_df_4$distance ==  '512 m'],
                    paired = TRUE)
          # p-value < 2.2e-16 reject null hypothesis that difference between groups is zero
          # Conclusion: The two receivers logged signifcantly different detections from one another
      
        ## Ratio for Mean hourly detection rates for 60m receiver and 512 m receiver 
          mean_detections_60m  = mean(detection_df_4$detections[detection_df_4$distance == '60 m'])
            # 4.879699
          mean_detections_508m = mean(detection_df_4$detections[detection_df_4$distance == '512 m'])
            # 7.669173
          mean_detections_508m/mean_detections_60m
            # receiver at 508m has 1.571649 x as many detections as receiver in 60 m
      
      ## Number of hours where deeper receiver detection rates exceeds shallow receiver 
        total_hour_bins = length(unique(detection_df_4$hour_bin))
        times_508m_exceeds_60m = 0
        times_508m_equal_to_60m = 0
        for(hour in unique(detection_df_4$hour_bin)){
          subset_detection_df_4 = filter(detection_df_4, hour_bin == hour)
          if(subset_detection_df_4$detections[subset_detection_df_4$distance == '512 m'] >
             subset_detection_df_4$detections[subset_detection_df_4$distance == '60 m']){
            times_508m_exceeds_60m = times_508m_exceeds_60m + 1
          }else if(subset_detection_df_4$detections[subset_detection_df_4$distance == '512 m'] ==
                   subset_detection_df_4$detections[subset_detection_df_4$distance == '60 m']){
            times_508m_equal_to_60m = times_508m_equal_to_60m + 1
          }
        }
        times_508m_exceeds_60m
          # 107
        times_508m_equal_to_60m
          # 9
      
    #### Analysis Cleanup -----
      ### Saving workspace in results folder
        setwd(results_dir)
        save.image(file = 'Depth and Distance Validation Test Results.R')
        
  ###### Experiment 5: Multipath Confirmation - Tank Experiment ----
        
      ## Analysis Ideas
        # Break up analysis into three groups.
        # 1.	Control
        # 2.	Multipath (CPDI Not Expexted)
        # 3.	Multipath (CPDI Expected)
        
      ## Anna's helpful analysis ideas:
          # 1. Create table Event (single event would be do, 2a, 150m, 15m)
          # 2. Catagorize each event (control, multipath-no CPDI, multipath - CPDI)
          # 3. Code by type (CPDI predicted vs. not predicted)
          # 4. Detected - number of actual detections at receiver
          # Flip to wide table
            # columns - 
              # 1. Event
              # 2. Type
              # 3. Detected
          # 6. Fit GLM (Detected ~ type + event, family = binomial)
        
        
    #### Loading in and Cleaning Data Files ----
        setwd(data_dir)
      ### Importing VR100 data 
        vr100_raw = read.csv(file.path(data_dir, 'VR100_10320_D2016.07.15T20.49.23.csv'))
      ### Importing VR2 data from VUE
        vr2 = load_vemco(file.path(data_dir, 'VUE_Export_CPDI_Tank_Test.csv'))
      ### Importing meta data regarding when transmissions were sent recorded from Brendan's laptop durring experiment
        condStoreDM = read.csv(paste(data_dir, 'condStoreDM.csv', sep = ""), header = FALSE)
        condStoreDO = read.csv(paste(data_dir, 'condStoreDO.csv', sep = ""), header = FALSE)
        
      ### Creating vr100 datafame for manipulation
        vr100 = vr100_raw
        vr100$datetime = as.POSIXct(paste(vr100$Date, vr100$Time, sep = " "), format = "%Y-%m-%d %H:%M:%S", tz = "GMT")
        attr(vr100$datetime, "tzone") <- "HST"
        colnames(vr100) = c("date", "time", "chan", "freq", "type", "code_space", "tag_id", "serial_number", "data1", "units1", "slope1", "init1", "data2", "units2", "slope2", "int2", "signal_db", "gain_mode", "gain_db", "lat", "lon", "comment", "datetime")
        
      # creating vr2 dataframe for manipulation
        vr2$datetime = as.POSIXct(vr2$datetime, tz = "GMT")
        attr(vr2$datetime, "tzone") <- "GMT"
        vr2 = vr2[vr2$datetime >= as.POSIXct("2016-07-12 16:30:00"), ]
        
      ### Creating dataframes for direct and multipath conditions
        colnames(condStoreDO) = c('pluse_train_time', 'tag_range', 'tag_depth', 'relative_arrival_time_0','case_number', 'sub_case_number')
        colnames(condStoreDM)[1:3] = c('pluse_train_time', 'tag_range', 'tag_depth')
        colnames(condStoreDM)[4:23] = paste('relative_arrival_time_', 0:19, sep = "")
        colnames(condStoreDM)[24:25] = c('case_number', 'sub_case_number')
        condStoreDO$dm = 'direct only' # designates whether multipath or not
        condStoreDM$dm = 'multipath' # designates whether multipath or not
        
      ### Predicting whether multipaths produced CPDI effects
        ## DO conditions only have a direct signal, all should be detected by receivers. thus they are marked as control
          condStoreDO$cpdi_type = 'Control'
        ## To determine if multipath signals cause CPDI, we compare the relative arrival times of the multipath signals to the blanking interval (0.260 seconds). 
          condStoreDM$cpdi_type = 'multipath_no_cpdi_predicted'
          for(i in 1:dim(condStoreDM)[1]){
            if(any(condStoreDM[i,4:23] != 0 & condStoreDM[i,4:23] > 0.260)){
              condStoreDM$cpdi_type[i] = 'multipath_cpdi_predicted'
            }
          }
        
      ### Combining direct and multipath condition dataframes into a single dataframe
        ## condStoreDO to do this, condStoreDO needs to be padded out to be same size as DM. This meens adding 20 additional columns
          condStoreDO = cbind(condStoreDO[ ,1:4], matrix(data = NA, nrow = dim(condStoreDO)[1], ncol = 19), condStoreDO[ ,5:8])
          colnames(condStoreDO)[5:23] = paste('relative_arrival_time_', 1:19, sep = "")
        ## Combining dataframes
          condition = data.frame(rbind(condStoreDO, condStoreDM)) # Combining into condition dataframe
        ## Converting matlab date nubmers to POSIX datetime objects    
          condition$datetime = datenum2posix(condition$pluse_train_time) 
        
      ### Ordering dataframes by increasing time of detection
        condition = condition[order(condition$datetime), ] # sorting by time
        vr2 = vr2[sort.list(x = vr2$datetime, decreasing = FALSE), ]
        vr100 = vr100[sort.list(x = vr100$datetime, decreasing = FALSE), ]
        rownames(vr100) = NULL
        rownames(vr2) = NULL
        
      ### Converting vr2$receiver and vr2$tag_id to numeric
        vr2$receiver = as.numeric(levels(vr2$receiver)[vr2$receiver])
        vr2$tag_id = as.numeric(levels(vr2$tag_id)[vr2$tag_id])
        
      ### Pulling in all syncs to see if there is major drift in receivers/vr100 clocks
        vr100_syncs = vr100[which(vr100$tag_id == 36817), ]
        vr2_syncs = vr2[which(vr2$tag_id == 36817), ]
        sync_data = data.frame(rbind(cbind(as.character(vr100_syncs$datetime), vr100_syncs$tag_id, 'vr100', vr100_syncs$signal_db),
                                     cbind(as.character(vr2_syncs$datetime), vr2_syncs$tag_id, vr2_syncs$receiver, NA)))
        colnames(sync_data) = c("datetime", "tag_id", "receiver", "signal_db")
        sync_data = sync_data[order(as.POSIXct(sync_data$datetime)), ]
        
      ### Removing tag detections prior to experiment day
        remove_prior_dates = as.POSIXct("2016-07-12 10:31:0 HST")
        vr2 = vr2[vr2$datetime >= remove_prior_dates, ]
        vr100 = vr100[vr100$datetime >= remove_prior_dates, ]
        
      ### Checking sync tag detected same number of times
        length(which(vr100$tag_id == 36817)) # 1
        length(which(vr2$tag_id == 36817))   # 1 
        
      ### Checking which receiver were part of experiment
        unique(vr2$receiver)
        # 129577
        
      ### Checking which tags were logged on receiver
        unique(vr2$tag_id)
        # 36814 36817 18266 47515
        unique(vr100$tag_id)
        # 18266 47515 36814 36817
        
      ### Separating two receivers
   
        ## For VR100
        vr100 = vr100[max(which(vr100$tag_id == 36817)):dim(vr100)[1], ]
        
      ### Adjusting vr100 dates to sync with vr2 dates. Used test tag on VR100 and VR2 in the lab for one detection
        # VR2 detected tag at "2016-03-16 00:14:16 GMT"
        # VR100 detected tag at"2016-03-16 00:13:40 GMT"
        
        ## Time offset  between VR2 and VR100
          vr2_vr100_time_diff = difftime(time1 = max(vr100$datetime[vr100$tag_id == 36817]),
                                       time2 = min(vr2$datetime[vr2$tag_id == 36817]))
              # Time difference of -1 secs
          vr100$datetime = vr100$datetime - vr2_vr100_time_diff
          
          all_data = data.frame(matrix(data = 0, nrow = dim(vr2)[1] + dim(vr100)[1], ncol = 4))
          colnames(all_data) = c("datetime", "tag_id", "receiver", "signal_db")
          all_data$datetime = c(vr2$datetime, vr100$datetime)
          all_data$tag_id = c(vr2$tag_id, vr100$tag_id)
          all_data$receiver = c(vr2$receiver, rep('vr100', dim(vr100)[1]))
          all_data$signal_db = c(rep(NA, length(vr2$receiver)), vr100$signal_db)
        
        ## Ordering by time
          all_data = all_data[order(all_data$datetime), ]
          rownames(all_data) = NULL
          all_data$diff_time = 0
          for(i in 2:length(all_data$diff_time)){
            all_data$diff_time[i] = difftime(all_data$datetime[i], all_data$datetime[i-1], units = 'secs')
          }
          all_data$exp = NA
          
        ## Determining where experiments begin and end
          all_data$exp[all_data$tag_id == 36817] = 'Sync'
        
          all_data = all_data[all_data$datetime >= min(all_data$datetime[all_data$tag_id == 18266]), ]
        
        ## How many detections did each receiver log?
          print(paste('Receiver 129577 detected', dim(all_data[which(all_data$receiver == 129577 & all_data$tag_id != 18266), ])[1], 'transmissions', sep = " "))
            # Receiver 129577 detected 594 transmissions"
          print(paste('The VR100 detected', dim(all_data[all_data$receiver == 'vr100', ])[1], 'transmissions', sep = " "))
            # The VR100 detected 929 transmissions
          
        
        ### Matching conditions to recorded detection data (all_data)
          rownames(all_data) = NULL
          all_data$tag_range = NA
          all_data$tag_depth = NA
          all_data$case_number = NA
          all_data$sub_case_number = NA
          all_data$dm = NA
          all_data$cpdi_type = NA
          all_data$datetime = all_data$datetime - 6 # Adjusted this by comparing when the first test started ~18:58 on all data to first transmission sent in condition
          
          for(i in 1:length(all_data$datetime)){
            for(r in 1:length(condition$datetime)){
              # If a detection occurred within an 8 second window
              if(all_data$datetime[i] - 4 <= condition$datetime[r] & all_data$datetime[i] + 4 >= condition$datetime[r]){ 
                all_data$tag_range[i] = condition$tag_range[r]
                all_data$tag_depth[i] = condition$tag_depth[r]
                all_data$case_number[i] = condition$case_number[r]
                all_data$sub_case_number[i] = condition$sub_case_number[r]
                all_data$dm[i] = condition$dm[r]
                all_data$cpdi_type[i] = condition$cpdi_type[r]
              }
            }
          }
        
        ## Changing all_data$sub_case_number from an integer to a character to match other formats
          all_data$sub_case_number = as.character(all_data$sub_case_number)
          all_data$sub_case_number[all_data$sub_case_number == '1'] = 'a'
          all_data$sub_case_number[all_data$sub_case_number == '2'] = 'b'
          all_data$sub_case_number[all_data$sub_case_number == '3'] = 'c'
          all_data$sub_case_number[all_data$sub_case_number == '0'] = ''
          all_data$exp = paste(all_data$dm, all_data$case_number, all_data$sub_case_number, sep = "")
          
        ### Creating a new dataframe for noting if each transmission was detected and by what (VR2 or vr100)
          exp_df = as.data.frame(matrix(0, nrow = length(condition$datetime), ncol = 3)) 
          colnames(exp_df) = c('datetime', 'v129577', 'vr100')
          exp_df$datetime = condition$datetime
          for(i in 1:length(all_data$datetime)){
            for(r in 1:length(exp_df$datetime)){
              if(all_data$datetime[i]-2 <= exp_df$datetime[r] & all_data$datetime[i]+2 >= exp_df$datetime[r]){
                if(all_data$receiver[i] == 129577){
                  exp_df$v129577[r] = 1
                }else if(all_data$receiver[i] == 'vr100'){
                  exp_df$vr100[r] = 1
                }
              }
            }
          }
        
        detection_df_5 = melt(exp_df, id = c('datetime'))
        model_data = merge(detection_df_5, condition[ ,c('tag_range', 'tag_depth', 'case_number', 'sub_case_number', 'dm', 'cpdi_type', 'datetime')], by = intersect(names(exp_df),names(condition[ ,c('tag_range', 'tag_depth', 'case_number', 'sub_case_number', 'dm', 'cpdi_type', 'datetime')])))
        colnames(model_data)[c(2, 3)] = c('receiver', 'detected')
        model_data$sub_case_number[model_data$sub_case_number == 0] = 1
        model_data$sub_case_number = letters[model_data$sub_case_number]
        model_data = model_data[model_data$receiver != 'vr100', ] # removing VR100 detections since we're not evaluation efficiency of VR100. 
       
        ## Need to adjust experiment 4 as its currently one case for both receiver conditions.
          model_data$sub_case_number[model_data$case_number == 4 & model_data$tag_range == 508] = 'b'
   
      #### Constructing a GLM Model for detection
        ### Data setup
        # 1. Create table Event (single event would be do, 2a, 150m, 15m)
          event_table = model_data
  
        # 2. Determining detection for predictions based on prediction type (1 if detection is predicted, 0 if not)
          event_table$event = paste(event_table$case_number, event_table$sub_case_number, sep = "")
          event_table$predicted_outcome = NA
          event_table$predicted_outcome[event_table$cpdi_type == 'multipath_cpdi_predicted' ] = 0
          event_table$predicted_outcome[event_table$cpdi_type == 'Control' ] = 1
          event_table$predicted_outcome[event_table$cpdi_type == 'multipath_no_cpdi_predicted' ] = 1
          
        # 3. Renaming variables for model fit
          colnames(event_table)[which(colnames(event_table) == 'receiver')] = 'predicted_observed'
          colnames(event_table)[which(colnames(event_table) == 'cpdi_type')] = 'cpdi_prediction'
          event_table$cpdi_prediction = as.factor(event_table$cpdi_prediction)
          colnames(event_table)[which(colnames(event_table) == 'event')] = 'experiment_analgoue'
          event_table$experiment_analgoue = as.factor(event_table$experiment_analgoue)
          
        ### Fitting GLMs to data
          ## Glm fit models whether transmission was detected as a function of our model's prediction
            cpdi_glm.interaction = glm(detected ~  cpdi_prediction + experiment_analgoue + cpdi_prediction*experiment_analgoue, data = event_table, family = binomial)
            summary(cpdi_glm.interaction)
            table_3a = summary(cpdi_glm.interaction)$coefficients
            # write.csv(table_3a, file = file.path(results_dir, 'table_3a.csv'))
            
          ## removing non-significant interaction terms
            cpdi_glm.no_interaction = glm(detected ~ cpdi_prediction + experiment_analgoue , data = event_table, family = binomial)
            summary(cpdi_glm.no_interaction)
            table_3b = summary(cpdi_glm.no_interaction)$coefficients
            write.csv(table_3b, file = file.path(results_dir, 'table_3.csv'))
            
          ### How many times did observations deviate from predictions?
            ## With control transmissions
              sum(abs(event_table$predicted_outcome - event_table$detected))
              length(event_table$predicted_outcome)
          
            ## Without control transmissions
            sum(abs(event_table$predicted_outcome[event_table$cpdi_prediction != "Control"] - event_table$detected[event_table$cpdi_prediction != "Control"]))
            length(event_table$predicted_outcome[event_table$cpdi_prediction != "Control"])
  
        #### Model Diagnostics performed with Anna on 26 May 2017
          library(MuMIn)
          ## Dredging model to determine most parsimonious model
          options(na.action=na.fail)
          dredge(cpdi_glm.no_interaction, extra="R^2")
          
          resp<-event_table[, "detected"]
          pred<-event_table[, c("cpdi_prediction", "experiment_analgoue")]
          ## Looking at model residuals
          plot(cpdi_glm.no_interaction)
          
          ## Comparing model to data
          library(car)
          marginalModelPlot(cpdi_glm.no_interaction)
          resp<-event_table[event_table$experiment_analgoue != '3b', "detected"]
          pred<-event_table[event_table$experiment_analgoue != '3b', c("cpdi_prediction", "experiment_analgoue")]
          
          ## Which model term explains how much variation?
          library(hier.part)
          hier.part(resp, pred, family=binomial)
          
          ## Tukey Contrasts between experimental analogue groups
          library('multcomp')
          tukey_comparison = glht(cpdi_glm.no_interaction, mcp(experiment_analgoue="Tukey"))
          plot(tukey_comparison)
          class(event_table$experiment_analgoue)
          event_table$experiment_analgoue = as.factor(event_table$experiment_analgoue)
          event_table$cpdi_prediction = as.factor(event_table$cpdi_prediction)
  
          ## Visualizing differences in glm factors 
          library(visreg)
          visreg(cpdi_glm.no_interaction)
          dredge(cpdi_glm.interaction)
          anova(cpdi_glm.no_interaction)
          dredge(cpdi_glm.interaction)
          .8*64
          
        ## Comparing glm with interaction to glm without interaction term
          anova(cpdi_glm.no_interaction, cpdi_glm.interaction, test="Chisq")
          
        ## Comparing best glm by stepwise addition of terms
          anova(cpdi_glm.no_interaction, test="Chisq")
          
        ## What happens if we pull out control cases that we know failed?
          cpdi_glm.no_interaction = glm(detected ~ cpdi_prediction + experiment_analgoue, data = subset(event_table, experiment_analgoue !="3a"), family = binomial)
          dredge(cpdi_glm.no_interaction) # Best model still includes both terms
          plot(cpdi_glm.no_interaction) # Checking out residuals
         # Is there equal variance across experimental analogues?
          plot(event_table$experiment_analgoue[event_table$experiment_analgoue != "3a"], residuals(cpdi_glm.no_interaction))
        
        ## What happens if we only remove 2b 
          cpdi_glm.no_interaction = glm(detected ~ cpdi_prediction + experiment_analgoue, data = subset(event_table, experiment_analgoue !="2b"), family = binomial)
          dredge(cpdi_glm.no_interaction)
          summary(cpdi_glm.no_interaction) # Best model still includes both terms
          visreg(cpdi_glm.no_interaction, scale="response") # Visualizing the differences in GLM factors
  
        #### Analysis Cleanup
          setwd(results_dir)
          save.image(file = 'HIMB Multipath Tank Test Results.R')
          
  ###### Miscelanious Plots and Figures
          
        label_size = 2.5
      #### Figure 2. Direct and First Multipath arrival comparision
        ### Setting up data repository and determining distance to calculate over
        direct_arrival_times = c()
        first_multipath_arrival_times = c()
        calc_distance = 2000
        simulated_sound_speed = 1530
        depth = 100
        
        ### Running arrival simulation
        for(dist in 1:calc_distance){
          arrival_times = get_surf_paths(surface_to_tag_distance = depth, 
                                         surface_to_receiver_distance = depth, 
                                         tag_to_receiver_horizontal_distance = dist, 
                                         bottom_depth = 500, 
                                         speed_of_sound = simulated_sound_speed,
                                         ave_max_detection_radius = 1000)
          direct_arrival_times = c(direct_arrival_times, arrival_times[1])
          first_multipath_arrival_times = c(first_multipath_arrival_times, arrival_times[2])
        }
        
        ### Plot Set-up
        setwd(figure_dir)
        png('Figure 2 - Direct and multipath arrival times.png', width = 2000, height = 1000)
        par(mar = mar.default + c(4, 4, 0, 0), oma = c(0, 1, 0 ,0)) 
        ### Plotting
        plot(y = direct_arrival_times, x = 1:calc_distance, 
             type = 'l', 
             cex = label_size,
             cex.lab = label_size,
             cex.axis = label_size,
             xlim = c(0, 1000), ylim = c(0, .65),
             xaxs = "i", yaxs = "i",
             # main = 'Arrival Time of Direct and First Surface Reflected Multipath \n
             #     Receiver Depth 250 m, Transmitter Depth = 250 m',
             ylab = 'Arrival Time (s)', xlab = 'Distance Between Tag and Receiver (m)', font = 1,  cex.lab = label_size, cex.axis = label_size, lwd = 5, col = 'black')
        # lines(y = first_multipath_arrival_times, x = 1:calc_distance, type = 'l', lty = 4)
        text(x = 45, y = 0.05, srt = 23, labels = "Direct Path", cex = label_size)
        
        for(depth in seq(50, 550, 100)){
        ### Running arrival simulation
          first_multipath_arrival_times = c()
          direct_arrival_times = c()
          
        for(dist in 1:calc_distance){
          arrival_times = get_surf_paths(surface_to_tag_distance = depth, 
                                         surface_to_receiver_distance = depth, 
                                         tag_to_receiver_horizontal_distance = dist, 
                                         bottom_depth = 500, 
                                         speed_of_sound = simulated_sound_speed,
                                         ave_max_detection_radius = 1000)
          direct_arrival_times = c(direct_arrival_times, arrival_times[1])
          first_multipath_arrival_times = c(first_multipath_arrival_times, arrival_times[2])
        }
          mp_index = c(!is.na(first_multipath_arrival_times), FALSE)
        lines(y = first_multipath_arrival_times[mp_index], x = c(1:calc_distance)[mp_index], type = 'l', lty = 3, lwd = 5)
        text(x = 30, y = first_multipath_arrival_times[40]+.024, labels = paste(depth, 'm'), cex = label_size)
        interference_index = which(first_multipath_arrival_times[mp_index] - direct_arrival_times[mp_index] <= 0.26)
        if(length(interference_index) > 0){
            lines(y = first_multipath_arrival_times[interference_index], x = c(1:calc_distance)[interference_index], lty = 1, lwd = 5)
          }
        }
        dev.off()
        
        
    #### Figure 8. Comparing Deep and Shallow Water Detection Profiles
      ### Plot Set-up
        mar.default <- c(5,8,4,2) + 0.1
        
        label_size = 2.75
        png('Figure 5 - Comparing deep and shallow detection profiles.png', width = 2550, height = 3300)
        par(mfrow = c(2, 1), mar = mar.default, mgp = c(40, 1, 0), oma = c(7,7,7,7)) 
  
        ### Plotting Deep Water Ranging Experiment Results
        exp_1_data_to_plot = detection_df_1[which(detection_df_1$height_of_tag_off_bottom == 1 & detection_df_1$height_of_receiver_off_bottom == 1 & detection_df_1$direction == "In" & detection_df_1$day_night == "Day"), ]
        plot((candidate_model_predictions_1$model_119$predicted_rates) ~ c(0:1000), 
             type = 'l', 
             xlim = c(0, 1000), xlab = "",
             ylim = c(0, 60), ylab = "",
             cex = 4,
             cex.lab = label_size,
             cex.axis = label_size,
             lwd = 10)
        points((exp_1_data_to_plot$detections) ~ exp_1_data_to_plot$distance, pch = 19, cex = label_size)
        
        ## Adding Axis labels
        mtext("Distance Between Tag and Receiver (m)", side=1, line=6, cex = 4)
        mtext("# of Transmissions Detected", side=2, line=6, cex = 4)
        
        
        ## Adding AMDR to plot
        points(x = candidate_model_summary_1$ave_max_distance[1], y = 3, pch = '*', cex = 10)
        text(x = candidate_model_summary_1$ave_max_distance[1] +50, y = 4, labels = paste('AMDR \n', candidate_model_summary_1$ave_max_distance[1]-1, 'm'), col = 'black', cex = label_size)
        
        ## Adding CPDI to plot
        if(candidate_model_summary_1$cpdi_extent[3] != 0){
          points(x = fivenum(candidate_model_summary_1$cpdi_extent)[3], y = max(candidate_model_predictions_1$model_119$predicted_rates), pch = '*', cex = 10)
          text(x = fivenum(candidate_model_summary_1$cpdi_extent)[3] + 15, y = 30, labels = paste('CPDI Extent \n ', candidate_model_summary_1$cpdi_extent[3], 'm'), cex = label_size)
        }
        
        ## Adding A Text
        text(x = 1000, y = 55, labels = "A", cex = 6)
        
        ### Plotting Shallow Water Ranging Experiment Results
        exp_2_data_to_plot = detection_df_2[which(detection_df_2$height_of_receiver_off_bottom == 1 & detection_df_2$direction == "in" & detection_df_2$day_night == "day"), ]
        plot((candidate_model_predictions_2$model_123$predicted_rates) ~ c(0:1200), 
             type = 'l', 
             xlim = c(0, 1200), xlab = "",
             ylim = c(0, 60), ylab = "",
             cex = 4,
             cex.lab = label_size,
             cex.axis = label_size,
             lwd = 10)
        points((exp_2_data_to_plot$detections) ~ exp_2_data_to_plot$distance, pch = 19, cex = label_size)
        
        
        mtext("Distance Between Tag and Receiver (m)", side=1, line=6, cex = 4)
        mtext("# Transmissions Detected", side=2, line=6, cex = 4)
        
        ## Adding AMDR to plot
        points(x = levels(candidate_model_summary_2$ave_max_distance[1])[candidate_model_summary_2$ave_max_distance[1]], y = 3, pch = '*', cex = 10)
        text(x = as.numeric(levels(candidate_model_summary_2$ave_max_distance[1])[candidate_model_summary_2$ave_max_distance[1]]) + 75, y = 4, labels = paste('AMDR \n', as.numeric(levels(candidate_model_summary_2$ave_max_distance[1])[candidate_model_summary_2$ave_max_distance[1]])-1, 'm'), col = 'black', cex = label_size)
        
        ## Adding B Text
        text(x = 1200, y = 55, labels = "B", cex = 6)
        
        dev.off()
        
        
      #### Figure 9. Comparing detections and pings
        ### Converting data classes to numeric for ggplot
          metadata_df$pings = as.number(metadata_df$pings)
          metadata_df$detections = as.number(metadata_df$detections)
        ### Creating a dataframe just for plotting
          meta_for_plot = as.data.frame(rbind(cbind("rec_50m", 'pings / 8', mean(metadata_df$pings[metadata_df$receiver == "receiver_50m"]) / 8), cbind("rec_50m", 'detections', mean(metadata_df$detections[metadata_df$receiver == "receiver_50m"])),
                                              cbind("rec_212m", 'pings / 8', mean(metadata_df$pings[metadata_df$receiver == "receiver_212m"]) / 8), cbind("rec_212m", 'detections', mean(metadata_df$detections[metadata_df$receiver == "receiver_212m"]))))
          colnames(meta_for_plot) = c('receiver', 'Metric', 'Detected')
          meta_for_plot$Detected = as.number(meta_for_plot$Detected)
        ### Plotting Figure
          setwd(figure_dir)
          png('exp_3_pings_vs_detections.png', height = 1800, width = 1800)
          (ggplot(meta_for_plot, aes(factor(receiver), Detected, fill = Metric)) +
              geom_bar(stat = "identity", position = "dodge") +
              scale_fill_grey(start = 0, end = .7) + theme_bw())
          dev.off()
          setwd(results_dir)
          
          
          #### Figure X. Comparing expected number of detections and number of detections observed
          plot(x = expected_detections$n_tags, y = expected_detections$hourly_detections_total, 
               pch = 19, 
               main = "Observed and Predicted Detections for deployed tags",
               xlim = c(0, 18), xlab = "# of Tags",
               ylim = c(0, 400), ylab = "# of Detections")
          sum_of_hourly_detections_by_ntags = aggregate(detection_df_1$detections, by = list(detection_df_1$hour_bin, detection_df_1$tags_per_hour), FUN = sum)
          colnames(sum_of_hourly_detections_by_ntags) = c("hour_bin", "n_tags", "n_detections")
          points(x = sum_of_hourly_detections_by_ntags$n_tags, y = sum_of_hourly_detections_by_ntags$n_detections, pch = 19, col = 'red')
          hourly_detections_by_distance = aggregate(detection_df_1$detections, by = list(detection_df_1$hour_bin, detection_df_1$tags_per_hour, detection_df_1$distance), FUN = sum)
          colnames(hourly_detections_by_distance) = c("hour_bin", "n_tags", "distance", "n_detections")
          plot(x = hourly_detections_by_distance$distance, y = hourly_detections_by_distance$n_detections)
  
          
  ##### Figure 4. Plotting a chart where all experiments took place
    #### Commented Out because loading in bathy file takes forever and a day. Uncomment if ever trying to produce these plots again
  # ### Loading in Data files for plotting
  #   ## Locations of each Experiment
  #     exp_data = read.csv('/Users/stephenscherrer/Desktop/Experiment Coordinates for Figure.csv', stringsAsFactors = FALSE)
  #   ## Importing Oahu 5m Bathy file (.xyz format exported from )
  #     oahu_5m_bathy = read.bathy('/Users/stephenscherrer/Google Drive/Weng Lab/Data/Oahu 5m Bathymetry/oahu_5m_bathy_lat_lon.xyz')
  #   ## Subsetting High Res Bathymetry File
  #     oahu_subset = subsetBathy(oahu_5m_bathy, x = c(-157.81,-157.9), y = c(21.24, 21.29), locator = FALSE)
  #   ## Loading in Low Res Bathyemtry for background plot
  #     oahu_1km_bathy = getNOAA.bathy(lon1 = -158.4, lon2 = -157.5, lat1 = 21.8, lat2 = 21.1, resolution = 1)
  #     
  #   ### Assigning Plot Parameters
  #     ## Creating color palettes for Bathymetry
  #     blues <- c("lightsteelblue4", "lightsteelblue3",
  #                "lightsteelblue2", "lightsteelblue1")
  #     greys <- c(grey(0.6), grey(0.93), grey(0.99))
  #     
  #     ## Plot Parameters for experiment locations
  #     plot_type = c(17, 19)
  #     plot_type_outline = c(2, 1)
  #     plot_col = c('yellow', 'red', 'green', 'purple')
  #   
  # #### Plotting Bathymetry Maps  
  # ### Low Res Bathymetry background
  #   ## Setting up plot parameters
  #   png(file.path(figure_dir, "oahu_1km_bathy.png"), height = 8, width = 8.5, units = "in", res = 224)
  #   ## Plotting bathymetric basemap
  #     plot.bathy(oahu_1km_bathy, image = TRUE, land = TRUE, lwd = 0.5, bpal = list(c(0, max(oahu_1km_bathy), 'white'), c(min(oahu_1km_bathy), 0, blues)), deep=min(oahu_1km_bathy), shallow=1, step=100,  drawlabel=TRUE)
  #     ## Making the coastline more visible
  #     plot(oahu_1km_bathy, deep = 0, shallow = 0, step = 0, lwd = 0.7, add = TRUE)
  #     ## Plotting points where experiment occurred
  #     points(exp_data$Lat ~ exp_data$Lon, pch = plot_type[as.factor(exp_data$Type..Tag.Receiver.)], col = plot_col[as.factor(exp_data$Experiment)])
  #     ## outlining those points
  #     points(exp_data$Lat ~ exp_data$Lon, pch = plot_type_outline[as.factor(exp_data$Type..Tag.Receiver.)], col = 'black')
  #     ## Adding scale bars
  #     scaleBathy(oahu_1km_bathy, deg = 0.1)
  #     dev.off()
  #     
  # ### High Res Inset Chart
  #     ## Setting up plot parameters
  #     png(file.path(figure_dir, "5m_inset.png"), height = 2.7253, width = 2.8, units = "in", res = 1000)
  #     ## Plotting bathymetric basemap
  #     plot.bathy(oahu_subset, image = TRUE, land = TRUE, lwd = 0.5, bpal = list(c(0, max(oahu_subset, na.rm = TRUE), 'white'), c(min(oahu_subset, na.rm = TRUE), 0, blues)), deep=min(oahu_subset, na.rm = TRUE), shallow=max(oahu_subset, na.rm = TRUE), step=100,  drawlabel=TRUE, cex = 2)
  #     ## Making the coastline more visible
  #     plot(oahu_subset, deep = 0, shallow = 0, step = 0, lwd = 0.7, add = TRUE)
  #     ## Plotting points where experiment occurred
  #     points(exp_data$Lat ~ exp_data$Lon, pch = plot_type[as.factor(exp_data$Type..Tag.Receiver.)], col = plot_col[as.factor(exp_data$Experiment)])
  #     ## outlining those points
  #     points(exp_data$Lat ~ exp_data$Lon, pch = plot_type_outline[as.factor(exp_data$Type..Tag.Receiver.)], col = 'black')
  #     ## Adding scale bars
  #     scaleBathy(oahu_1km_bathy, deg = 0.01)
  #     dev.off()
   
  ###### Script Cleanup ----
  save.image('completed_script.R')
      # load(file.path(results_dir, 'completed_script.R'))
  run_time = proc.time() - script_timer
  send_push(user = 'uGEHvA4hr37tsrCCtpSv4sUUxVuTqN', message = paste("CPDI Analysis Completed in", round(run_time[3]), 'seconds!'))
