#### Mechanistic Model for Predicting the Occurrence of CPDI

#### Written by Stephen Scherrer
#### Written 24 April 2015 
#### Last Modified: 22 June 2017
#### Email: scherrer@hawaii.edu
#### All rights preserved, all wrongs traversed

#### Predicting the occurance of interfering bounce paths for Vemco
#### Acoustic hardware over a range of vertical and horizontal tag
#### and receiver depth distributions producing CPDI (the doughnut effect). 

#### See Example at end of script. 

#### This script includes two user functions.
 ### predict_cpdi_interference()
  ##    This function produces a matrix where rows corrospond to transmitter depth and columns corrospond to transmitter distance from receiver.
  ##      Use of the plot argument will produce an image for easy viewing
  ##    inputs: 
  ##            bottom_depth - Depth in meters of study site. 
  ##            ave_max_detection_radius - Distance at which tags are no longer detected by receivers. Determined by range testing.
  ##            receiver_depth - Depth of receiver  in meters relative to the surface. (ie: If receiver is on rope 3 m above seafloor, receiver depth = bottom depth - 3)
  ##            speed_of_sound - Speed of sound in enviornment. Assummed to be 1530 m/sec unless specified otherwise
  ##            max_horizontal_dist - Maximum distance model should search over
  ##            evaluation_interval - To reduce computational speed, model can bin search area. By default, model makes a prediction for every 1 meter.  
  ##            blanking_interval - The receiver's blanking interval in seconds (default is 0.260)
  ##            plot - Whether or not to produce a plot showing where detections are predicted to occur. Light grey indicates positions where a tag can be heard while dark grey are positions where tag is not heard
  ##            ...  - Additonal plot options found in plot heat map function. Notably, save_file = TRUE will save a plot to the working directory
 
 ### rank_receiver_depths()
  ##    This function generates a list of the optimum heights for placing a receiver based on the number of total  positions that a tag can be detected
  ##    Inputs:
  ##            bottom_depth - Depth in meters of study site. 
  ##            max_horizontal_dist - Maximum distance model should search over
  ##            ave_max_detection_radius - Distance at which tags are no longer detected by receivers. Determined by range testing.
  ##            evaluation_interval - To reduce computational speed, model can bin search area. By default, model makes a prediction for every 1 meter.  
  ##            speed_of_sound - Speed of sound in enviornment. Assummed to be 1530 m/sec unless specified otherwise
  ##            blanking_interval - The receiver's blanking interval in seconds (default is 0.260)
  ##            tag_depth_range - An optional vector of depths a tag may appear of the format min_depth:max_depth. If a fish is known to be present at depths of 120-300 m, it makes no sense to evaluate positions shallower than 120 m.  


#### Model Functions ------------------------------------------------

## Calculating arrival time of higher order transmission paths
get_surf_paths = function(surface_to_tag_distance, 
                          surface_to_receiver_distance, 
                          tag_to_receiver_horizontal_distance = x, 
                          bottom_depth, 
                          speed_of_sound,
                          ave_max_detection_radius){
  
  R = 0 # Number of Reflections in bounce path
  
  ## Path can be thought of as the geometric summation of a series of similar
    ## triangles which can be unfolded into a large right triangle (Pascal's triangle) 
    ## where the base is the horizontal distance between the tag and receiver, and the 
    ## leg is the sum of the vertical distance between the tag and receiver, as well 
    ## as the depth of the water column for any bounce paths that are effected by
    ## substrate. 
  
  surf_paths = (1/speed_of_sound) * # (m/s) to convert distance to seconds
    sqrt( # Implementation of Pathagorean theorum
      tag_to_receiver_horizontal_distance ^ 2 + # Horizontal distance tag to receiver
                          (surface_to_tag_distance + # Vertical distance of tag 
                           #to acoustically reflective surface
                           (R-1) * # Term to constrain addition of bottom depth term
                                    # and the number of times it is implemented
                             bottom_depth + # Bottom depth term only included when
                                    # number of reflections in path 
                             abs(((R-1) %% 2 * bottom_depth) - # Path either
                                  # encounters receiver from below or above. When
                                  # number of reflections that begin at surface is
                                  # even, reflection at receiver is from below.
                                  # Modulo is 1 when R is odd, which makes this term
                                  # drop out. When even, this term is equivilent to
                                  # bottom depth. 
                                   surface_to_receiver_distance))^2) # Because path
                                  # encounter from below height is only equivilent 
                                  # height of receiver from bottom and not full water
                                  # column, subtracting the surface to receiver
                                  # distance gives the benthos to receiver distance.
                                  # Absolute value of the whole term is taken so 
                                  # when modulo evaluates 0, the added height to 
                                  # triangle leg is the surface to tag distance and
                                  # should be positive. 
  
  ## Until the length of higher order reflective paths is greater than the maximum
    ## path length (linearly transformed by speed of sound)
  while (surf_paths[length(surf_paths)] <= (ave_max_detection_radius * (1/speed_of_sound))){
    ## Increment order of reflections by one
    R = R + 1 
    ## Then re-run surf paths formula, and add to surf_paths output
    surf_paths = c(surf_paths, 
                  (1/speed_of_sound) * sqrt(tag_to_receiver_horizontal_distance^2 + 
                                 (surface_to_tag_distance + (R-1) * bottom_depth + 
                                    abs(((R-1) %% 2 * bottom_depth) - 
                                          surface_to_receiver_distance))^2))
  }
  ## When greatest surf path is larger than maximum surf path, stop running loop
    ## and remove impossibly long surf path from output
  surf_paths = surf_paths[-(length(surf_paths))]
  ## Then return
  return (surf_paths)
}

get_benthos_paths = function(surface_to_tag_distance = dt, 
                             surface_to_receiver_distance = dr,
                             tag_to_receiver_horizontal_distance = x, 
                             bottom_depth = db, 
                             speed_of_sound = v,
                             ave_max_detection_radius = amdr){
  R = 1
  ### See get_surf_paths() function for explination of geometry and equation below
  ### But reverse 
  benth_paths = (1/speed_of_sound) * sqrt(tag_to_receiver_horizontal_distance^2 + 
                               (bottom_depth-surface_to_tag_distance + 
                                  (R-1) * bottom_depth + 
                                  abs(((R) %% 2 *bottom_depth) - 
                                        surface_to_receiver_distance))^2)

  while (benth_paths[length(benth_paths)] <= (ave_max_detection_radius * (1/speed_of_sound))){
    R = R + 1
    benth_paths = c(benth_paths, 
                    (1/speed_of_sound) * sqrt(tag_to_receiver_horizontal_distance^2 + 
                               (bottom_depth-surface_to_tag_distance + (R-1) * bottom_depth + abs(((R) %% 2 *bottom_depth) - surface_to_receiver_distance))^2)) 
  }
  benth_paths = benth_paths[-(length(benth_paths))]
  return (benth_paths)
}

evaluate_interference = function(surface_to_tag_distance, 
                                 surface_to_receiver_distance, 
                                 tag_to_receiver_horizontal_distance, 
                                 bottom_depth, 
                                 speed_of_sound,
                                 ave_max_detection_radius,
                                 blanking_interval,
                                 plot = TRUE){
  
  ## Takes arrival times of surface path and benthic paths, and subtracts the time
    ## of the direct path. Then evaluates this against the blanking interval to 
    ## determine how many signal paths interfere with the blanking time for a
    ## choosen variable set. Returns number of conflicting signals. Will return 1 if 
    ## tag is out of range of the receiver
  
  arrival_times = unique(c(get_surf_paths(surface_to_tag_distance, 
                                   surface_to_receiver_distance, 
                                   tag_to_receiver_horizontal_distance, 
                                   bottom_depth, speed_of_sound, ave_max_detection_radius), 
                    get_benthos_paths(surface_to_tag_distance, 
                                      surface_to_receiver_distance, 
                                      tag_to_receiver_horizontal_distance, 
                                      bottom_depth, speed_of_sound, ave_max_detection_radius)))
  if(length(arrival_times) != 0){
    ## Getting arrival time of multipath signals relative to the direct path
    arrival_differences = arrival_times - min(arrival_times)
    ## Return the number of multipaths that exceed the length of the blanking interval
    return(sum((arrival_differences > blanking_interval)))
  }else{
    ## Return one if tag is out of range of the receiver
    return(1)
  }
}

plot_heat_map = function(matrix_to_plot, 
                         color_palette = sort(gray.colors(n = 2), decreasing = TRUE), 
                         compress_x_labels = 1, 
                         compress_y_labels = 1, 
                         save_file = FALSE, 
                         save_description = NULL, 
                         main = NULL){
  matrix_to_plot[matrix_to_plot > 0] = 1
  if(save_file == TRUE){
    if(is.null(save_description) == FALSE){
      png(save_description)
    }else{
      title = 'CPDI'
      png(paste(as.character(Sys.time()),'CPDI Simulation Results.png'))
    }
  }
  image(x = t(apply(matrix_to_plot, 2, rev)), 
        col = rev(color_palette), 
        xaxt = 'n', yaxt = 'n', 
        xlab = 'Distance Between Tag and Receiver (m)', 
        ylab = 'Depth (m)',
        main = main)
  
  ## scaling and positioning x axis labels 
  if(is.null(compress_x_labels) == TRUE){
    compress_x_labels = 1
  }
  xaxis_n = seq(1, ncol(matrix_to_plot), by = compress_x_labels)
  xaxis_labels = colnames(matrix_to_plot)[xaxis_n]
  xaxis_pos = xaxis_n / ncol(matrix_to_plot)
  axis(side = 1, labels = xaxis_labels, at = xaxis_pos)
  
  
  ## scaling and positioning y axis labels
  if(is.null(compress_y_labels) == TRUE){
    compress_y_labels = 1 
  }
  yaxis_n = seq(1, nrow(matrix_to_plot), by = compress_y_labels)
  yaxis_labels = rev(rownames(matrix_to_plot)[yaxis_n])
  yaxis_pos = yaxis_n / nrow(matrix_to_plot)
  axis(side = 2, labels = yaxis_labels, at = yaxis_pos)
  
  if(save_file == TRUE){
    dev.off()
  }
}

predict_cpdi_interference = function(bottom_depth,
                                       ave_max_detection_radius,
                                       receiver_depth,
                                       speed_of_sound = 1530,
                                       max_horizontal_dist = 1000,
                                       evaluation_interval = 1,
                                       blanking_interval = .260,
                                       plot = TRUE,
                                       ...){
  ## Models CPDI simulating a tag and fixed depth receiver 
  ## Runs evaluate interference function over a series of distances for a receiver
  ## at a fixed height to calculate a matrix where the rows represent various transmitter
  ## heights  and the columns are distances from the receiver. Values
  ## are the number of conflicting signals for each variable set
  mod_matrix = matrix(0, 
                   1 + ceiling(bottom_depth)/evaluation_interval, 
                   1 + ceiling(max_horizontal_dist)/evaluation_interval)
  for (mat_row in 1:nrow(mod_matrix)){
    transmitter_depth = (mat_row - 1) * evaluation_interval
    for (mat_col in 1:ncol(mod_matrix)){
      transmitter_distance = (mat_col - 1) * evaluation_interval
      transmitter_height_eval_score = 0
      ### Calculate unique multipath arrival times
      score = evaluate_interference(surface_to_tag_distance = transmitter_depth, 
                                       surface_to_receiver_distance = receiver_depth, 
                                       tag_to_receiver_horizontal_distance = transmitter_distance, 
                                       bottom_depth = bottom_depth, 
                                       speed_of_sound = speed_of_sound,
                                       ave_max_detection_radius = ave_max_detection_radius,
                                       blanking_interval = blanking_interval,
                                       plot = FALSE)

      ## If no detection is predicted (interference or out of range), mark it as such
        if(score > 0){
          mod_matrix[mat_row, mat_col] = 0 # zero = no detection
        }else{ ## Otherwise
          mod_matrix[mat_row, mat_col] = 1 # one = detection
        }
      }
  }
  colnames(mod_matrix) = as.character(0:(dim(mod_matrix)[2]-1))
  rownames(mod_matrix) = as.character(0:(dim(mod_matrix)[1]-1))
  if(plot == TRUE){
    plot_heat_map(mod_matrix, ...)
  }
  return(mod_matrix)
}
  

rank_receiver_depths = function(bottom_depth, 
                                      max_horizontal_dist,
                                      ave_max_detection_radius,
                                      evaluation_interval = 1,
                                      speed_of_sound = 1530,
                                      blanking_interval = .260,
                                      tag_depth_range = FALSE){
  ## Models CPDI simulating a fixed depth 
  ## Runs evaluate interference function over a series of distances receiver
    ## heights to calculate a matrix where the rows represent various receiver
    ## height configurations and the rows are distances from the receiver. Values
    ## are the number of conflicting signals for each variable set
  receiver_position_score = c()
  for(receiver_depth in seq(from = 1, to = bottom_depth, by = evaluation_interval)){
    receiver_interference = predict_cpdi_interference(bottom_depth = bottom_depth,
                                                        ave_max_detection_radius = ave_max_detection_radius,
                                                        receiver_depth = receiver_depth,
                                                        speed_of_sound = speed_of_sound,
                                                        max_horizontal_dist = max_horizontal_dist,
                                                        evaluation_interval = evaluation_interval,
                                                        blanking_interval = blanking_interval,
                                                        plot = FALSE)
  
  ## If tag depth range argument is provided, any tag positions not in the depth range are ignored
    if(tag_depth_range[1] != FALSE){
      tag_depth_range = min(tag_depth_range):max(tag_depth_range)
      tag_depth_range = tag_depth_range / evaluation_interval
      ## Note: If depth range is not easily divisible by evaluation interval, depth range is rounded
      tag_depth_range = unique(round(tag_depth_range))
      tag_depth_range = tag_depth_range[tag_depth_range <= nrow(receiver_interference)]
      receiver_interference = receiver_interference[tag_depth_range, ]
      }
    
  ## Receiver Position Score is the total number of positions where detections are thought 
    ## to occur for a receiver at a given height.
    ## A larger receiver position score indicates a large number of detectable positions
  receiver_position_score = c(receiver_position_score, sum(receiver_interference))
  }
  ## Reorder receivers to produce a list of receiver deployment positions from best to worst
  ideal_receiver_heights = order(receiver_position_score, decreasing = TRUE) * evaluation_interval
  return(ideal_receiver_heights)
}


######## Usage Example ######## 
test = function(){
## Determining if, and how far CPDI extends from my receiver
cpdi_presence = predict_cpdi_interference(bottom_depth = 300,
                                    ave_max_detection_radius = 846,
                                    receiver_depth = 299,
                                    speed_of_sound = 1530,
                                    max_horizontal_dist = 1000,
                                    evaluation_interval = 5,
                                    blanking_interval = 0.260,
                                    plot = TRUE,
                                    save_file = FALSE)

## Generating a list of best receiver placements at my study site
ideal_receiver_depth = rank_receiver_depths(bottom_depth = 300,
                                              max_horizontal_dist = 1000,
                                              ave_max_detection_radius = 846,
                                              evaluation_interval = 50,
                                              speed_of_sound = 1530,
                                              blanking_interval = 0.260,
                                              tag_depth_range = 75:400)
}
