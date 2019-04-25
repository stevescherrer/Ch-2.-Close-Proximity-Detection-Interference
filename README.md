# Ch2-Close-Proximity-Detection-Interference
Analysis and Manuscript associated with Scherrer et. al. 2018

Welcome to the README file for our R implementation of the Mechanistic Model for Predicting CPDI
The source code for this project can be found at the following address:

https://github.com/stevescherrer/CPDI-Submission-Repository/blob/master/Code/Mechanistic%20CPDI%20Model.R

This repository contains scripts used for analysis for the manuscript "Depth-and range-dependent variation in the performance of aquatic telemetry systems: understanding and predicting the susceptibility of acoustic tag–receiver pairs to close proximity detection interference." PeerJ 6 (2018): e4249.
  - 'Analysis/Code/Analysis and Figures - CPDI Paper.R'
  
As well as script files containing functions for implementing predictive CPDI models in both R and matlab.
  - R script: 'Analysis/Code/Mechanistic Model Implemented - R/Mechanistic CPDI Model.R'
  - MatLab file: 'Analysis/Code/Mechanistic Model Implementation - MATLAB/Mechanistic CPDI Model.m'

Recreating manuscript analysis should be relatively straight forward after installing 

As noted in the scripts comment section, there are two main functions for users such as yourself.

1. "Analysis/Code/predict_cpdi_interference()
2. rank_receiver_depths()


Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

Recreating Manuscript Analysis:
Prerequisites
Prior to running code in this repository, you will need to download and install the following software:
	- R (https://www.r-project.org)
  - RStudio (http://rstudio.com)

Once R is installed and running, the following packages that must be installed using the install.packages('package_name') function (note: replace package_name with the name of the package you're trying to install)
This can be accomplished in one of 2 ways. 
  - The easy way: With the script 'Analysis and Figures - CPDI Paper.R' open, use find and replace to find all instances of "## install.packages(" and replace with "(install.packages(", then proceed running the script as normal.
  - individually install the following packages using the install.packages() function
    1.  geosphere
    2.  reshape
    3.  MuMIn
    4.  dplyr
    5.  doParallel
    6.  lubridate
    7.  beepr
    8.  notifyR
    9.  ggplot2
    10. data.table
    11. marmap
    12. mcgv
  
Before running the script, you must also update the 'project_directory' variable on line 14 to reflect your local system path


Implementation of mechanistic CPDI model

R Impementation:
script file: 'Analysis/Code/Mechanistic Model Implemented - R/Mechanistic CPDI Model.R'

As noted in the script's comment section, there are two main functions for users such as yourself.
  1. predict_cpdi_interference()
  2. rank_receiver_depths()

We have also included the function test() with example values to test that the script is working properly. 

1. predict_cpdi_interference
This function takes a number of arguments and produces a matrix with rows corresponding to transmitter depth and columns corresponding to transmitter distance from the receiver. If the argument plot = TRUE, 

The function takes the following arguments:
  — bottom_depth - numeric value. Depth in meters of the study site. In this implementation, it is assumed the study site has a uniform bottom depth. 
  — average_max_detection_radius - numeric value. Maximum distance in meters from a receiver that a tag can be detected. 
  — receiver_depth - numeric value. Depth of receiver  in meters relative to the surface. (ie: If receiver is on rope 3 m above seafloor, receiver depth = bottom depth - 3)
  —  speed_of_sound - numeric value. Speed of sound in enviornment. Assummed to be 1530 m/sec unless specified otherwise
  —  max_horizontal_dist - numeric value. Maximum distance model should simulate
  — evaluation_interval - numeric value.  To reduce computational speed, model can bin search area. By default, model makes a prediction for every 1 meter.  
  — blanking_interval - numeric value. The receiver's blanking interval in seconds (default is 0.260)
  — plot - TRUE or FALSE. Whether or not to produce a plot showing where detections are predicted to occur. Light grey indicates positions where a tag can be heard while dark grey are positions where tag is not heard
  save_file - TRUE or FALSE. If plot argument is true and save_file argument is true, plot will be saved to working directory.
  — ...  - Additonal plot options found in plot heat map function. Notably, save_file = TRUE will save a plot to the working directory
 
2. rank_receiver_depths()
This function takes a number of arguments and generates a list of the optimum heights for placing a receiver based on the number of total  positions that a tag can be detected
Inputs:
  — bottom_depth - numeric value. Depth in meters of study site. 
  — max_horizontal_dist - numeric value . Maximum distance from the receiver model should simulate
  — ave_max_detection_radius - numeric value. Maximum distance in meters from a receiver that a tag can be detected. 
  — evaluation_interval - numeric value.  To reduce computational speed, model can bin search area. By default, model makes a prediction for every 1 meter.  
  —  speed_of_sound - numeric value. Speed of sound in enviornment. Assummed to be 1530 m/sec unless specified otherwise
  — blanking_interval - numeric value. The receiver's blanking interval in seconds (default is 0.260)
  — tag_depth_range - An optional vector of depths a tag may appear of the format min_depth:max_depth or c(min_depth, max_depth). If a fish is known to be present at depths of 120-300 m, it makes no sense to evaluate positions shallower than 120 m.  



Matlab Implementation:
script file: 'Analysis/Code/Mechanistic Model Implementation - MATLAB/Mechanistic CPDI Model.m'

There are two scripts for our MATLAB implementation. 
  1. testCPDI.m
  2. Mechanistic CPDI Model.m

TestCPDI.m tests the parameter set provided for CPDI.  It relies on MOI_XYZ.m to provide the underlying model. 

Both functions take the following arguments.

—m = a vector of model parameters in the following order.
	xSource - Numeric Value. The horizontal position of the source transmitter (tag) in meters
	ySource - 0. Is not used in this model, should be set to zero 
	zSource - Numeric Value - Corresponds to the tag depth, measured positive downward from surface. Ex: To measure a tag 1m above a seafloor of 300 m depth use 299. In meters
	cBartSource - 0. Is not relevant for implementing this model and can be set to zero
	cBardOffset - 0. Is not relevant for implementing this model and can be set to zero
	xPhone - Numeric value. The horizontal position of the hydrophone (receiver’s grid location x axis)
	yPhone - 0.  Is not relevant for implementing this model and should be set to zero
	zPhone - Numeric value. Corresponds to the receiver depth, measured positive downward from surface. Ex: To measure a receiver 1m above a seafloor of 300 m depth use 299.	
	dWater - Numeric value. The depth of the water column
	cWater - Numeric value The sound speed in meters/second
	cConstant - 0, scaling value used when this propagation model is used for localization. Not important for this model

arrivePath      (binary vector indicating which arrival paths arrival times must be calculated). This is a vector of 20 elements corresponding to the first 20 arrival paths. The direct path is the first element, followed by a single bottom reflection(B). The third element is a single surface reflection(S). The fourth element is reflected twice, first off the bottom and then surface (BS), the fifth is reflected twice, first the surface bottom reflection (SB), The sixth is reflected three times, the bottom surface bottom reflection (BSB) followed by the seventh element reflected three times first surface, then bottom, then surface reflections (SBS)… etc… 

For example, a water column of 220 m. Tag is positioned 20 m above seafloor (200 m depth) and Receiver is located 25 m above seafloor (195 m depth). Sound speed in water is 1530 m/s. Receiver is said to be at x = 0 with tag a relative distance of 122 m away.

ex: m = [122, 0, 200, 0, 0, 0, 0, 195, 220, 1530, 0] 

If we were to determine CPDI from the direct path and first 4 multi paths we would declare arrivePath as such:
ex: arrivePath = [1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]

After loading both function files into our path, we run the function: testCPDI(m, arrivePath)

To Run the function over a set of tag depths, put function in a for loop and update zSource with each loop iteration. 

To determine best receiver placement, nest for loops so inner loop 
