# five-arm-maze-analysis
Here, we provide matlab functions for analysing behavioural data acquired in a spatial navigation task. The spatial navigation task was a modified version of the "starmaze"-task (Rondi-Reig et al. 2006, Igloi et al. 2009). 

Matlab R2020b



A. fiveArmMaze_analysisScript

the script and functions were written to analyze x-y-coordinates acquired in a maze with the shape of a star/ five arms (version 190819). 

Input:
The script (fiveArmMaze_analysisScript requires .csv as input-files.
Input-file "trial_results" should contain info about the order of track-files and goal/-task information
Input-files "trackepoint_movment_filename" should contain timestamp ("time"), x-position (pos_x), y-position (pos_z)
 BE AWARE: In case you change the input-format/ structure of input-files,check which coloumns do contain your data and adjust the script accordingly
 
 Output: table.mat, optional testfigure

Script starts with input:
   1. Provide the span of subjects you would like to analyse, start with your first subjectnumber, end with your last subjectnumber (Attention:
   the relevant session folders  (e.g. S001 or S002) has to be in a folder that is named after thesubjectnumber), the script will save the existing IDs into
   group-arrays. The group-arrays will be called later during the actual analysis
   ALTERNATIVE: Provide an array with the numbers of your participants
   Optional: in case you would like to add participants later provide the length/ subjectnumber already analysed

   2. results-folder/ files will be created, choose your name

   3. a m-table called ytable (your table) either exists (due to a previous analysis) or is created which will be saved in the folders provided by you, the table will contain the relevant data

  4. the actual analysis is embedded in three loops, 1st session loop (sN = max number of sessions), 2nd group-loop (gN = max number of groups), 3rd participant loop (called from group-arrays)

  Optional: add analysis functions into loop

  5. data is saved in m-table



B. fiveArmMaze_trajectoryDissimilarity

  The 2nd script was written to compare the shape of two trajectories with different lengths by using dynamical time warping and procrustes 
  
  Input:
  x & y coordinates
  
  Output: 
  Dissimilarity-ratio & optional testfigure
