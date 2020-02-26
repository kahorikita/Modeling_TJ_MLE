# Modeling_TJ_MLE

Estimate preparation time for target jump at baseline condition.

Main: optim_X_1D_BSL.m
  - Input data: 1:simulated velocity, 2:mean of all subs, 3:bootstrap data, 4:mean of each sub
  - Fitting method: 1:MLE, 2:bads
  
function
- optim_X_1D_BSL.m: main program to run
- get_error_X1D_BSL.m: calculate error between actual and optimal trajectories
- sim_vel_X1D_BSL.m: get optimal trajectory 
- get_trajectory.m: simulate trajectory if choose input 1
- get_bootstrap.m: generate bootstrap data if choose input 3
