
fprintf('\n \n'); fprintf('Solving the TDSE with Conditional wavefunctions... \n');
time_index = 0;


I_init_traj; 


time_start = clock;
time_step = time_start;
for t = 1:t_end
    if (mod(t,time_int) == 0 || t == 1) 
        fprintf('Estimated Total Time: %18.6f \n',etime(clock(), time_step)*num_saved_points)
        fprintf('Elapsed Time: %18.6f \n',etime(clock(), time_start))
        fprintf('Total number of trajectories: %18.6f \n',N_traj)
        fprintf('\n')   

        time_step = clock;                    
        saving_data_cond
    end        
    J_Propagate_CWF_traj  
end
fprintf('Total Cond Time: %18.6f \n',etime(clock(), time_start));



