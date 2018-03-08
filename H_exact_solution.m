% Solving the TDSE with Runge Kutta:
fprintf('Solving the TDSE... \n');
time_index = 0;

PES = load('PES.txt');
PES = reshape(PES,Dim_nuc,Dim_ele);

time_start = clock;
time_step = time_start;
for t = 1:t_end
    if (mod(t,time_int) == 0 || t == 1)
        fprintf('Estimated Total Time: %18.6f \n',etime(clock(), time_step)*num_saved_points)
        fprintf('Elapsed Time: %18.6f \n',etime(clock(), time_start))
        fprintf('Total number of trajectories: %18.6f \n',N_traj)
        fprintf('\n')  
        
        time_step = clock;
        saving_data
        
%         figure(20)
%         pcolor(xe_axis,xn_axis,reshape(abs(phi).^2,Dim_nuc,Dim_ele))
%         shading interp   
%         hold on
%         scatter(xe,xn,'w')
%         hold off        
    end   
    J_Propagate_exact
end
fprintf('Total Scho Time: %18.6f \n',etime(clock(), time_start));


