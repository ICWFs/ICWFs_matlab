
fprintf('\n \n'); fprintf('Solving the TDSE with Conditional wavefunctions... \n');
time_index = 0;


I_init_traj_IPW
pot_separable 

time_start = clock;
time_step = time_start;
for t = 1:t_end
    if (mod(t,time_int) == 0 || t == 1) 
        fprintf('Estimated Total Time: %18.6f \n',etime(clock(), time_step)*num_saved_points)
        fprintf('Elapsed Time: %18.6f \n',etime(clock(), time_start))
        fprintf('Total number of trajectories: %18.6f \n',N_traj)
        fprintf('\n')   
        
%         disp(singular)

        time_step = clock;                    
        saving_data_IPW
        
%         figure(17)
%         plot(abs(C))
%         
%         figure(18)
%         plot(xe_axis,abs(phi_e).^2)
%         
%         figure(19)
%         plot(xn_axis,abs(phi_n).^2)
% 
%         figure(20)
%         pcolor(xe_axis,xn_axis,abs(phi_recon).^2)
%         shading interp
%         hold on
%         scatter(xe,xn,'w')
%         hold off
        
    end        
    J_Propagate_IPW
end
fprintf('Total Cond Time: %18.6f \n',etime(clock(), time_start));



