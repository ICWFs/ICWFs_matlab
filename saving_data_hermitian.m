% SAVING DATA %
time_index = time_index + 1;

%% TRAJECTORY POSITIONS AND VELOCITIES
str = strcat('xe_hermitian',int2str(time_index),'.txt'); xe_cond = fopen(str,'w');
str = strcat('xn_hermitian',int2str(time_index),'.txt'); xn_cond = fopen(str,'w');
str = strcat('ve_hermitian',int2str(time_index),'.txt'); ve_cond = fopen(str,'w');
str = strcat('vn_hermitian',int2str(time_index),'.txt'); vn_cond = fopen(str,'w');
fwrite(xe_cond,xe,'double');
fwrite(xn_cond,xn,'double');
if t == 1
    fwrite(ve_cond,ve,'double');
    fwrite(vn_cond,vn,'double');
else
    fwrite(ve_cond,(1/6)*(v1_e + 2*v2_e + 2*v3_e + v4_e),'double');
    fwrite(vn_cond,(1/6)*(v1_n + 2*v2_n + 2*v3_n + v4_n),'double');
end
fclose(xe_cond);
fclose(xn_cond);
fclose(ve_cond);
fclose(vn_cond);


%% REDUCED NUCLEAR/ELECTRONIC DENSITIES 
time_dependent_wavepacket_n = abs(reshape(phi_n,Dim_nuc,N_traj)./repmat(sqrt(sum(abs(reshape(phi_n,Dim_nuc,N_traj)).^2,1)),Dim_nuc,1)).^2*prob_traj;
time_dependent_wavepacket_e = abs(reshape(phi_e,Dim_ele,N_traj)./repmat(sqrt(sum(abs(reshape(phi_e,Dim_ele,N_traj)).^2,1)),Dim_ele,1)).^2*prob_traj;
str = strcat('red_n_hermitian',int2str(time_index),'.txt'); CWF_n = fopen(str,'w');
str = strcat('red_e_hermitian',int2str(time_index),'.txt'); CWF_e = fopen(str,'w');
fwrite(CWF_n,time_dependent_wavepacket_n,'double');
fwrite(CWF_e,time_dependent_wavepacket_e,'double');
fclose(CWF_n);
fclose(CWF_e);


%% REDUCED ADIABATIC QUANTITIES
sum_fft_phi_n_wave = zeros(Dim_nuc,Dim_ele);
sum_proj_gr_n_wave = zeros(Dim_nuc,Dim_ele);
sum_proj_ex1_n_wave = zeros(Dim_nuc,Dim_ele);
for ri = 1:Dim_ele
    for alpha = 1:N_traj
        if mesh_e(alpha) == ri
            sum_proj_gr_n_wave(:,ri) = sum_proj_gr_n_wave(:,ri) + abs(BOEIGSTATE_GR(:,ri).*(phi_n(:,alpha)/sqrt(sum(abs(phi_n(:,alpha)).^2,1)*dx_n))).^2*prob_traj(alpha);%1/N_traj;
            sum_proj_ex1_n_wave(:,ri) = sum_proj_ex1_n_wave(:,ri) + abs(BOEIGSTATE_EX1(:,ri).*(phi_n(:,alpha)/sqrt(sum(abs(phi_n(:,alpha)).^2,1)*dx_n))).^2*prob_traj(alpha);%1/N_traj;
            sum_fft_phi_n_wave(:,ri) = sum_fft_phi_n_wave(:,ri) + abs(x_p_transform*(phi_n(:,alpha)/sqrt(sum(abs(phi_n(:,alpha)).^2,1)*dx_n))*dx_n).^2*prob_traj(alpha);%1/N_traj;
        end
    end
end
comp_gr_n_wave = sum(sum_proj_gr_n_wave,2)*dx_e^2;    
comp_ex1_n_wave = sum(sum_proj_ex1_n_wave,2)*dx_e^2;  
dist_P_n_wave = sum(sum_fft_phi_n_wave,2)*dx_e/(2*pi);

comp_gr_e_wave = zeros(Dim_nuc,1);
comp_ex1_e_wave = zeros(Dim_nuc,1);
for Ri = 1:Dim_nuc
    for alpha = 1:N_traj
        if mesh_n(alpha) == Ri
            comp_gr_e_wave(Ri) = comp_gr_e_wave(Ri) + abs(BOEIGSTATE_GR(Ri,:)*(phi_e(:,alpha)/sqrt(sum(abs(phi_e(:,alpha)).^2,1)*dx_e))*dx_e).^2*prob_traj(alpha);%1/N_traj;     %ok
            comp_ex1_e_wave(Ri) = comp_ex1_e_wave(Ri) + abs(BOEIGSTATE_EX1(Ri,:)*(phi_e(:,alpha)/sqrt(sum(abs(phi_e(:,alpha)).^2,1)*dx_e))*dx_e).^2*prob_traj(alpha);%1/N_traj;  %ok
        end
    end
end


% COND_E_WAVE:
str = strcat('comp_gr_e_hermitian',int2str(time_index),'.txt'); CWF_e = fopen(str,'w');
fwrite(CWF_e,comp_gr_e_wave,'double');
fclose(CWF_e);
str = strcat('comp_ex1_e_hermitian',int2str(time_index),'.txt'); CWF_e = fopen(str,'w');
fwrite(CWF_e,comp_ex1_e_wave,'double');
fclose(CWF_e);


% COND_N_WAVE:
str = strcat('comp_gr_n_hermitian',int2str(time_index),'.txt'); CWF_e = fopen(str,'w');
fwrite(CWF_e,comp_gr_n_wave,'double');
fclose(CWF_e);
str = strcat('comp_ex1_n_hermitian',int2str(time_index),'.txt'); CWF_e = fopen(str,'w');
fwrite(CWF_e,comp_ex1_n_wave,'double');
fclose(CWF_e);


% NUCLEAR MOMENTUM DISTRIBUTION
str = strcat('momentum_n_hermitian',int2str(time_index),'.txt'); CWF_e = fopen(str,'w');
fwrite(CWF_e,dist_P_n_wave,'double');
fclose(CWF_e);

% DECOHERENCE
str = strcat('decoh_e_hermitian',int2str(time_index),'.txt'); CWF_e = fopen(str,'w');
fwrite(CWF_e,comp_gr_e_wave.'*(comp_ex1_e_wave.*time_dependent_wavepacket_n),'double');
fclose(CWF_e);

str = strcat('decoh_n_hermitian',int2str(time_index),'.txt'); CWF_e = fopen(str,'w');
fwrite(CWF_e,comp_gr_n_wave.'*(comp_ex1_n_wave.*time_dependent_wavepacket_n),'double');
fclose(CWF_e);



% % ENERGY
% k_energy_e = 0;
% k_energy_n = 0;
% p_energy_e = 0;
% p_energy_n = 0;
% for alpha = 1:N_traj 
%     PES_e_1 = 1./abs(xn(alpha) - R_0) + 1./abs(xn(alpha) - R_2) ...
%             - (error_function_over_r(full(abs(xe_axis - R_0)),R_c_fixed_l) ...
%                + error_function_over_r(full(abs(xe_axis - R_2)),R_c_fixed_r) ...
%                + error_function_over_r(full(abs(xe_axis - xn(alpha))),R_c));    
% 
%     PES_n_1 = 1./abs(xn_axis - R_0) + 1./abs(xn_axis - R_2) ...
%             - (error_function_over_r(full(abs(xe(alpha) - R_0)),R_c_fixed_l) ...
%                + error_function_over_r(full(abs(xe(alpha) - R_2)),R_c_fixed_r) ...
%                + error_function_over_r(full(abs(xe(alpha) - xn_axis)),R_c));    
%            
%            
%     phi_e_aux_renorm = phi_e(:,alpha)/sqrt(sum(abs(phi_e(:,alpha)).^2)*dx_e);               
%     k_energy_e = k_energy_e + (phi_e_aux_renorm'*((e_cte_kinetic*e_laplacian)*phi_e_aux_renorm)*dx_e).'*1/N_traj;%prob_traj(alpha)/sum(prob_traj);
%     p_energy_e = p_energy_e + (phi_e_aux_renorm'*((diag(PES_e_1))*phi_e_aux_renorm)*dx_e).'*1/N_traj;%prob_traj(alpha)/sum(prob_traj);    
%     
%     phi_n_aux_renorm = phi_n(:,alpha)/sqrt(sum(abs(phi_n(:,alpha)).^2)*dx_n);               
%     k_energy_n = k_energy_n + (phi_n_aux_renorm'*((n_cte_kinetic*n_laplacian)*phi_n_aux_renorm)*dx_n).'*1/N_traj;%prob_traj(alpha)/sum(prob_traj);
%     p_energy_n = p_energy_n + (phi_n_aux_renorm'*((diag(PES_n_1))*phi_n_aux_renorm)*dx_n).'*1/N_traj;%prob_traj(alpha)/sum(prob_traj);    
% end
%   
% str = strcat('nuclear_k_energy_cond',int2str(time_index),'.txt'); 
% energy_exact_e = fopen(str,'w');
% fwrite(energy_exact_e,k_energy_n,'double');
% fclose(energy_exact_e);
% 
% str = strcat('nuclear_p_energy_cond',int2str(time_index),'.txt'); 
% energy_exact_e = fopen(str,'w');
% fwrite(energy_exact_e,p_energy_n,'double');
% fclose(energy_exact_e);
% 
% 
% str = strcat('electronic_k_energy_cond',int2str(time_index),'.txt'); 
% energy_exact_e = fopen(str,'w');
% fwrite(energy_exact_e,k_energy_e,'double');
% fclose(energy_exact_e);
% 
% str = strcat('electronic_p_energy_cond',int2str(time_index),'.txt'); 
% energy_exact_e = fopen(str,'w');
% fwrite(energy_exact_e,p_energy_e,'double');
% fclose(energy_exact_e);
