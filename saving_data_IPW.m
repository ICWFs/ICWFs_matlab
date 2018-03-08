% SAVING DATA %
time_index = time_index + 1;

%% TRAJECTORY POSITIONS AND VELOCITIES
str = strcat('xe_IPW',int2str(time_index),'.txt'); xe_cond = fopen(str,'w');
str = strcat('xn_IPW',int2str(time_index),'.txt'); xn_cond = fopen(str,'w');
str = strcat('ve_IPW',int2str(time_index),'.txt'); ve_cond = fopen(str,'w');
str = strcat('vn_IPW',int2str(time_index),'.txt'); vn_cond = fopen(str,'w');
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
time_dependent_recon_wavepacket_n = sum(abs(phi_recon).^2,2)*dx_e;
time_dependent_recon_wavepacket_e = sum(abs(phi_recon).^2,1)*dx_n;
str = strcat('red_n_IPW_recon',int2str(time_index),'.txt'); CWF_n = fopen(str,'w');
str = strcat('red_e_IPW_recon',int2str(time_index),'.txt'); CWF_e = fopen(str,'w');
fwrite(CWF_n,time_dependent_recon_wavepacket_n,'double');
fwrite(CWF_e,time_dependent_recon_wavepacket_e,'double');
fclose(CWF_n);
fclose(CWF_e);

% time_dependent_wavepacket_n = abs(reshape(phi_n,Dim_nuc,N_traj)./repmat(sqrt(sum(abs(reshape(phi_n,Dim_nuc,N_traj)).^2,1)),Dim_nuc,1)).^2*prob_traj;
% time_dependent_wavepacket_e = abs(reshape(phi_e,Dim_ele,N_traj)./repmat(sqrt(sum(abs(reshape(phi_e,Dim_ele,N_traj)).^2,1)),Dim_ele,1)).^2*prob_traj;
% str = strcat('time_dependent_IPW_wavepacket_n',int2str(time_index),'.txt'); CWF_n = fopen(str,'w');
% str = strcat('time_dependent_IPW_wavepacket_e',int2str(time_index),'.txt'); CWF_e = fopen(str,'w');
% fwrite(CWF_n,time_dependent_wavepacket_n,'double');
% fwrite(CWF_e,time_dependent_wavepacket_e,'double');
% fclose(CWF_n);
% fclose(CWF_e);

phi_ee = phi_e'*phi_e*dx_e;
phi_nn = phi_n'*phi_n*dx_n;
time_dependent_wavepacket_n_from_c = zeros(Dim_nuc,1);
time_dependent_wavepacket_e_from_c = zeros(Dim_ele,1);
aux = zeros(Dim_nuc,1);
phi_n_transf = x_p_transform*phi_n*dx_n;
for i = 1:N_traj
    for k = 1:N_traj
        time_dependent_wavepacket_n_from_c = time_dependent_wavepacket_n_from_c + conj(C(k))*C(i)*conj(phi_n(:,k)).*phi_n(:,i)*phi_ee(k,i);
        time_dependent_wavepacket_e_from_c = time_dependent_wavepacket_e_from_c + conj(C(k))*C(i)*conj(phi_e(:,k)).*phi_e(:,i)*phi_nn(k,i);
        aux = aux + conj(C(k))*C(i)*phi_ee(k,i)*conj(phi_n_transf(:,k)).*(phi_n_transf(:,i));
    end
end
dist_P_n = abs(aux)/(2*pi);

str = strcat('red_n_IPW',int2str(time_index),'.txt'); CWF_n = fopen(str,'w');
str = strcat('red_e_IPW',int2str(time_index),'.txt'); CWF_e = fopen(str,'w');
fwrite(CWF_n,time_dependent_wavepacket_n_from_c,'double');
fwrite(CWF_e,time_dependent_wavepacket_e_from_c,'double');
fclose(CWF_n);
fclose(CWF_e);

str = strcat('momentum_n_IPW',int2str(time_index),'.txt'); CWF_e = fopen(str,'w');
fwrite(CWF_e,dist_P_n,'double');
fclose(CWF_e);

%% REDUCED ADIABATIC QUANTITIES
BOEIGSTATE_GR = reshape(BOEIGSTATE_GR,Dim_nuc,Dim_ele);
BOEIGSTATE_EX1 = reshape(BOEIGSTATE_EX1,Dim_nuc,Dim_ele);
gr_comp = zeros(Dim_nuc,1);
ex1_comp = zeros(Dim_nuc,1);
for i = 1:Dim_nuc
    for alpha = 1:N_traj
        gr_comp(i) = gr_comp(i) + C(alpha)*phi_n(i,alpha)*(conj(BOEIGSTATE_GR(i,:))*phi_e(:,alpha)*dx_e);
        ex1_comp(i) = ex1_comp(i) + C(alpha)*phi_n(i,alpha)*(conj(BOEIGSTATE_EX1(i,:))*phi_e(:,alpha)*dx_e);
    end
end
str = strcat('gr_comp_IPW',int2str(time_index),'.txt'); CWF_n = fopen(str,'w');
fwrite(CWF_n,abs(gr_comp).^2,'double');
fclose(CWF_n);
str = strcat('ex1_comp_IPW',int2str(time_index),'.txt'); CWF_n = fopen(str,'w');
fwrite(CWF_n,abs(ex1_comp).^2,'double');
fclose(CWF_n);


decoh = gr_comp.'*(ex1_comp.*time_dependent_recon_wavepacket_n)*dx_n;
str = strcat('decoh_IPW',int2str(time_index),'.txt'); CWF_n = fopen(str,'w');
fwrite(CWF_n,decoh,'double');
fclose(CWF_n);







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
%     phi_e_1_aux_renorm = phi_e_1(:,alpha)/sqrt(sum(abs(phi_e_1(:,alpha)).^2)*dx_e);               
%     k_energy_e = k_energy_e + (phi_e_1_aux_renorm'*((e_cte_kinetic*e_laplacian)*phi_e_1_aux_renorm)*dx_e).'*1/N_traj;%prob_traj(alpha)/sum(prob_traj);
%     p_energy_e = p_energy_e + (phi_e_1_aux_renorm'*((diag(PES_e_1))*phi_e_1_aux_renorm)*dx_e).'*1/N_traj;%prob_traj(alpha)/sum(prob_traj);    
%     
%     phi_n_1_aux_renorm = phi_n_1(:,alpha)/sqrt(sum(abs(phi_n_1(:,alpha)).^2)*dx_n);               
%     k_energy_n = k_energy_n + (phi_n_1_aux_renorm'*((n_cte_kinetic*n_laplacian)*phi_n_1_aux_renorm)*dx_n).'*1/N_traj;%prob_traj(alpha)/sum(prob_traj);
%     p_energy_n = p_energy_n + (phi_n_1_aux_renorm'*((diag(PES_n_1))*phi_n_1_aux_renorm)*dx_n).'*1/N_traj;%prob_traj(alpha)/sum(prob_traj);    
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
