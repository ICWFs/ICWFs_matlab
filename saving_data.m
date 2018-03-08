% SAVING DATA %

time_index = time_index + 1;

%% BOHMIAN TRAJECTORIES
str = strcat('xe_exact',int2str(time_index),'.txt'); xe_exact = fopen(str,'w');
str = strcat('xn_exact',int2str(time_index),'.txt'); xn_exact = fopen(str,'w');
str = strcat('ve_exact',int2str(time_index),'.txt'); ve_exact = fopen(str,'w');
str = strcat('vn_exact',int2str(time_index),'.txt'); vn_exact = fopen(str,'w');
fwrite(xe_exact,xe,'double');
fwrite(xn_exact,xn,'double');
fwrite(ve_exact,ve,'double');
fwrite(vn_exact,vn,'double');
fclose(xe_exact);
fclose(xn_exact);
fclose(ve_exact);
fclose(vn_exact);


%% PROBABILITY DENSITY
str = strcat('exact_wavepacket',int2str(time_index),'.txt'); full_wavepacket = fopen(str,'w');
fwrite(full_wavepacket,abs(phi).^2,'double');
fclose(full_wavepacket);


phi_aux = reshape(phi,Dim_nuc,Dim_ele);
BOEIGSTATE_GR = reshape(BOEIGSTATE_GR,Dim_nuc,Dim_ele);
BOEIGSTATE_EX1 = reshape(BOEIGSTATE_EX1,Dim_nuc,Dim_ele);
gr_comp = zeros(Dim_nuc,1);
ex1_comp = zeros(Dim_nuc,1);
for i = 1:Dim_nuc
    gr_comp(i) = abs(conj(BOEIGSTATE_GR(i,:))*phi_aux(i,:).'*dx_e).^2;
    ex1_comp(i) = abs(conj(BOEIGSTATE_EX1(i,:))*phi_aux(i,:).'*dx_e).^2;
end
str = strcat('gr_comp_exact',int2str(time_index),'.txt'); CWF_n = fopen(str,'w');
fwrite(CWF_n,gr_comp,'double');
fclose(CWF_n);
str = strcat('ex1_comp_exact',int2str(time_index),'.txt'); CWF_n = fopen(str,'w');
fwrite(CWF_n,ex1_comp,'double');
fclose(CWF_n);


red_n = sum(abs(phi_aux).^2,2)*dx_e;
decoh = gr_comp.'*(ex1_comp.*red_n)*dx_n;
str = strcat('decoh_exact',int2str(time_index),'.txt'); CWF_n = fopen(str,'w');
fwrite(CWF_n,decoh,'double');
fclose(CWF_n);

% %% NUCLERA MOMENTUM DISTRIBUTION
dist_P = zeros(Dim_nuc,1);
for j = 1:Dim_nuc
    dist_P(j) = (1/(2*pi))*sum(abs(exp(-1i*pn_axis(j)*xn_axis).'*reshape(phi,Dim_nuc,Dim_ele)*dx_n).^2)*dx_e;
end
str = strcat('dist_P',int2str(time_index),'.txt'); full_wavepacket = fopen(str,'w');
fwrite(full_wavepacket,dist_P,'double');
fclose(full_wavepacket);



%% ENERGY
% total_energy_wave = phi'*(hamiltonian*phi)*dx_e*dx_n;
% str = strcat('total_energy',int2str(time_index),'.txt'); total_energy = fopen(str,'w');
% fwrite(total_energy,total_energy_wave,'double');
% fclose(total_energy);
% 
% potential_energy_wave = phi'*(spdiags(PES(:),0,Dim_ele*Dim_nuc,Dim_ele*Dim_nuc)*phi)*dx_e*dx_n;
% str = strcat('potential_energy',int2str(time_index),'.txt'); total_energy = fopen(str,'w');
% fwrite(total_energy,potential_energy_wave,'double');
% fclose(total_energy);
% 
% kinetic_energy_e = phi'*(e_cte_kinetic*e_lapl_s*phi)*dx_e*dx_n;
% str = strcat('kinetic_energy_e',int2str(time_index),'.txt'); total_energy = fopen(str,'w');
% fwrite(total_energy,kinetic_energy_e,'double');
% fclose(total_energy);
