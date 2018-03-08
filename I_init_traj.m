%% ------------------------------------------------------------------%
% INITIALIZING VARIABLES THAT DEPEND ON THE NUMBER OF TRAJECTORIES:  %
% -------------------------------------------------------------------%
xe = xe_ini;
xn = xn_ini;
xe_old = xe_ini;
xn_old = xn_ini;
ve = ve_ini;
vn = vn_ini;
ve_old = ve;
vn_old = vn;
mesh_e = mesh_e_ini;
mesh_n = mesh_n_ini;
mesh_pn = floor(vn*n_mass/dp_n - pn_axis(1)/dp_n) + 1;


stoped_particles = 0; 
non_coh_part = 0; 


phi_ini = reshape(phi_ini,Dim_nuc,Dim_ele);
phi_e = phi_ini(mesh_n,:).';
phi_n = phi_ini(:,mesh_e); 
% RENORMALIZATION
for alpha = 1:N_traj
    phi_e(:,alpha) = phi_e(:,alpha)/sqrt(sum(abs(phi_e(:,alpha)).^2)*dx_e);
    phi_n(:,alpha) = phi_n(:,alpha)/sqrt(sum(abs(phi_n(:,alpha)).^2)*dx_n);
end


PES = load('PES.txt');
PES = reshape(PES,Dim_nuc,Dim_ele);

