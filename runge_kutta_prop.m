%% RUNGE-KUTTA PROPAGATION


%% K1
v1_e = zeros(N_traj,1);
v1_n = zeros(N_traj,1);
phi_e_1 = zeros(Dim_ele,N_traj);
phi_n_1 = zeros(Dim_nuc,N_traj);
for i = 1:N_traj
    
    %% PROPAGATE VELOCITIES AND TRAJECTORIES   
    xe_aux = xe(i);
    xn_aux = xn(i);  
    
    four_point_stencil = imag(e_gradient*phi_e(:,i)./phi_e(:,i));
    v1_e(i) = real(four_point_stencil(mesh_e(i))/e_mass);
    
    four_point_stencil = imag(n_gradient*phi_n(:,i)./phi_n(:,i));
    v1_n(i) = real(four_point_stencil(mesh_n(i))/n_mass);      
    
    
    %% PROPAGATE CONDITIONAL WAVEFUNCTIONS
    PES_e = PES(mesh_n(i),:).';
    PES_n = PES(:,mesh_e(i));

    PROP_chi = n_cte_kinetic*n_laplacian + diag(PES_n);
    PROP_psi = e_cte_kinetic*e_laplacian + diag(PES_e);
    
    phi_e_1(:,i) = -(1i)*PROP_psi*phi_e(:,i);
    phi_n_1(:,i) = -(1i)*PROP_chi*phi_n(:,i);

end



%% k2
v2_e = zeros(N_traj,1);
v2_n = zeros(N_traj,1);
phi_e_2 = zeros(Dim_ele,N_traj);
phi_n_2 = zeros(Dim_nuc,N_traj);
phi_e_aux = phi_e + phi_e_1*dt/2;
phi_n_aux = phi_n + phi_n_1*dt/2;
for i = 1:N_traj
 
    
    %% PROPAGATE VELOCITIES AND TRAJECTORIES    
    xe_aux = xe(i) + v1_e(i)*dt/2;
    xn_aux = xn(i) + v1_n(i)*dt/2;    
    mesh_e_aux = floor(xe_aux/dx_e - xe_axis(1)/dx_e) + 1;
    mesh_n_aux = floor(xn_aux/dx_n - xn_axis(1)/dx_n) + 1;
 
    four_point_stencil = imag((e_gradient*phi_e_aux(:,i))./phi_e_aux(:,i));
    v2_e(i) = real(four_point_stencil(mesh_e_aux)/e_mass);
    
    four_point_stencil = imag((n_gradient*phi_n_aux(:,i))./phi_n_aux(:,i));
    v2_n(i) = real(four_point_stencil(mesh_n_aux)/n_mass);  
   
        
    %% PROPAGATE CONDITIONAL WAVEFUNCTIONS
    PES_e = PES(mesh_n_aux,:).';
    PES_n = PES(:,mesh_e_aux);   

    PROP_chi = n_cte_kinetic*n_laplacian + diag(PES_n);
    PROP_psi = e_cte_kinetic*e_laplacian + diag(PES_e);
    
    phi_e_2(:,i) = -(1i)*PROP_psi*phi_e_aux(:,i);
    phi_n_2(:,i) = -(1i)*PROP_chi*phi_n_aux(:,i);
 
end



%% k3
v3_e = zeros(N_traj,1);
v3_n = zeros(N_traj,1);
phi_e_3 = zeros(Dim_ele,N_traj);
phi_n_3 = zeros(Dim_nuc,N_traj);
phi_e_aux = phi_e + phi_e_2*dt/2;
phi_n_aux = phi_n + phi_n_2*dt/2;
for i = 1:N_traj
  
    
    %% PROPAGATE VELOCITIES AND TRAJECTORIES    
    xe_aux = xe(i) + v2_e(i)*dt/2;
    xn_aux = xn(i) + v2_n(i)*dt/2;    
    mesh_e_aux = floor(xe_aux/dx_e - xe_axis(1)/dx_e) + 1;
    mesh_n_aux = floor(xn_aux/dx_n - xn_axis(1)/dx_n) + 1;   

    four_point_stencil = imag((e_gradient*phi_e_aux(:,i))./phi_e_aux(:,i));
    v3_e(i) = real(four_point_stencil(mesh_e_aux)/e_mass);
    
    four_point_stencil = imag((n_gradient*phi_n_aux(:,i))./phi_n_aux(:,i));
    v3_n(i) = real(four_point_stencil(mesh_n_aux)/n_mass);  

    
    %% PROPAGATE CONDITIONAL WAVEFUNCTIONS
    PES_e = PES(mesh_n_aux,:).';
    PES_n = PES(:,mesh_e_aux); 

    PROP_chi = n_cte_kinetic*n_laplacian + diag(PES_n);
    PROP_psi = e_cte_kinetic*e_laplacian + diag(PES_e);
    
    phi_e_3(:,i) = -(1i)*PROP_psi*phi_e_aux(:,i);
    phi_n_3(:,i) = -(1i)*PROP_chi*phi_n_aux(:,i);
          
end


%% k4
v4_e = zeros(N_traj,1);
v4_n = zeros(N_traj,1);
phi_e_4 = zeros(Dim_ele,N_traj);
phi_n_4 = zeros(Dim_nuc,N_traj);
phi_e_aux = phi_e + phi_e_3*dt;
phi_n_aux = phi_n + phi_n_3*dt;
for i=1:N_traj
  
    %% PROPAGATE VELOCITIES AND TRAJECTORIES    
    xe_aux = xe(i) + v3_e(i)*dt;
    xn_aux = xn(i) + v3_n(i)*dt;    
    mesh_e_aux = floor(xe_aux/dx_e - xe_axis(1)/dx_e) + 1;
    mesh_n_aux = floor(xn_aux/dx_n - xn_axis(1)/dx_n) + 1;

    four_point_stencil = imag((e_gradient*phi_e_aux(:,i))./phi_e_aux(:,i));
    v4_e(i) = real(four_point_stencil(mesh_e_aux)/e_mass);
    
    four_point_stencil = imag((n_gradient*phi_n_aux(:,i))./phi_n_aux(:,i));
    v4_n(i) = real(four_point_stencil(mesh_n_aux)/n_mass);  

    
    %% PROPAGATE CONDITIONAL WAVEFUNCTIONS
    PES_e = PES(mesh_n_aux,:).';
    PES_n = PES(:,mesh_e_aux); 

    PROP_chi = n_cte_kinetic*n_laplacian + diag(PES_n);
    PROP_psi = e_cte_kinetic*e_laplacian + diag(PES_e);
    
    phi_e_4(:,i) = -(1i)*PROP_psi*phi_e_aux(:,i);
    phi_n_4(:,i) = -(1i)*PROP_chi*phi_n_aux(:,i);

end


%% EVOLVED CONDITIONAL WAVEFUNCTION AND TRAJECTORIES
phi_e = phi_e + (dt/6)*(phi_e_1 + 2*phi_e_2 + 2*phi_e_3 + phi_e_4);
phi_n = phi_n + (dt/6)*(phi_n_1 + 2*phi_n_2 + 2*phi_n_3 + phi_n_4);

xe = xe + (dt/6)*(v1_e + 2*v2_e + 2*v3_e + v4_e);
xn = xn + (dt/6)*(v1_n + 2*v2_n + 2*v3_n + v4_n);



