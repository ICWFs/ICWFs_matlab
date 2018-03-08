
PD = reshape(abs(phi_ini).^2,Dim_nuc,Dim_ele);
PD_aux = PD;

% PD_e = sum(PD,1)*dx_n;
% PD_n = sum(PD,2)*dx_e;
% abovethreshold_e = (PD_e > threshold).*ones(1,Dim_ele);
% abovethreshold_n = (PD_n > threshold).*ones(Dim_nuc,1);
% abovethreshold = kron(abovethreshold_e,abovethreshold_n);
% PD_aux = abovethreshold;
if montecarlo
    
    %% MONTE CARLO SAMPLING:
    [mesh_n_ini,mesh_e_ini] = mcsampling(PD_aux,N_traj,threshold);
    N_traj = length(mesh_n_ini);
    rand_vals = rand(N_traj,1);
    xe_ini = (mesh_e_ini-1)*dx_e + xe_axis(1) + (rand_vals - 0.5)*dx_e;
    rand_vals = rand(N_traj,1);
    xn_ini = (mesh_n_ini-1)*dx_n + xn_axis(1) + (rand_vals - 0.5)*dx_n;
    
    prob_traj = zeros(N_traj,1);
    for i = 1:N_traj
        prob_traj(i) = PD(mesh_n_ini(i),mesh_e_ini(i));
    end
    
    xe = xe_ini;
    xn = xn_ini;
    xe_old = xe;
    xn_old = xn;
    mesh_e = mesh_e_ini;
    mesh_n = mesh_n_ini;
    four_point_stencil = imag(e_gradent_s*phi_ini(:)./phi_ini(:));
    ve_ini = real(four_point_stencil((mesh_e - 1)*Dim_nuc + mesh_n)/e_mass);
    four_point_stencil = imag(n_gradent_s*phi_ini(:)./phi_ini(:));
    vn_ini = real(four_point_stencil((mesh_e - 1)*Dim_nuc + mesh_n)/n_mass);
    ve = ve_ini;
    vn = vn_ini;
    
else
    
    %% HOMOGENEOUS SAMPLING
    N_traj = 0;
    PD = PD(:);
    for l=1:Dim_ele*Dim_nuc
        if PD(l) >= threshold
            N_traj = N_traj + 1;
        end
    end
    prob_traj = zeros(N_traj,1);
    xe_ini = zeros(N_traj,1);
    xn_ini = zeros(N_traj,1);
    mesh_e_ini = zeros(N_traj,1);
    mesh_n_ini = zeros(N_traj,1);
    
    % Initializing variables:
    N_traj = 0;
    for l=1:Dim_ele*Dim_nuc
        if PD(l) >= threshold
            N_traj = N_traj + 1;
            prob_traj(N_traj) = PD(l);
            
            rand_vals = rand(1);
            xe_ini(N_traj) = xe_grid(l) + (rand_vals - 0.5)*dx_e;
            rand_vals = rand(1);
            xn_ini(N_traj) = xn_grid(l) + (rand_vals - 0.5)*dx_n;
            
            mesh_e_ini(N_traj) = floor(xe_ini(N_traj)/dx_e - xe_axis(1)/dx_e) + 1;
            mesh_n_ini(N_traj) = floor(xn_ini(N_traj)/dx_n - xn_axis(1)/dx_n) + 1;            
        end
    end
    
    xe = xe_ini;
    xn = xn_ini;
    xe_old = xe;
    xn_old = xn;
    mesh_e = mesh_e_ini;
    mesh_n = mesh_n_ini;
    four_point_stencil = imag(e_gradent_s*phi_ini(:)./phi_ini(:));
    ve_ini = real(four_point_stencil((mesh_e - 1)*Dim_nuc + mesh_n)/e_mass);
    four_point_stencil = imag(n_gradent_s*phi_ini(:)./phi_ini(:));
    vn_ini = real(four_point_stencil((mesh_e - 1)*Dim_nuc + mesh_n)/n_mass);
    ve = ve_ini;
    vn = vn_ini;    
    
    fprintf('Number of trajectories: %18.6f \n',N_traj);
    
end





