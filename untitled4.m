
%% QUANTUM POTENTIAL
rho = zeros(Dim_nuc,Dim_ele);
for i = 1:N_traj_real
    rho = rho + abs(phi_n_1(:,i))*abs(phi_e_1(:,i).');%*prob_traj(i);
end
rho = rho/norm(rho(:));
% phase = zeros(Dim_nuc,Dim_ele);
% for i = 1:N_traj_real
%     phase = phase + (repmat(angle(phi_n_1(:,i)),1,Dim_ele) + repmat(angle(phi_e_1(:,i)).',Dim_nuc,1))*prob_traj(i);
% end
% reconstr_phi = rho;%.*exp(1i*phase);
% reconstr_phi = reconstr_phi/norm(reconstr_phi);

K_e = (reshape(e_cte_kinetic*e_lapl_s*rho(:),Dim_nuc,Dim_ele));
K_n = (reshape(n_cte_kinetic*n_lapl_s*rho(:),Dim_nuc,Dim_ele));

% figure
% mesh(xe_axis,xn_axis,abs(reconstr_phi).^2-abs(phi_ini).^2)%,'or')
% figure
% mesh(xe_axis,xn_axis,abs(phi_ini).^2)%,'-k')
% 
% figure
% mesh(xe_axis,xn_axis,angle(reconstr_phi)-angle(phi_ini))%,'or')
% figure
% mesh(xe_axis,xn_axis,angle(phi_ini))%,'-k')

% pause


for i=1:N_traj_real
    
    PES_n_1 = PES(:,mesh_e(i)); PES_n_1 = diag(PES_n_1);
    PES_e_1 = PES(mesh_n(i),:).'; PES_e_1 = diag(PES_e_1);
    PROP_chi = n_cte_kinetic*n_laplacian + PES_n_1;
    PROP_psi = e_cte_kinetic*e_laplacian + PES_e_1;
    
    phi_e_1_aux = phi_e_1(:,i);
    phi_n_1_aux = phi_n_1(:,i);

    
    %% PROPAGATE TRAJECTORIES    
    ve_old(i) = ve(i);
    vn_old(i) = vn(i);
    
    if (abs(phi_e_1_aux(mesh_e(i)))^2 < 1E-12)
        bohm_ve = 0;
        non_coh_part = non_coh_part + 1;
    else
        four_point_stencil = imag(e_gradient*phi_e_1_aux./phi_e_1_aux);
        bohm_ve = real(four_point_stencil(mesh_e(i))/e_mass);
    end
    
    if (abs(phi_n_1_aux(mesh_n(i)))^2 < 1E-12);
        bohm_vn = 0;
        non_coh_part = non_coh_part + 1;
    else
        four_point_stencil = imag(n_gradient*phi_n_1_aux./phi_n_1_aux);
        bohm_vn = real(four_point_stencil(mesh_n(i))/n_mass);
    end
    
    xe_old(i) = xe(i);
    xn_old(i) = xn(i);
    
    xe(i) = xe(i) + dt*bohm_ve;
    xn(i) = xn(i) + dt*bohm_vn;
    
    ve(i) = bohm_ve;
    vn(i) = bohm_vn;


   
    %% GLOBAL PHASE
%    global_phase_chi = e_cte_kinetic*(e_laplacian*phi_e_1_aux)./phi_e_1_aux ...
%                       + 1i*((e_gradient*phi_e_1_aux)./phi_e_1_aux)*bohm_ve;  
%    global_phase_chi_accum = global_phase_chi_accum + global_phase_chi(mesh_e(i))*dt;
%    global_phase_chi = diag(repmat(global_phase_chi(mesh_e(i)),Dim_nuc,1));
%     
%    global_phase_psi = n_cte_kinetic*(n_laplacian*phi_n_1_aux)./phi_n_1_aux ...
%                       + 1i*((n_gradient*phi_n_1_aux)./phi_n_1_aux)*bohm_vn;
%    global_phase_psi_accum = global_phase_psi_accum + global_phase_chi(mesh_n(i))*dt;   
%    global_phase_psi = diag(repmat(global_phase_psi(mesh_n(i)),Dim_ele,1)); 
    
    
    %% PROPAGATE NUCLEAR CWFs
    ysize = size(phi_n_1_aux);
    kv = zeros(ysize(1),nvec);
    for jj =1:nvec
        y0 = phi_n_1_aux;
        for kk = 1:jj-1
            y0 = y0 + ah(jj,kk)*kv(:,kk);
        end
        kv(:,jj) = -(dt*1i)*(PROP_chi)*y0 -(dt*1i)*K_e(:,mesh_e(i)).*exp(1i*angle(y0)); %global_phase_chi
    end
    y1 = phi_n_1_aux;
    for jj =1:nvec
        y1 = y1 + bh(jj)*kv(:,jj);
    end
    phi_n_1(:,i) = y1;
    
    if norm(phi_n_1(:,i)) > norm_n(i)+0.01 || norm(phi_n_1(:,i)) < norm_n(i)-0.01
        fprintf('norm out of control')
        disp(norm(phi_n_1(i)))
        disp(i)
        stop
    end
    
    %% PROPAGATE ELECTRONIC CWFs
    ysize = size(phi_e_1_aux);
    kv = zeros(ysize(1),nvec);
    for jj =1:nvec
        y0 = phi_e_1_aux;
        for kk = 1:jj-1
            y0 = y0 + ah(jj,kk)*kv(:,kk);
        end
        kv(:,jj) = -(dt*1i)*(PROP_psi)*y0 -(dt*1i)*K_n(mesh_n(i),:).'.*exp(1i*angle(y0)); %global_phase_psi
    end
    y1 = phi_e_1_aux;
    for jj =1:nvec
        y1 = y1 + bh(jj)*kv(:,jj);
    end
    phi_e_1(:,i) = y1;
    
    if norm(phi_e_1(:,i)) > norm_e(i)+0.01 || norm(phi_e_1(:,i)) < norm_e(i)-0.01
        fprintf('norm out of control')
        disp(norm(phi_e_1(:,i)))
        disp(i)
        stop
    end
     
    phi_n_1_phase(:,i) = phi_n_1(:,i).*abs(phi_e_1(mesh_e(i),i).^2)./abs(phi_e_1_ini(mesh_e_ini(i),i)).^2 ;  %exp(-1i*global_phase_chi_accum);
    phi_e_1_phase(:,i) = phi_e_1(:,i).*abs(phi_n_1(mesh_n(i),i).^2)./abs(phi_n_1_ini(mesh_n_ini(i),i)).^2;  %exp(-1i*global_phase_psi_accum);
end


indices = (xe < xe_axis(1));
xe(indices) = xe_old(indices);
stoped_particles = stoped_particles + sum(sum(indices));
indices = (xe > xe_axis(end));
xe(indices) = xe_old(indices);
stoped_particles = stoped_particles + sum(sum(indices));

indices = (xn < xn_axis(1));
xn(indices) = xn_old(indices);
stoped_particles = stoped_particles + sum(sum(indices));
indices = (xn > xn_axis(end));
xn(indices) = xn_old(indices);
stoped_particles = stoped_particles + sum(sum(indices));

mesh_e = floor(xe/dx_e - xe_axis(1)/dx_e) + 1;
mesh_n = floor(xn/dx_n - xn_axis(1)/dx_n) + 1;

