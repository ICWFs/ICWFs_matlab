






%ve = zeros(N_traj,1);
%vn = zeros(N_traj,1);
%xe = zeros(N_traj,1);
%xn = zeros(N_traj,1);
phi_e_1_new = zeros(Dim_ele,N_traj);
phi_n_1_new = zeros(Dim_nuc,N_traj);


%parfor i=1:N_traj
for i=1:N_traj
    
    PROP_chi = 0;
    PROP_psi = 0;
    split_V_n = 0;
    split_K_n = 0;
    split_V_e = 0;
    split_K_e = 0;
    
    %% PROPAGATE TRAJECTORIES    
    phi_e_aux = phi_e_1(:,i);
    phi_n_aux = phi_n_1(:,i);
    xn_aux = xn(i);
    xe_aux = xe(i);
    

    
    % VERLET ALGORITHM
    grad_PES = (L/2 - xn_aux)./abs(L/2 - xn_aux).^3 ...
        - (L/2 + xn_aux)./abs(L/2 + xn_aux).^3 ...
        - (2/sqrt(pi))*(exp(-(abs(xn_aux-xe_aux)/R_c).^2).*(xn_aux-xe_aux))./(R_c*abs(xn_aux-xe_aux).^2) ...
        + (erf(full(abs(xn_aux - xe_aux)/R_c)).*(xn_aux-xe_aux))./abs(xn_aux-xe_aux).^3;
%    nabla_n_Q_e = n_gradient*Q_e(:,mesh_e(i));    
    nabla_n_K_e_class = n_gradient*(K_e_class(:,mesh_e(i)) + A_e(:,mesh_e(i))*ve(i));        
    n_acc_old = -(grad_PES + nabla_n_K_e_class(mesh_n(i)))/n_mass; % + nabla_n_Q_e(mesh_n(i))

   
    
    xn(i) = xn_aux + vn(i)*dt + 0.5*n_acc_old*dt^2;
    xn_aux = xn(i);    
    
    grad_PES = (L/2 - xn_aux)./abs(L/2 - xn_aux).^3 ...
        - (L/2 + xn_aux)./abs(L/2 + xn_aux).^3 ...
        - (2/sqrt(pi))*(exp(-(abs(xn_aux-xe_aux)/R_c).^2).*(xn_aux-xe_aux))./(R_c*abs(xn_aux-xe_aux).^2) ...
        + (erf(full(abs(xn_aux - xe_aux)/R_c)).*(xn_aux-xe_aux))./abs(xn_aux-xe_aux).^3;
%    nabla_n_Q_e = n_gradient*Q_e(:,mesh_e(i));    
    nabla_n_K_e_class = n_gradient*(K_e_class(:,mesh_e(i)) + A_e(:,mesh_e(i))*ve(i));        
    n_acc = -(grad_PES + nabla_n_K_e_class(mesh_n(i)))/n_mass; % + nabla_n_Q_e(mesh_n(i))
    
    vn(i) = vn(i) + 0.5*(n_acc + n_acc_old)*dt;

  
    % Bohmian ve
    four_point_stencil = imag(e_gradient*phi_e_aux./phi_e_aux);
    ve(i) = real(four_point_stencil(mesh_e(i))/e_mass);
%     nabla_e_Q_n = e_gradient*Q_n(mesh_n(i),:).';    
%     nabla_e_K_n_class = e_gradient*K_n_class(mesh_n(i),:).';        
%     e_acc_extra = -(nabla_e_Q_n(mesh_e(i)) + nabla_e_K_n_class(mesh_e(i)))/e_mass;    
    xe(i) = xe_aux + dt*ve(i);% + dt^2*0.5*e_acc_extra;
    
    
%     four_point_stencil = imag(e_gradient*phi_e_aux./phi_e_aux);
%     ve(i) = real(four_point_stencil(mesh_e(i))/e_mass);
%     
%     four_point_stencil = imag(n_gradient*phi_n_aux./phi_n_aux);
%     vn(i) = real(four_point_stencil(mesh_n(i))/n_mass);
%      
%     xe(i) = xe_old(i) + ve(i)*dt;% + nabla_e_Q_n(mesh_e(i))*dt^2;
%     xn(i) = xn_old(i) + vn(i)*dt;% + nabla_n_Q_e(mesh_n(i))*dt^2;
    

    %% PROPAGATE CONDITIONAL WAVEFUNCTIONS
    xe_aux = xe_old(i);
    xn_aux = xn_old(i);
    
    PES_e_1 = 1./abs(xn_aux - R_0) + 1./abs(xn_aux - R_2) ...
        - (error_function_over_r(full(abs(xe_axis - R_0)),R_c_fixed_l) ...
        + error_function_over_r(full(abs(xe_axis - R_2)),R_c_fixed_r) ...
        + error_function_over_r(full(abs(xe_axis - xn_aux)),R_c));
    
    PES_n_1 = 1./abs(xn_axis - R_0) + 1./abs(xn_axis - R_2) ...
        - (error_function_over_r(full(abs(xe_aux - R_0)),R_c_fixed_l) ...
        + error_function_over_r(full(abs(xe_aux - R_2)),R_c_fixed_r) ...
        + error_function_over_r(full(abs(xe_aux - xn_axis)),R_c));
    
   
        
        
    if runge_kutta
        PROP_chi = n_cte_kinetic*n_laplacian + diag(PES_n_1);
        PROP_psi = e_cte_kinetic*e_laplacian + diag(PES_e_1);
    else
        split_V_n = exp(-1i*(dt/2)*(PES_n_1));% + Q_e(:,mesh_e(i)) + S_e(:,mesh_e(i)) + A_e(:,mesh_e(i))*ve(i)));% + 1i*K_im_e_cross(:,mesh_e(i)) + 1i*K_im_e_phase(:,mesh_e(i)) + 1i*A_im_e(:,mesh_e(i))*ve(i))); 
        split_K_n = exp(-1i*dt*kn.^2/(2*n_mass)).';
        
        split_V_e = exp(-1i*(dt/2)*(PES_e_1));% + Q_n(mesh_n(i),:).' + S_n(mesh_n(i),:).' + A_n(mesh_n(i),:).'*vn(i)));% + 1i*K_im_n_cross(mesh_n(i),:).' + 1i*K_im_n_phase(mesh_n(i),:).' + 1i*A_im_n(mesh_n(i),:).'*vn(i)));
        split_K_e = exp(-1i*dt*ke.^2/(2*e_mass)).';
    end
    
    
    if runge_kutta
        % NUCLEI
        ysize = size(phi_n_aux);
        kv = zeros(ysize(1),nvec);
        for jj =1:nvec
            y0 = phi_n_aux;
            for kk = 1:jj-1
                y0 = y0 + ah(jj,kk)*kv(:,kk);
            end
            kv(:,jj) = -(dt*1i)*PROP_chi*y0;
        end
        y1 = phi_n_aux;
        for jj =1:nvec
            y1 = y1 + bh(jj)*kv(:,jj);
        end
        phi_n_1_new(:,i) = y1;

        
        % ELECTRONS
        ysize = size(phi_e_aux);
        kv = zeros(ysize(1),nvec);
        for jj =1:nvec
            y0 = phi_e_aux;
            for kk = 1:jj-1
                y0 = y0 + ah(jj,kk)*kv(:,kk);
            end
            kv(:,jj) = -(dt*1i)*PROP_psi*y0;
        end
        y1 = phi_e_aux;
        for jj =1:nvec
            y1 = y1 + bh(jj)*kv(:,jj);
        end
        phi_e_1_new(:,i) = y1;
 
    else
        % NUCLEI
        phi_n_aux = split_V_n.*phi_n_aux;
        c = fftshift(fft2(ifftshift(phi_n_aux))); 
        c = split_K_n.*c; 
        phi_n_aux = fftshift(ifft2(ifftshift(c))); 
        phi_n_1_new(:,i) = split_V_n.*phi_n_aux;
%%%%%%%%%%%%%%%%
%        phi_n_1_new(:,i) = phi_n_1_new(:,i)/sqrt(sum(abs(phi_n_1_new(:,i)).^2)*dx_n);        
        
        % ELECTRONS
        phi_e_aux = split_V_e.*phi_e_aux;
        c = fftshift(fft2(ifftshift(phi_e_aux))); 
        c = split_K_e.*c; 
        phi_e_aux = fftshift(ifft2(ifftshift(c))); 
        phi_e_1_new(:,i) = split_V_e.*phi_e_aux;
%%%%%%%%%%%%%%%%
%        phi_e_1_new(:,i) = phi_e_1_new(:,i)/sqrt(sum(abs(phi_e_1_new(:,i)).^2)*dx_e);        
    end
end
phi_e_1 = phi_e_1_new;
phi_n_1 = phi_n_1_new;


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

xe_old = xe;
xn_old = xn;

mesh_e = floor(xe/dx_e - xe_axis(1)/dx_e) + 1;
mesh_n = floor(xn/dx_n - xn_axis(1)/dx_n) + 1;
mesh_pn = floor(vn*n_mass/dp_n - pn_axis(1)/dp_n) + 1;

