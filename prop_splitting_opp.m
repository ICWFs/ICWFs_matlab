

%% PROPAGATE TRAJECTORIES
four_point_stencil = imag(e_gradient*phi_e_aux./phi_e_aux);
ve(i) = real(four_point_stencil(mesh_e(i))/e_mass);
xe(i) = xe_old(i) + ve(i)*dt;

four_point_stencil = imag(n_gradient*phi_n_aux./phi_n_aux);
vn(i) = real(four_point_stencil(mesh_n(i))/n_mass);
xn(i) = xn_old(i) + vn(i)*dt;


%% PROPAGATE CONDITIONAL WAVEFUNCTIONS
xe_aux = xe_old(i);
xn_aux = xn_old(i);

PES_e = (1./abs(xn_aux - R_0) + 1./abs(xn_aux - R_2) ...
    - (error_function_over_r(full(abs(xe_axis - R_0)),R_c_fixed_l) ...
    + error_function_over_r(full(abs(xe_axis - R_2)),R_c_fixed_r) ...
    + error_function_over_r(full(abs(xe_axis - xn_aux)),R_c)));

PES_n = (1./abs(xn_axis - R_0) + 1./abs(xn_axis - R_2) ...
    - (error_function_over_r(full(abs(xe_aux - R_0)),R_c_fixed_l) ...
    + error_function_over_r(full(abs(xe_aux - R_2)),R_c_fixed_r) ...
    + error_function_over_r(full(abs(xe_aux - xn_axis)),R_c)));


pot_e = squeeze(pot_orb_e(mesh_n(i),:,:)).'*a(:,i);
pot_n = squeeze(pot_orb_n(mesh_e(i),:,:)).'*A(:,i);

Kin_e = kin_orb_e.'*a(:,i);
Kin_n = kin_orb_n.'*A(:,i);

Kin_corr_e = n_cte_kinetic*c(:,i);
Kin_corr_n = e_cte_kinetic*C(:,i);

Adv_corr_e = 1i*vn(i)*b(:,i);
Adv_corr_n = 1i*ve(i)*B(:,i);





split_V_n = exp(-1i*(dt/2)*(PES_n));% + Q_e(:,mesh_e(i)) + S_e(:,mesh_e(i)) + A_e(:,mesh_e(i))*ve(i)));% + 1i*K_im_e_cross(:,mesh_e(i)) + 1i*K_im_e_phase(:,mesh_e(i)) + 1i*A_im_e(:,mesh_e(i))*ve(i)));
split_K_n = exp(-1i*dt*kn.^2/(2*n_mass)).';

split_V_e = exp(-1i*(dt/2)*(PES_e));% + Q_n(mesh_n(i),:).' + S_n(mesh_n(i),:).' + A_n(mesh_n(i),:).'*vn(i)));% + 1i*K_im_n_cross(mesh_n(i),:).' + 1i*K_im_n_phase(mesh_n(i),:).' + 1i*A_im_n(mesh_n(i),:).'*vn(i)));
split_K_e = exp(-1i*dt*ke.^2/(2*e_mass)).';



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
