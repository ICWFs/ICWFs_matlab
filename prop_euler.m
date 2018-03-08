
    
%% PROPAGATE TRAJECTORIES
for i=1:N_traj
    
    phi_e_aux = eigst_e*a_e(:,i);
    phi_n_aux = eigst_n*a_n(:,i);
    
    four_point_stencil = imag(e_gradient*phi_e_aux./phi_e_aux);
    ve(i) = real(four_point_stencil(mesh_e(i))/e_mass);
    xe(i) = xe_old(i) + ve(i)*dt;
    
    four_point_stencil = imag(n_gradient*phi_n_aux./phi_n_aux);
    vn(i) = real(four_point_stencil(mesh_n(i))/n_mass);
    xn(i) = xn_old(i) + vn(i)*dt;
    
    
    %% PROPAGATE CONDITIONAL WAVEFUNCTIONS
    pot_e = squeeze(pot_orb_e(mesh_n(i),:,:)).'*a_e(:,i);
    pot_n = squeeze(pot_orb_n(mesh_e(i),:,:)).'*a_n(:,i);
    
    Kin_e = kin_orb_e.'*a_e(:,i);
    Kin_n = kin_orb_n.'*a_n(:,i);
    
    Kin_corr_e = n_cte_kinetic*c_e(:,i);
    Kin_corr_n = e_cte_kinetic*c_n(:,i);
    
    Adv_corr_e = 1i*vn(i)*b_e(:,i);
    Adv_corr_n = 1i*ve(i)*b_n(:,i);
    
    a_e(:,i) = a_e(:,i)-(1i*dt)*(pot_e + Kin_e + Kin_corr_e + Adv_corr_e);
    a_n(:,i) = a_n(:,i)-(1i*dt)*(pot_n + Kin_n + Kin_corr_n + Adv_corr_n);
    
end
