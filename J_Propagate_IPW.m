
%% PROPAGATE CONDITIONAL WAVEFUNCTIONS
runge_kutta_prop_IPW

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
mesh_pn = floor(vn*n_mass/dp_n - pn_axis(1)/dp_n) + 1;


% RECONSTRUCTING THE FULL WAVEFUNCTION
phi_recon = zeros(Dim_nuc,Dim_ele);
for alpha = 1:N_traj
    phi_recon = phi_recon + (phi_n(:,alpha)*phi_e(:,alpha).')*C(alpha);
end
phi_recon = phi_recon/sqrt(sum(abs(phi_recon(:)).^2)*dx_e*dx_n);

% REDEFINING THE CONDITIONAL WAVEFUNCTIONS
% phi_e_new = phi_recon(mesh_n,:).';
% phi_n_new = phi_recon(:,mesh_e); 
% for alpha = 1:N_traj
%    phi_e_new(:,alpha) = phi_e_new(:,alpha)/sqrt(sum(abs(phi_e_new(:,alpha)).^2)*dx_e);
%    phi_n_new(:,alpha) = phi_n_new(:,alpha)/sqrt(sum(abs(phi_n_new(:,alpha)).^2)*dx_n);
% end
% phi_e = phi_e_new;
% phi_n = phi_n_new;
% 
% M = (phi_n'*phi_n).*(phi_e'*phi_e)*dx_n*dx_e;
% G = zeros(N_traj,1);
% for alpha = 1:N_traj
%     G(alpha) = phi_n(:,alpha)'*(reshape(phi_recon,Dim_nuc,Dim_ele)*conj(phi_e(:,alpha)))*dx_e*dx_n;
% end
% C = pinv(M)*G;%gmres(M,G); %M\G

