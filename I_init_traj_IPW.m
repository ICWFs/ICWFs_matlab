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
for alpha = 1:N_traj
    phi_e(:,alpha) = phi_e(:,alpha)/sqrt(sum(abs(phi_e(:,alpha)).^2)*dx_e);
    phi_n(:,alpha) = phi_n(:,alpha)/sqrt(sum(abs(phi_n(:,alpha)).^2)*dx_n);
end


PES = load('PES.txt');
PES = reshape(PES,Dim_nuc,Dim_ele);


% COMPUTE INITIAL Cs
M = (phi_n'*phi_n).*(phi_e'*phi_e)*dx_n*dx_e;
G = zeros(N_traj,1);
for alpha = 1:N_traj
    G(alpha) = phi_n(:,alpha)'*(reshape(phi_ini,Dim_nuc,Dim_ele)*conj(phi_e(:,alpha)))*dx_e*dx_n;
end
C = pinv(M)*G;%gmres(M,G); %M\G
% C = zeros(N_traj,1);
% for alpha = 1:N_traj
%     C(alpha) = 1/N_traj;%PD(mesh_n_ini(alpha),mesh_e_ini(alpha));
% end
% C = C./norm(C);

% C = zeros(N_traj,1);
% G = zeros(N_traj,1);
% for alpha = 1:N_traj
%     G(alpha) = phi_n(:,alpha)'*(reshape(phi_ini,Dim_nuc,Dim_ele)*phi_e(:,alpha))*dx_e*dx_n;
% end
% M = (phi_n'*phi_n)*(phi_e'*phi_e)*dx_n*dx_e;
% [Xsub,idx] = licols(M);
% M_new = vertcat(Xsub);
% C(idx) = M_new\G;

% C = C/norm(C);

% M = (phi_n'*phi_n)*(phi_e'*phi_e)*dx_n*dx_e;
% if rank(M) < floor(N_traj/2)+1
%    singular = true;
% else
%    singular = false;
% end

phi_recon = zeros(Dim_nuc,Dim_ele);
for i = 1:N_traj
    phi_recon = phi_recon + (phi_n(:,i)*phi_e(:,i).')*C(i);
end
phi_recon = phi_recon/sqrt(sum(abs(phi_recon(:)).^2)*dx_e*dx_n);




% figure(15)
% surf(xe_axis,xn_axis,reshape(abs(phi_recon).^2,Dim_nuc,Dim_ele))
% shading interp
% 
% 
% figure(7)
% surf(xe_axis,xn_axis,reshape(abs(phi_ini).^2-abs(phi_recon).^2,Dim_nuc,Dim_ele))
% shading interp
% 
% figure(8)
% surf(xe_axis,xn_axis,reshape(angle(phi_ini)-angle(phi_recon),Dim_nuc,Dim_ele))
% shading interp
% 
% disp(sum(abs(phi_recon(:)).^2)*dx_e*dx_n)
% disp(sum(abs(phi_ini(:)).^2)*dx_e*dx_n)
% 
% figure(17)
% plot(abs(C))

