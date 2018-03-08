%% RUNGE-KUTTA PROPAGATION OF THE IPW


%% K1
v1_e = zeros(N_traj,1);
v1_n = zeros(N_traj,1);
dphi_e_1 = zeros(Dim_ele,N_traj);
dphi_n_1 = zeros(Dim_nuc,N_traj);

grad_aux_e = e_gradient*phi_e; 
grad_aux_n = n_gradient*phi_n;
PES_e = PES(mesh_n,:).';
PES_n = PES(:,mesh_e);
for i = 1:N_traj
    
    %% PROPAGATE VELOCITIES AND TRAJECTORIES      
    aux_tot = (phi_n(mesh_n(i),:).*phi_e(mesh_e(i),:))*C;
    
    aux_der = (phi_n(mesh_n(i),:).*grad_aux_e(mesh_e(i),:))*C;
    v1_e(i) = imag(aux_der/aux_tot)/e_mass;
        
    aux_der = (phi_e(mesh_e(i),:).*grad_aux_n(mesh_n(i),:))*C;
    v1_n(i) = imag(aux_der/aux_tot)/n_mass;
        
    
    %% PROPAGATE CONDITIONAL WAVEFUNCTIONS    
    PROP_chi = n_cte_kinetic*n_laplacian + diag(PES_n(:,i));
    PROP_psi = e_cte_kinetic*e_laplacian + diag(PES_e(:,i));
    
    dphi_e_1(:,i) = -(1i)*PROP_psi*phi_e(:,i);
    dphi_n_1(:,i) = -(1i)*PROP_chi*phi_n(:,i);

end
%% PROPAGATE Cs
phi_nn = phi_n'*phi_n;
phi_ee = phi_e'*phi_e;
W = phi_nn.*(phi_e'*(Ve.*phi_e))*dx_n*dx_e ...
   +phi_ee.*(phi_n'*(Vn.*phi_n))*dx_n*dx_e;% ...
%   +(phi_n'*(Venn.*phi_n)).*(phi_e'*(Vene.*phi_e))*dx_e*dx_n;
% aux_n = reshape(repmat(conj(phi_n),1,N_traj),Dim_nuc,N_traj,N_traj);
% aux_e = reshape(repmat(conj(phi_e),1,N_traj),Dim_ele,N_traj,N_traj);
% for k = 1:N_traj
%     W(k,:) = W(k,:) + sum((aux_n(:,:,k).*phi_n).*(Ve1n1*(aux_e(:,:,k).*phi_e)),1)*dx_n*dx_e;
% end
for k = 1:N_traj
    W(k,:) = W(k,:) + sum((repmat(conj(phi_n(:,k)),1,N_traj).*phi_n).*(Ve1n1*(repmat(conj(phi_e(:,k)),1,N_traj).*phi_e)),1)*dx_n*dx_e;
end
We = phi_nn.*(phi_e'*(PES_e.*phi_e))*dx_n*dx_e;
Wn = phi_ee.*(phi_n'*(PES_n.*phi_n))*dx_n*dx_e;
M = phi_nn.*phi_ee*dx_n*dx_e;
C_dot_1 = pinv(M)*(-1i*(W-We-Wn)*C);


%% k2
v2_e = zeros(N_traj,1);
v2_n = zeros(N_traj,1);
dphi_e_2 = zeros(Dim_ele,N_traj);
dphi_n_2 = zeros(Dim_nuc,N_traj);

xe_aux = xe + v1_e*dt/2;
xn_aux = xn + v1_n*dt/2;
mesh_e_aux = floor(xe_aux/dx_e - xe_axis(1)/dx_e) + 1;
mesh_n_aux = floor(xn_aux/dx_n - xn_axis(1)/dx_n) + 1;
phi_e_aux = phi_e + dphi_e_1*dt/2;
phi_n_aux = phi_n + dphi_n_1*dt/2;
C_aux = C + C_dot_1*dt/2;
grad_aux_e = e_gradient*phi_e_aux;   
grad_aux_n = n_gradient*phi_n_aux;    
PES_e = PES(mesh_n_aux,:).';
PES_n = PES(:,mesh_e_aux);
for i = 1:N_traj
    
    
    %% PROPAGATE VELOCITIES AND TRAJECTORIES   
    aux_tot = (phi_n_aux(mesh_n_aux(i),:).*phi_e_aux(mesh_e_aux(i),:))*C_aux;
    
    aux_der = (phi_n_aux(mesh_n_aux(i),:).*grad_aux_e(mesh_e_aux(i),:))*C_aux;
    v2_e(i) = imag(aux_der/aux_tot)/e_mass;
    
    aux_der = (phi_e_aux(mesh_e_aux(i),:).*grad_aux_n(mesh_n_aux(i),:))*C_aux;
    v2_n(i) = imag(aux_der/aux_tot)/n_mass;
    
    
    %% PROPAGATE CONDITIONAL WAVEFUNCTIONS   
    PROP_chi = n_cte_kinetic*n_laplacian + diag(PES_n(:,i));
    PROP_psi = e_cte_kinetic*e_laplacian + diag(PES_e(:,i));
    
    dphi_e_2(:,i) = -(1i)*PROP_psi*(phi_e_aux(:,i));
    dphi_n_2(:,i) = -(1i)*PROP_chi*(phi_n_aux(:,i));
        
end
%% PROPAGATE Cs
phi_nn = phi_n_aux'*phi_n_aux;
phi_ee = phi_e_aux'*phi_e_aux;
W = phi_nn.*(phi_e_aux'*(Ve.*phi_e_aux))*dx_n*dx_e ...
   +phi_ee.*(phi_n_aux'*(Vn.*phi_n_aux))*dx_n*dx_e;% ...
%   +(phi_n_aux'*(Venn.*phi_n_aux)).*(phi_e_aux'*(Vene.*phi_e_aux))*dx_e*dx_n;
% aux_n = reshape(repmat(conj(phi_n_aux),1,N_traj),Dim_nuc,N_traj,N_traj);
% aux_e = reshape(repmat(conj(phi_e_aux),1,N_traj),Dim_ele,N_traj,N_traj);
% for k = 1:N_traj
%     W(k,:) = W(k,:) + sum((aux_n(:,:,k).*phi_n_aux).*(Ve1n1*(aux_e(:,:,k).*phi_e_aux)),1)*dx_n*dx_e;
% end
for k = 1:N_traj
    W(k,:) = W(k,:) + sum((repmat(conj(phi_n_aux(:,k)),1,N_traj).*phi_n_aux).*(Ve1n1*(repmat(conj(phi_e_aux(:,k)),1,N_traj).*phi_e_aux)),1)*dx_n*dx_e;
end
We = phi_nn.*(phi_e_aux'*(PES_e.*phi_e_aux))*dx_n*dx_e;
Wn = phi_ee.*(phi_n_aux'*(PES_n.*phi_n_aux))*dx_n*dx_e;
M = phi_nn.*phi_ee*dx_n*dx_e;
C_dot_2 = pinv(M)*(-1i*(W - We - Wn)*C_aux);%gmres(M,(-1i*(W - We - Wn)*(C+C_dot_1*dt/2)));


%% k3
v3_e = zeros(N_traj,1);
v3_n = zeros(N_traj,1);
dphi_e_3 = zeros(Dim_ele,N_traj);
dphi_n_3 = zeros(Dim_nuc,N_traj);

xe_aux = xe + v2_e*dt/2;
xn_aux = xn + v2_n*dt/2;
mesh_e_aux = floor(xe_aux/dx_e - xe_axis(1)/dx_e) + 1;
mesh_n_aux = floor(xn_aux/dx_n - xn_axis(1)/dx_n) + 1;
phi_e_aux = phi_e + dphi_e_2*dt/2;
phi_n_aux = phi_n + dphi_n_2*dt/2;
C_aux = C + C_dot_2*dt/2;
grad_aux_e = e_gradient*phi_e_aux;   
grad_aux_n = n_gradient*phi_n_aux; 
PES_e = PES(mesh_n_aux,:).';
PES_n = PES(:,mesh_e_aux);
for i = 1:N_traj
    
    
    %% PROPAGATE VELOCITIES AND TRAJECTORIES   
    aux_tot = (phi_n_aux(mesh_n_aux(i),:).*phi_e_aux(mesh_e_aux(i),:))*C_aux;
    
    aux_der = (phi_n_aux(mesh_n_aux(i),:).*grad_aux_e(mesh_e_aux(i),:))*C_aux;
    v3_e(i) = imag(aux_der/aux_tot)/e_mass;
    
    aux_der = (phi_e_aux(mesh_e_aux(i),:).*grad_aux_n(mesh_n_aux(i),:))*C_aux;
    v3_n(i) = imag(aux_der/aux_tot)/n_mass;
    
    
    %% PROPAGATE CONDITIONAL WAVEFUNCTIONS    
    PROP_chi = n_cte_kinetic*n_laplacian + diag(PES_n(:,i));
    PROP_psi = e_cte_kinetic*e_laplacian + diag(PES_e(:,i));
    
    dphi_e_3(:,i) = -(1i)*PROP_psi*(phi_e_aux(:,i));
    dphi_n_3(:,i) = -(1i)*PROP_chi*(phi_n_aux(:,i));
    
end
%% PROPAGATE Cs
phi_nn = phi_n_aux'*phi_n_aux;
phi_ee = phi_e_aux'*phi_e_aux;
W = phi_nn.*(phi_e_aux'*(Ve.*phi_e_aux))*dx_n*dx_e ...
   +phi_ee.*(phi_n_aux'*(Vn.*phi_n_aux))*dx_n*dx_e;% ...
%   +(phi_n_aux'*(Venn.*phi_n_aux)).*(phi_e_aux'*(Vene.*phi_e_aux))*dx_e*dx_n;
% aux_n = reshape(repmat(conj(phi_n_aux),1,N_traj),Dim_nuc,N_traj,N_traj);
% aux_e = reshape(repmat(conj(phi_e_aux),1,N_traj),Dim_ele,N_traj,N_traj);
% for k = 1:N_traj
%     W(k,:) = W(k,:) + sum((aux_n(:,:,k).*phi_n_aux).*(Ve1n1*(aux_e(:,:,k).*phi_e_aux)),1)*dx_n*dx_e;
% end
for k = 1:N_traj
    W(k,:) = W(k,:) + sum((repmat(conj(phi_n_aux(:,k)),1,N_traj).*phi_n_aux).*(Ve1n1*(repmat(conj(phi_e_aux(:,k)),1,N_traj).*phi_e_aux)),1)*dx_n*dx_e;
end
We = phi_nn.*(phi_e_aux'*(PES_e.*phi_e_aux))*dx_n*dx_e;
Wn = phi_ee.*(phi_n_aux'*(PES_n.*phi_n_aux))*dx_n*dx_e;
M = phi_nn.*phi_ee*dx_n*dx_e;
C_dot_3 = pinv(M)*(-1i*(W - We - Wn)*C_aux);%gmres(M,(-1i*(W - We - Wn)*(C+C_dot_2*dt/2)));


%% k4
v4_e = zeros(N_traj,1);
v4_n = zeros(N_traj,1);
dphi_e_4 = zeros(Dim_ele,N_traj);
dphi_n_4 = zeros(Dim_nuc,N_traj);

xe_aux = xe + v3_e*dt;
xn_aux = xn + v3_n*dt;
mesh_e_aux = floor(xe_aux/dx_e - xe_axis(1)/dx_e) + 1;
mesh_n_aux = floor(xn_aux/dx_n - xn_axis(1)/dx_n) + 1;
phi_e_aux = phi_e + dphi_e_3*dt;
phi_n_aux = phi_n + dphi_n_3*dt;
C_aux = C + C_dot_3*dt;
grad_aux_e = e_gradient*phi_e_aux;   
grad_aux_n = n_gradient*phi_n_aux; 
PES_e = PES(mesh_n_aux,:).';
PES_n = PES(:,mesh_e_aux);
for i=1:N_traj
    
    %% PROPAGATE VELOCITIES AND TRAJECTORIES    
    aux_tot = (phi_n_aux(mesh_n_aux(i),:).*phi_e_aux(mesh_e_aux(i),:))*C_aux;
    
    aux_der = (phi_n_aux(mesh_n_aux(i),:).*grad_aux_e(mesh_e_aux(i),:))*C_aux;
    v4_e(i) = imag(aux_der/aux_tot)/e_mass;
    
    aux_der = (phi_e_aux(mesh_e_aux(i),:).*grad_aux_n(mesh_n_aux(i),:))*C_aux;
    v4_n(i) = imag(aux_der/aux_tot)/n_mass;
    
    
    
    %% PROPAGATE CONDITIONAL WAVEFUNCTIONS    
    PROP_chi = n_cte_kinetic*n_laplacian + diag(PES_n(:,i));
    PROP_psi = e_cte_kinetic*e_laplacian + diag(PES_e(:,i));
    
    dphi_e_4(:,i) = -(1i)*PROP_psi*(phi_e_aux(:,i));
    dphi_n_4(:,i) = -(1i)*PROP_chi*(phi_n_aux(:,i));
    
end
%% PROPAGATE Cs
phi_nn = phi_n_aux'*phi_n_aux;
phi_ee = phi_e_aux'*phi_e_aux;
W = phi_nn.*(phi_e_aux'*(Ve.*phi_e_aux))*dx_n*dx_e ...
   +phi_ee.*(phi_n_aux'*(Vn.*phi_n_aux))*dx_n*dx_e;% ...
%   +(phi_n_aux'*(Venn.*phi_n_aux)).*(phi_e_aux'*(Vene.*phi_e_aux))*dx_e*dx_n;
% aux_n = reshape(repmat(conj(phi_n_aux),1,N_traj),Dim_nuc,N_traj,N_traj);
% aux_e = reshape(repmat(conj(phi_e_aux),1,N_traj),Dim_ele,N_traj,N_traj);
% for k = 1:N_traj
%     W(k,:) = W(k,:) + sum((aux_n(:,:,k).*phi_n_aux).*(Ve1n1*(aux_e(:,:,k).*phi_e_aux)),1)*dx_n*dx_e;
% end
for k = 1:N_traj
    W(k,:) = W(k,:) + sum((repmat(conj(phi_n_aux(:,k)),1,N_traj).*phi_n_aux).*(Ve1n1*(repmat(conj(phi_e_aux(:,k)),1,N_traj).*phi_e_aux)),1)*dx_n*dx_e;
end
We = phi_nn.*(phi_e_aux'*(PES_e.*phi_e_aux))*dx_n*dx_e;
Wn = phi_ee.*(phi_n_aux'*(PES_n.*phi_n_aux))*dx_n*dx_e;
M = phi_nn.*phi_ee*dx_n*dx_e;
C_dot_4 = pinv(M)*(-1i*(W - We - Wn)*C_aux);%gmres(M,(-1i*(W - We - Wn)*(C+C_dot_3*dt)));


%% EVOLVED CONDITIONAL WAVEFUNCTION AND TRAJECTORIES
phi_e = phi_e + (dt/6)*(dphi_e_1 + 2*dphi_e_2 + 2*dphi_e_3 + dphi_e_4);
phi_n = phi_n + (dt/6)*(dphi_n_1 + 2*dphi_n_2 + 2*dphi_n_3 + dphi_n_4);

xe = xe + (dt/6)*(v1_e + 2*v2_e + 2*v3_e + v4_e);
xn = xn + (dt/6)*(v1_n + 2*v2_n + 2*v3_n + v4_n);

C = C + (dt/6)*(C_dot_1 + 2*C_dot_2 + 2*C_dot_3 + C_dot_4);


