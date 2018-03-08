% Defining the Operators:
fprintf('\n \n'); fprintf('Defining the Operators... \n');

%% OPERATORS IN THE ELECTRONIC SPACE:
e_unit_nx = speye(Dim_ele);
e_position = zeros(Dim_ele);
e_gradient = zeros(Dim_ele);
e_laplacian = zeros(Dim_ele);

xh = (Dim_ele-1)/2;
for jj=1:Dim_ele
    e_position(jj,jj) = (jj-1-xh)*dx_e;
end
e_position = sparse(e_position);

for jj=1:Dim_ele
    e_gradient(jj,jj) = 0;
    if (jj < Dim_ele),
        e_gradient(jj,jj+1) = 1/(dx_e)*( +5.0/6.0);
    end
    if (jj > 1),
        e_gradient(jj,jj-1) = 1/(dx_e)*( -5.0/6.0);
    end
    if (jj < Dim_ele - 1),
        e_gradient(jj,jj+2) = 1/(dx_e)*( -5.0/21.0);
    end
    if (jj > 2),
        e_gradient(jj,jj-2) = 1/(dx_e)*( +5.0/21.0);
    end
    if (jj < Dim_ele - 2),
        e_gradient(jj,jj+3) = 1/(dx_e)*( +5.0/84.0);
    end
    if (jj > 3),
        e_gradient(jj,jj-3) = 1/(dx_e)*( -5.0/84.0);
    end
    if (jj < Dim_ele - 3),
        e_gradient(jj,jj+4) = 1/(dx_e)*( -5.0/504.0);
    end
    if (jj > 4),
        e_gradient(jj,jj-4) = 1/(dx_e)*( +5.0/504.0);
    end    
    if (jj < Dim_ele - 4),
        e_gradient(jj,jj+5) = 1/(dx_e)*( +1.0/1260.0);
    end
    if (jj > 5),
        e_gradient(jj,jj-5) = 1/(dx_e)*( -1.0/1260.0);
    end      
end
e_gradient = sparse(e_gradient);

for jj=1:Dim_ele
    e_laplacian(jj,jj) = 1/(dx_e^2)*( -73766.0/25200.0);
    if (jj < Dim_ele),
        e_laplacian(jj,jj+1) = 1/(dx_e^2)*( 5.0/3.0);
    end
    if (jj > 1),
        e_laplacian(jj,jj-1) = 1/(dx_e^2)*( 5.0/3.0);
    end
    if (jj < Dim_ele - 1),
        e_laplacian(jj,jj+2) = 1/(dx_e^2)*( -5.0/21.0);
    end
    if (jj > 2),
        e_laplacian(jj,jj-2) = 1/(dx_e^2)*( -5.0/21.0);
    end
    if (jj < Dim_ele - 2),
        e_laplacian(jj,jj+3) = 1/(dx_e^2)*( 5.0/126.0);
    end
    if (jj > 3),
        e_laplacian(jj,jj-3) = 1/(dx_e^2)*( 5.0/126.0);
    end
    if (jj < Dim_ele - 3),
        e_laplacian(jj,jj+4) = 1/(dx_e^2)*( -5.0/1008.0);
    end
    if (jj > 4),
        e_laplacian(jj,jj-4) = 1/(dx_e^2)*( -5.0/1008.0);
    end
    if (jj < Dim_ele - 4),
        e_laplacian(jj,jj+5) = 1/(dx_e^2)*( +1.0/3150.0);
    end
    if (jj > 5),
        e_laplacian(jj,jj-5) = 1/(dx_e^2)*( +1.0/3150.0);
    end    
end
e_laplacian = sparse(e_laplacian);


%% OPERATORS IN THE NUCLEAR SPACE:
n_unit_nx = speye(Dim_nuc);
n_position = zeros(Dim_nuc);
n_gradient = zeros(Dim_nuc);
n_laplacian = zeros(Dim_nuc);

xh = (Dim_nuc-1)/2;
for jj=1:Dim_nuc
    n_position(jj,jj) = (jj-1-xh)*dx_n;
end
n_position = sparse(n_position);

for jj=1:Dim_nuc
    n_gradient(jj,jj) = 0;
    if (jj < Dim_nuc),
        n_gradient(jj,jj+1) = 1/(dx_n)*( +5.0/6.0);
    end
    if (jj > 1),
        n_gradient(jj,jj-1) = 1/(dx_n)*( -5.0/6.0);
    end
    if (jj < Dim_nuc - 1),
        n_gradient(jj,jj+2) = 1/(dx_n)*( -5.0/21.0);
    end
    if (jj > 2),
        n_gradient(jj,jj-2) = 1/(dx_n)*( +5.0/21.0);
    end
    if (jj < Dim_nuc - 2),
        n_gradient(jj,jj+3) = 1/(dx_n)*( +5.0/84.0);
    end
    if (jj > 3),
        n_gradient(jj,jj-3) = 1/(dx_n)*( -5.0/84.0);
    end
    if (jj < Dim_nuc - 3),
        n_gradient(jj,jj+4) = 1/(dx_n)*( -5.0/504.0);
    end
    if (jj > 4),
        n_gradient(jj,jj-4) = 1/(dx_n)*( +5.0/504.0);
    end    
    if (jj < Dim_nuc - 4),
        n_gradient(jj,jj+5) = 1/(dx_n)*( +1.0/1260.0);
    end
    if (jj > 5),
        n_gradient(jj,jj-5) = 1/(dx_n)*( -1.0/1260.0);
    end      
end
n_gradient = sparse(n_gradient);

for jj=1:Dim_nuc
    n_laplacian(jj,jj) = 1/(dx_n^2)*( -73766.0/25200.0);
    if (jj < Dim_nuc),
        n_laplacian(jj,jj+1) = 1/(dx_n^2)*( 5.0/3.0);
    end
    if (jj > 1),
        n_laplacian(jj,jj-1) = 1/(dx_n^2)*( 5.0/3.0);
    end
    if (jj < Dim_nuc - 1),
        n_laplacian(jj,jj+2) = 1/(dx_n^2)*( -5.0/21.0);
    end
    if (jj > 2),
        n_laplacian(jj,jj-2) = 1/(dx_n^2)*( -5.0/21.0);
    end
    if (jj < Dim_nuc - 2),
        n_laplacian(jj,jj+3) = 1/(dx_n^2)*( 5.0/126.0);
    end
    if (jj > 3),
        n_laplacian(jj,jj-3) = 1/(dx_n^2)*( 5.0/126.0);
    end
    if (jj < Dim_nuc - 3),
        n_laplacian(jj,jj+4) = 1/(dx_n^2)*( -5.0/1008.0);
    end
    if (jj > 4),
        n_laplacian(jj,jj-4) = 1/(dx_n^2)*( -5.0/1008.0);
    end
    if (jj < Dim_nuc - 4),
        n_laplacian(jj,jj+5) = 1/(dx_n^2)*( +1.0/3150.0);
    end
    if (jj > 5),
        n_laplacian(jj,jj-5) = 1/(dx_n^2)*( +1.0/3150.0);
    end    
end
n_laplacian = sparse(n_laplacian);

xe_axis = full(diag(e_position));
xn_axis = full(diag(n_position));

nn = -(Dim_nuc-1)/2:(Dim_nuc-1)/2;
pn_axis = 2*nn.'*pi/(Dim_nuc*dx_n);
dp_n = abs(pn_axis(2)-pn_axis(1));
x_p_transform = exp(-1i*pn_axis*xn_axis.');

%% DEFINING OPERATORS IN THE FULL SPACE:
e_gradent_s = kron(e_gradient,n_unit_nx);
n_gradent_s = kron(e_unit_nx,n_gradient);

e_lapl_s = kron(e_laplacian,n_unit_nx);
n_lapl_s = kron(e_unit_nx,n_laplacian);


%% DEFINING THE FULL HAMILTONIAN:
[xe_grid,xn_grid] = meshgrid(xe_axis,xn_axis);
    
R_c = 5;
R_c_fixed_l = 4;
R_c_fixed_r = 3.1;
R_0 = -9.5;
R_2 = +9.5;
L = 19;

Vn1n0 = 1./abs(xn_grid - R_0);
Vn1n2 = 1./abs(xn_grid - R_2);
Ve1n0 = error_function_over_r(full(abs(xe_grid - R_0)),R_c_fixed_l);
Ve1n2 = error_function_over_r(full(abs(xe_grid - R_2)),R_c_fixed_r);
Ve1n1 = error_function_over_r(full(abs(xe_grid - xn_grid)),R_c); 
PES = spdiags(Vn1n0(:) + Vn1n2(:) - Ve1n0(:) - Ve1n1(:) - Ve1n2(:),0,Dim_ele*Dim_nuc,Dim_ele*Dim_nuc); 


% Ve = 5E-2*xe_grid.^2;
% Vn = 5E-2*xn_grid.^2;
% Ven = lambda*(xe_grid.*xn_grid);
% PES = spdiags(Ve(:) + Vn(:) + Ven(:),0,Dim_ele*Dim_nuc,Dim_ele*Dim_nuc); 



if exact
    if runge_kutta
        hamiltonian = e_cte_kinetic*e_lapl_s + n_cte_kinetic*n_lapl_s + PES;
    else
        % SINCE THE POTENTIAL AND KINETIC OPERATORS ARE DIAGONAL IN THE POSITION AND MOMENTUM SPACE RESPECTIVELY, EXP = EXPM IN THOSE SPACES.
        % DIAG(EXPM(A)) = DIAG(EXP(A))
        hamiltonian = e_cte_kinetic*e_lapl_s + n_cte_kinetic*n_lapl_s + PES;
        split_V = exp(-1i*(dt/2)*full(diag(PES)));
        split_K = exp(-1i*dt*(Ee(:).^2/(2*e_mass)+En(:).^2/(2*n_mass)));
    end
end

%% BO SURFACES & EIGENSTATES
BOPES = zeros(Dim_nuc,1);
BOEIGSTATE_GR = zeros(Dim_nuc,Dim_ele);
BOEIGSTATE_EX1 = zeros(Dim_nuc,Dim_ele);
BOEIGSTATE_EX2 = zeros(Dim_nuc,Dim_ele);

PES_aux = reshape(full(diag(PES)),Dim_nuc,Dim_ele);
for R_aux = 1:Dim_nuc
    PES_eigen_e = PES_aux(R_aux,:);
    PES_eigen_e = sparse(diag(PES_eigen_e(:)));
    
    hamiltonian_aux = e_cte_kinetic*e_laplacian + PES_eigen_e;
    
    [eigv_e,eig_e] = eigs(hamiltonian_aux,6,'SA');
    BOPES(R_aux,1) = eig_e(1,1);
    BOPES(R_aux,2) = eig_e(2,2);
    BOPES(R_aux,3) = eig_e(3,3);
    BOEIGSTATE_GR(R_aux,:) = eigv_e(:,1);
    BOEIGSTATE_EX1(R_aux,:) = eigv_e(:,2);
    BOEIGSTATE_EX2(R_aux,:) = eigv_e(:,3);
    
    BOEIGSTATE_GR(R_aux,:) = BOEIGSTATE_GR(R_aux,:)/sqrt(sum(abs(BOEIGSTATE_GR(R_aux,:)).^2)*dx_e);
    BOEIGSTATE_EX1(R_aux,:) = BOEIGSTATE_EX1(R_aux,:)/sqrt(sum(abs(BOEIGSTATE_EX1(R_aux,:)).^2)*dx_e);
    BOEIGSTATE_EX2(R_aux,:) = BOEIGSTATE_EX2(R_aux,:)/sqrt(sum(abs(BOEIGSTATE_EX2(R_aux,:)).^2)*dx_e);
end
clear eig_e
clear eigv_e
clear PES_aux
clear n_unit_nx
clear e_unit_nx
clear PES_eigen_e
clear hamiltonian_aux

limit = 0.05;
for R_aux = 1:Dim_nuc
    if R_aux > 1
        aux = BOEIGSTATE_GR(R_aux-1,:)*BOEIGSTATE_GR(R_aux,:).';
        if aux < limit; BOEIGSTATE_GR(R_aux,:) = -real(BOEIGSTATE_GR(R_aux,:)); end
        aux = BOEIGSTATE_EX1(R_aux-1,:)*BOEIGSTATE_EX1(R_aux,:).';
        if aux < limit; BOEIGSTATE_EX1(R_aux,:) = -BOEIGSTATE_EX1(R_aux,:); end
        aux = BOEIGSTATE_EX2(R_aux-1,:)*BOEIGSTATE_EX2(R_aux,:).';
        if aux < limit; BOEIGSTATE_EX2(R_aux,:) = -BOEIGSTATE_EX2(R_aux,:); end
    end
end

%% PLOTING AND SAVING EXACT POTENTIAL ENERGY LANDSCAPE
if save_data
    aux = full(xe_axis);
    save eix_e.txt aux -ascii
    
    aux = full(xn_axis);
    save eix_n.txt aux -ascii
    
    aux = full(pn_axis);
    save eix_pn.txt aux -ascii    
    
    aux = full(reshape(diag(PES),Dim_nuc,Dim_ele));
    save PES.txt aux -ascii
    
    aux = BOPES;
    save BOPES.txt aux -ascii
    
    aux = BOEIGSTATE_GR(:);
    save BO_eigst_gr.txt aux -ascii;
    
    aux = BOEIGSTATE_EX1(:);
    save BO_eigst_ex1.txt aux -ascii;
    
    aux = BOEIGSTATE_EX2(:);
    save BO_eigst_ex2.txt aux -ascii;
end

if plot_PES
    figure(figure_index)
    aux = reshape(full(diag(PES)),Dim_nuc,Dim_ele);
    surf(xe_axis,xn_axis,aux)
    shading interp
    xlabel('electronic coordinates (m)')
    ylabel('nuclear coordinates (m)')
    zlabel('Potential Energy (eV)')
    title('Potential Energy Landscape')
    figure_index = figure_index + 1;
    
    for energ_eigen_e = 1:3
        figure(figure_index)
        hold on
        plot(xn_axis,BOPES(:,energ_eigen_e),'-ob')
        xlabel('nuclear coordinates (m)')
        ylabel('Potential Energy (eV)')
        title('BO PES at t=0')
    end
    figure(figure_index+1);mesh(xe_axis,xn_axis,squeeze(real(BOEIGSTATE_GR(:,:))));title('GS BO Conditional Electronic Eigenstate')
    figure(figure_index+2);mesh(xe_axis,xn_axis,squeeze(real(BOEIGSTATE_EX1(:,:))));title('EX1 BO Conditional Electronic Eigenstate')
    figure(figure_index+3);mesh(xe_axis,xn_axis,squeeze(real(BOEIGSTATE_EX2(:,:))));title('EX2 BO Conditional Electronic Eigenstate')
    figure_index = figure_index + 4;
end

clear BOPES
clear PES



