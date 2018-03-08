% Defining the Initial Wavefunction:

fprintf('\n \n'); fprintf('Defining the Initial Wavefunction... \n');

gr = false; ex1 = true; ex2 = false;

% Nuclear Initial Wavefunction:
n_u0 = -4;%-10;%-4;
n_k0 = 0;
n_dispersion = 1/sqrt(2.85); %0.2
n_cte0 =  sqrt(sqrt(2)/(n_dispersion*sqrt(pi)));
phi_nn = n_cte0*exp(-(xn_axis-n_u0).^2/n_dispersion^2).*exp(1i*n_k0*(xn_axis-n_u0));



% Electronic Initial Wavefunction:
if gr == true
    phi_ee = BOEIGSTATE_GR(:);
elseif ex1 == true
    phi_ee = BOEIGSTATE_EX1(:);
else
    phi_ee = BOEIGSTATE_EX2(:);
end

% e_u0 = -10;
% e_k0 = 0;
% e_dispersion = 1/sqrt(0.2); 
% e_cte0 =  sqrt(sqrt(2)/(e_dispersion*sqrt(pi)));
% phi_ee = e_cte0*exp(-(xe_axis-e_u0).^2/e_dispersion^2).*exp(1i*e_k0*(xe_axis-e_u0));


% Electron-Nuclear Initial Wave Function:
phi_ini = phi_ee.*repmat(phi_nn.',1,Dim_ele).';%kron(phi_ee,phi_nn);
phi_ini = phi_ini/sqrt(sum(abs(phi_ini(:)).^2)*dx_e*dx_n);
phi = phi_ini;

if plot_STATE0
    figure(figure_index)
    surf(xe_axis,xn_axis,reshape(abs(phi).^2,Dim_nuc,Dim_ele))
    shading interp
    grid on
    xlabel('x_{ele}')
    ylabel('x_{nuc}')
    zlabel('Normalized Probability')
  
end
figure_index = figure_index + 1;

if save_data
    aux = phi;
    save phi_en_ini.txt aux -ascii
    clear aux
end
clear phi_ee




