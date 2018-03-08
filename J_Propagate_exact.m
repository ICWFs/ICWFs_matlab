%% Direct Integration Bohmian Trajectory Evolution:
% ------------------------------------------------        
four_point_stencil = imag(e_gradent_s*phi./phi);
ve = real(four_point_stencil((mesh_e - 1)*Dim_nuc + mesh_n)/e_mass);

four_point_stencil = imag(n_gradent_s*phi./phi);
vn = real(four_point_stencil((mesh_e - 1)*Dim_nuc + mesh_n)/n_mass);

xe = xe_old + dt*ve;
xn = xn_old + dt*vn;

indices = (xe < xe_axis(1));
xe(indices) = xe_old(indices);
indices = (xe > xe_axis(end));
xe(indices) = xe_old(indices);

indices = (xn < xn_axis(1));
xn(indices) = xn_old(indices);
indices = (xn > xn_axis(end));
xn(indices) = xn_old(indices);

xe_old = xe;
xn_old = xn;

mesh_e = floor(xe/dx_e - xe_axis(1)/dx_e) + 1;
mesh_n = floor(xn/dx_n - xn_axis(1)/dx_n) + 1;
mesh_pn = floor(vn/dp_n - pn_axis(1)/dp_n) + 1;



%% WaveFunction Propagation:
% --------------------------
if runge_kutta
    ysize = size(phi);
    kv = zeros(ysize(1),nvec);
    for jj =1:nvec
        y0 = phi;
        for kk = 1:jj-1
            y0 = y0 + ah(jj,kk)*kv(:,kk);
        end
        kv(:,jj) = -(dt*1i)*(hamiltonian)*y0;
    end
    y1 = phi;
    for jj =1:nvec
        y1 = y1 + bh(jj)*kv(:,jj);
    end
    phi = y1;
else
    phi = split_V.*phi;                                % Solve dt/2 non-constant part of 2D-LSE
    c = fftshift(fft2(ifftshift(reshape(phi,Dim_nuc,Dim_ele))));  % Take 2D-Fourier transform
    c = split_K.*c(:);                                 % Advance dt in Fourier space
    phi = fftshift(ifft2(ifftshift(reshape(c,Dim_nuc,Dim_ele)))); % Return to physical space
    phi = split_V.*phi(:);                             % Solve dt/2 non-constant part of 2D-LSE
    % phi = phi/sqrt(sum(abs(phi).^2)*dx_e*dx_n); %renormalize
end


if sqrt(sum(abs(phi(:)).^2)*dx_e*dx_n) > 1.01 || sqrt(sum(abs(phi(:)).^2)*dx_e*dx_n) < 0.99
    fprintf('norm out of control')
    disp(sqrt(sum(abs(phi(:)).^2)*dx_e*dx_n))
    stop
end


