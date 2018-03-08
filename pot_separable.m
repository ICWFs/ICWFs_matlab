

% Ve = 5E-2*xe_axis.^2;
% Ve = repmat(Ve,1,N_traj);
% Vn = 5E-2*xn_axis.^2;
% Vn = repmat(Vn,1,N_traj);
% Venn = sqrt(lambda)*(xn_axis);
% Venn = repmat(Venn,1,N_traj);
% Vene = sqrt(lambda)*(xe_axis);
% Vene = repmat(Vene,1,N_traj);


Vn1n0 = 1./abs(xn_axis - R_0);
Vn1n2 = 1./abs(xn_axis - R_2);
Ve1n0 = error_function_over_r(full(abs(xe_axis - R_0)),R_c_fixed_l);
Ve1n2 = error_function_over_r(full(abs(xe_axis - R_2)),R_c_fixed_r);

Ve = -Ve1n0 - Ve1n2;
Ve = repmat(Ve,1,N_traj);
Vn = +Vn1n0 + Vn1n2;
Vn = repmat(Vn,1,N_traj);
Ve1n1 = -error_function_over_r(full(abs(xe_grid - xn_grid)),R_c); 

