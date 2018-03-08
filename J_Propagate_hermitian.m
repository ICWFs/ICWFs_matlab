
runge_kutta_prop


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

