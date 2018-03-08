clear all;
close all;
%delete *.txt

%% ---------%
%   MODEL   %
% ----------%
lambda = 8E-2;%8E-2;%0.111803;%1;

IPW = true;

n_mass = 1836;%1;%1836;          % 1.67E-27;
e_mass = 1;            % 9.1E-31; 

Dim_ele = 403;
Dim_nuc = 401;
dx_e = 0.35;%0.45;%0.5;
dx_n = 0.045;%0.45;%0.09;        % 5.29E-11 = 1Bohr
dt = 0.1;%0.012;%0.15;             % 2.4E-17s

N_traj = 250;%40;%25;
threshold = 1E-2;

t_end = 32;%1;%32;           %femtoseconds
t_end = t_end*1E-15/(dt*2.4E-17);

num_saved_points = 200;
time_int = floor(t_end/num_saved_points)+1;
num_saved_points = floor(t_end/time_int);
aux = (1:time_int:t_end)*dt*2.4E-17;
save time.txt aux -ascii


% TDSE CONSTANTS %
e_cte_kinetic = -1/(2*e_mass);
n_cte_kinetic = -1/(2*n_mass);

% OPTIONS:
save_data = true;
plot_PES = false;
plot_STATE0 = false;

% Trajectories Propagation:
controlled = false;
direct = true;

% Wave Function Propagation:
euler = false;
runge_kutta = false;
splitting_op = true;
montecarlo = true;


% RK4 %
nvec = 4;
ah(2,1) = 0.5;
ah(3,2) = 0.5;
ah(4,3) = 1.0;
bh(1) = 1.0/6.0;
bh(2) = 1.0/3.0;
bh(3) = 1.0/3.0;
bh(4) = 1.0/6.0;
ch(1) = 0.0;
ch(2) = 0.5;
ch(3) = 0.5;
ch(4) = 1.0;

% SPLIT-OPERATOR %
ne = -(Dim_ele-1)/2:(Dim_ele-1)/2;
nn = -(Dim_nuc-1)/2:(Dim_nuc-1)/2;
ke = 2*ne*pi/(Dim_ele*dx_e);
kn = 2*nn*pi/(Dim_nuc*dx_n);
[Ee,En] = meshgrid(ke,kn);



figure_index = 1;
exact = true;

%% -------- %
% Operators %
% --------- %
B_set_operators %% FUNCTION

%% ------------------- %
% Initial wavefunction %
% -------------------- %
C_set_initial_wavefunction  %% FUNCTION

%% ------------------------- %
% Sampling with Trajectories %
% -------------------------- %
E_initial_trajectory_sampling  

% %% ------------------- %
% % Solution of the TDSE %
% % -------------------- %
H_exact_solution

% %% --------------------%
% % IPW %
% % ---------------------%
H_IPW   %% FUNCTION


% %% --------------------%
% % HERMITIAN %
% % ---------------------%
H_hermitian   %% FUNCTION
