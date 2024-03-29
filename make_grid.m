function [N,Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz,Ux,Uy,Uz,Vx,Vy,Vz,grid] = make_grid

%%% Initialize memory %%%
%[DEFAULT] Setup Grid: (Boundary Grid):
grid.Nx = 2*5120; % Only specified here 
grid.xmin = 0;
grid.xmax = 1.0;
grid.dx = (grid.xmax - grid.xmin)/grid.Nx;
grid.time = 0;
grid.dt = 0.1;
grid.t_max = 200;
grid.NT = ceil(grid.t_max/grid.dt);
grid.Output_interval = 100;
grid.moving_frame = 0;

%[DEFAULT] Dictates Solve type
% Options: "FDTD"; "Muscl";
grid.solve_type_field = "FDTD"; %"Muscl"; %"FDTD"; %"Muscl"; 
grid.laser_envelope_model = "false"; %"false"; %"true";

%[DEFAULT] Constants, updated in IC.m
grid.c = 1;
grid.mu_0 = 1;
grid.eps_0 = 1;
grid.iter = 1;
grid.m0 = 1;
grid.e0 = -1; %Electrons

%[DEFAULT] Grids, updated in IC.m
grid.dx = (grid.xmax - grid.xmin)/grid.Nx;
grid.time = 0;
grid.cfl = 0.45; %clf = udt/dx <= C_max
grid.dt = 0.98*grid.dx/grid.c;
grid.NT = ceil(grid.t_max/grid.dt);


%Restrict:
Nx = grid.Nx;

%Boundaries, B (subdomain):
grid.Bx_i_0 = 1;
grid.Bx_i_end = Nx;

grid.By_i_0 = 1;
grid.By_i_end = Nx-1;

grid.Bz_i_0 = 1;
grid.Bz_i_end = Nx-1;

%Boundaries, E (subdomain):
grid.Ex_i_0 = 1;
grid.Ex_i_end = Nx-1;

grid.Ey_i_0 = 2;
grid.Ey_i_end = Nx-1;

grid.Ez_i_0 = 2;
grid.Ez_i_end = Nx-1;

% If there are externally applied fields
grid.external_fields = "Specified";

%Build grids:
% 1: Bx, Ey, Ez, U.. J.. V..
% 1: Ex. By, Bz
grid.x1 = linspace(grid.xmin,grid.xmax,Nx);
grid.x2 = linspace(grid.xmin+grid.dx/2,grid.xmax-grid.dx/2,Nx-1);


%Main Grids
Nx = grid.Nx;
if grid.solve_type_field == "FDTD"
    Ex = zeros(1,Nx-1);
    Jx = zeros(1,Nx-1);
    By = zeros(1,Nx-1);
    Bz = zeros(1,Nx-1);
else
    Ex = zeros(1,Nx);
    Jx = zeros(1,Nx);
    By = zeros(1,Nx);
    Bz = zeros(1,Nx);
end

Ey = zeros(1,Nx);
Ez = zeros(1,Nx);

Jy = zeros(1,Nx);
Jz = zeros(1,Nx);
Bx = zeros(1,Nx);

Ux = zeros(1,Nx);
Uy = zeros(1,Nx);
Uz = zeros(1,Nx);
N = zeros(1,Nx);
Vx = zeros(1,Nx);
Vy = zeros(1,Nx);
Vz = zeros(1,Nx);

% WFA
% grid.problem_name = "WFA";
% grid.BC_cond = "Non_Periodic";
% grid.BC_type = "WFA";
% grid.WFA_type = 1;
% grid.IC_type = grid.BC_type;

% WFA2: Energy diagnostic 
grid.problem_name = "WFA";
grid.BC_cond = "Non_Periodic";
grid.BC_type = "WFA";
grid.WFA_type = 2;
grid.IC_type = grid.BC_type;

%fluid diagnostic
% grid.problem_name = "Fluid Diagnostic";
% grid.BC_cond = "Non_Periodic";
% grid.BC_type = "fluid_only_diagnostic";
% grid.IC_type = grid.BC_type;


% BCs and ICs extensions
% grid.BC_cond = "Periodic";
% grid.BC_type = "Periodic";
% grid.IC_type = grid.BC_type;
% grid.problem_name = "Periodic_Photon";

% JE9: Cuttoff EC
%Turn on options: i.e. If there are externally applied fields, based on 
%the name
% grid.problem_name = "JE9";
% grid.BC_cond = "Non_Periodic";
% grid.BC_type = "Tunneling through an electron-cyclotron cutoff layer";
% grid.IC_type = grid.BC_type;

% JE8: Cuttoff Plasma Wave beach
% grid.problem_name = "JE8";
% grid.BC_cond = "Non_Periodic";
% grid.BC_type = "Propagation into a plasma wave beach";
% grid.IC_type = grid.BC_type;


%Make the file/ delete if already exists:
grid.filename = "Output/out.txt";
if exist(grid.filename, 'file')==2
  delete(grid.filename);
end

end