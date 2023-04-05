function [N,Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz,Ux,Uy,Uz,Vx,Vy,Vz,grid] = make_grid

%%% Initialize memory %%%
%Setup Grid: (Boundary Grid):
grid.Nx = 200;
grid.xmin = 0;
grid.xmax = 1.0;
grid.dx = (grid.xmax - grid.xmin)/grid.Nx;
grid.time = 0;
grid.dt = 0.1;
grid.t_max = 1000;
grid.NT = ceil(grid.t_max/grid.dt);
grid.Output_interval = 1000;

%Constants
grid.c = 1;
grid.mu_0 = 1;
grid.eps_0 = 1;
grid.iter = 1;
grid.m0 = 1;
grid.e0 = -1; %Electrons


grid.dx = (grid.xmax - grid.xmin)/grid.Nx;
grid.time = 0;
grid.cfl = 0.98; %clf = udt/dx <= C_max
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
grid.x2 = interp_edge_to_center(grid.x1);


%Total Energy
grid.Total_Energy_E_field = zeros(1,grid.NT);
grid.Total_Energy_B_field = zeros(1,grid.NT);
grid.Total_Energy_field = zeros(1,grid.NT);
grid.Total_Energy_ptcls = zeros(1,grid.NT);
grid.Total_Momentum_ptcls = zeros(3,grid.NT);
grid.Total_Momentum_fields = zeros(3,grid.NT);
grid.Total_Momentum_Magnitude_ptcls = zeros(1,grid.NT);
grid.Total_Momentum_Magnitude_fields = zeros(1,grid.NT);

%Main Grids
Nx = grid.Nx;
Ex = zeros(1,Nx-1);
Ey = zeros(1,Nx);
Ez = zeros(1,Nx);
Jx = zeros(1,Nx-1);
Jy = zeros(1,Nx);
Jz = zeros(1,Nx);
Bx = zeros(1,Nx);
By = zeros(1,Nx-1);
Bz = zeros(1,Nx-1);
Ux = zeros(1,Nx-1);
Uy = zeros(1,Nx);
Uz = zeros(1,Nx);
N = zeros(1,Nx);
Vx = zeros(1,Nx-1);
Vy = zeros(1,Nx);
Vz = zeros(1,Nx);

%BCs and ICs extensions
%Standing Photon grid.BC_cond == "Periodic"
grid.BC_cond = "Periodic";
grid.BC_type = "Periodic";
grid.IC_type = grid.BC_type;
grid.problem_name = "Periodic_Photon";

%JE9: Cuttoff
%Turn on options: i.e. If there are externally applied fields, based on 
%the name
%grid.problem_name = "JE9";
%grid.BC_cond = "Non_Periodic";
%grid.BC_type = "Tunneling through an electron-cyclotron cutoff layer";
%grid.IC_type = grid.BC_type;


%Make the file/ delete if already exists:
grid.filename = "Output/out.txt";
if exist(grid.filename, 'file')==2
  delete(grid.filename);
end

end