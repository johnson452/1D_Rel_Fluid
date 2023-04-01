function grid = make_grid

%%% Initialize memory %%%
%Setup Grid: (Boundary Grid):
grid.Nx = 200;
grid.xmin = 0;
grid.xmax = grid.Nx;
grid.dx = 1;
grid.time = 0;
grid.dt = 0.1;
grid.t_max = 100;
grid.NT = ceil(grid.t_max/grid.dt);

%Constants
grid.c = 1;
grid.mu_0 = 1;
grid.eps_0 = 1;
grid.iter = 1;
grid.m0 = 1;
grid.e0 = -1; %Electrons


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

%Offset
%grid.Bx_i_offset = 0.;
%grid.By_i_offset = 0.5;
%grid.Bz_i_offset = 0.5;

%grid.Ex_i_offset = 0.5;
%grid.Ey_i_offset = 0.;
%grid.Ez_i_offset = 0.;

%Build grids: 
% 1: Bx, Ey, Ez, U.. J.. V.. 
% 1: Ex. By, Bz
grid.x1 = linspace(grid.xmin,grid.xmax,Nx);
grid.x2 = interp_edge_to_center(grid.x1);


%Total Energy
grid.Total_Energy_E_field = zeros(1,grid.NT);
grid.Total_Energy_B_field = zeros(1,grid.NT);
grid.Total_Energy_field = zeros(1,grid.NT);

%BCs and ICs extensions
grid.BC_cond = "Non_Periodic";
grid.BC_type = "Tunneling through an electron-cyclotron cutoff layer"; 
grid.IC_type = grid.BC_type;

end