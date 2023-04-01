%Intial Conditions
function [N,Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz,Ux,Uy,Uz,grid] = IC(N,Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz,Ux,Uy,Uz,grid)


%Grab inital size:
Nx = grid.Nx;

%Initial Density Profile (Uniform)
if grid.BC_type == "Tunneling through an electron-cyclotron cutoff layer"


    %Constants
    grid.c = 299792458.;
    grid.mu_0 = 4*pi*1e-7;
    grid.eps_0 = 8.85418781762039e-12;
    grid.iter = 1;
    grid.m0 = 9.1093837e-31; %Electrons
    grid.e0 = -1.60217663e-19; %Electrons

    %Additional inputs [ SI ]
    grid.fd = 15e9;
    grid.t0 = 10/grid.fd;
    grid.xc = 0.04;
    grid.xmax = 0.14;
    grid.xmin = 0.0;
    grid.R0 = 5e-3;
    grid.B0 = 0.536;
    grid.Te = 1e-2; %Unused [eV]
    grid.J0 = 1.0;

    grid.dx = (grid.xmax - grid.xmin)/grid.Nx;
    grid.time = 0;
    grid.cfl = 0.98; %clf = udt/dx <= C_max
    grid.dt = 0.98*grid.dx/grid.c;
    grid.t_max = 1.4e-9;
    grid.NT = ceil(grid.t_max/grid.dt);

    %Density
    N = N + 1e17;

    %New grids
    grid.x1 = linspace(grid.xmin,grid.xmax,Nx);
    grid.x2 = interp_edge_to_center(grid.x1);

    %Magnetic Field Profile
    x_Bz = grid.x2;
    Bz = grid.B0 * (grid.R0+grid.xc)./(grid.R0+x_Bz);

end


%Intial Current
%Jx = zeros(Nx-1,1);
%Jy = zeros(Nx,1);
%Jz = zeros(Nx,1);


%Create a photon ( E only)
%for i = (grid.Ey_i_0-1):grid.Ey_i_end
%    Ey(i) = 1.0*sin(3*2*(pi)*(i-1)/(grid.Ey_i_end));
%end

%Lastly Print Stability / Stats
fprintf("Stability Requires: (cdt/dx) we have C = %1.3f of max 1.0\n",grid.c*grid.dt/grid.dx);
fprintf("Grid: Nx: %d, NT: %d\n",grid.Nx,grid.NT);
fprintf("Grid-Spacing: dx: %g, dT: %g\n",grid.dx,grid.dt);

end
