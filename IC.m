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
    grid.t_max = 5e-8;
    grid.NT = ceil(grid.t_max/grid.dt);

    %Density
    N = N*0.0 + 1e17;

    %New grids
    grid.x1 = linspace(grid.xmin,grid.xmax,Nx);
    grid.x2 = interp_edge_to_center(grid.x1);

    %Magnetic Field Profile
    x_Bz = grid.x2;
    grid.external_Bz = grid.B0 * (grid.R0+grid.xc)./(grid.R0+x_Bz);


    % External quantities
    grid.external_Bx = 0;
    grid.external_By = 0;

    grid.external_Ex = 0;
    grid.external_Ey = 0;
    grid.external_Ez = 0;
end


%Intial Current
%Jx = zeros(Nx-1,1);
%Jy = zeros(Nx,1);
%Jz = zeros(Nx,1);


%Create a photon ( E only)
L = (grid.xmax - grid.xmin);
if grid.BC_type == "Periodic"

    % Frequency (omega = C_frac*wp)
    N = N*0.0 + 1;
    wp = sqrt(grid.mu_0*grid.e0*grid.e0*mean(N)/(grid.m0));
    omega_o_wave = 50.0* wp;


    % Phase error
    k_bar = (1/grid.c)*sqrt(omega_o_wave^2 - wp^2);
    phase = atan(imag(k_bar)/real(k_bar));
    k = real(k_bar);
    ik = imag(k_bar);
    K = sqrt(k^2 + ik^2);

    %Redo the spatial grid:
    grid.lambda = 2*pi/K;
    grid.xmax = 3*grid.lambda;
    grid.xmin = 0.0;
    grid.dx = (grid.xmax - grid.xmin)/grid.Nx;
    grid.time = 0;
    grid.cfl = 0.98; %clf = udt/dx <= C_max
    grid.dt = 0.98*grid.dx/grid.c;
    grid.t_max = 1000*( 1/(omega_o_wave/(2*pi)) );
    grid.NT = ceil(grid.t_max/grid.dt);

    %New grids
    grid.x1 = linspace(grid.xmin,grid.xmax,Nx);
    grid.x2 = interp_edge_to_center(grid.x1);

    % E and B
    E0 = 1.0;
    Ey = E0*sin(2*K*grid.x1/L);
    Bz = (E0*K/omega_o_wave)*sin(2*K*grid.x2/L+phase);

    %Initial Vy, Uy for current density with specified N
    Bz_interp = interp_center_to_edge(Bz);
    Vy = (1/(mean(N)*grid.mu_0*grid.e0))*( ...
        (1/grid.c^2)*(Ey*omega_o_wave) - ...
        (Bz_interp * k) );
    if max(Vy)/grid.c > 1 
        fprintf("Invalid IC\n");
        pause(1000)
    end
    gamma = 1./sqrt(1-Vy.*Vy/(grid.c^2));
    Uy = gamma.*Vy;
end

%Lastly Print Stability / Stats
fprintf("Stability Requires: (cdt/dx) we have C = %1.3f of max 1.0\n",grid.c*grid.dt/grid.dx);
fprintf("Grid: Nx: %d, NT: %d\n",grid.Nx,grid.NT);
fprintf("Grid-Spacing: dx: %g, dT: %g\n",grid.dx,grid.dt);
fprintf("Average: wp: %f\n",sqrt(grid.mu_0*grid.e0*grid.e0*mean(N)/(grid.m0)));

end
