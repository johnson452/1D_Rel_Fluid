%Intial Conditions
function [N,Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz,Ux,Uy,Uz,grid] = IC(N,Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz,Ux,Uy,Uz,grid)


%Grab inital size:
Nx = grid.Nx;

%Wake Field Acceleration 1D
if grid.BC_type == "WFA"
    %Constants
    grid.c = 299792458.;
    grid.mu_0 = 4*pi*1e-7;
    grid.eps_0 = 8.85418781762039e-12;
    grid.iter = 1;
    grid.m0 = 9.1093837e-31; %Electrons
    grid.e0 = -1.60217663e-19; %Electrons

    %Additional inputs [ SI ]
    grid.xmin = -100.e-6;
    grid.xmax = 20.e-6;
    grid.xmax0 = grid.xmax;

    grid.L = (grid.xmax - grid.xmin);
    grid.dx = grid.L/grid.Nx;
    grid.time = 0;
    grid.cfl = 0.45; %clf = udt/dx <= C_max
    grid.dt = grid.cfl*grid.dx/grid.c;
    grid.NT = 20000;
    grid.t_max = grid.NT*grid.dt;

    %New grids
    grid.x1 = linspace(grid.xmin,grid.xmax,Nx);
    grid.x2 = linspace(grid.xmin+grid.dx/2,grid.xmax-grid.dx/2,Nx-1);

    % External quantities
    grid.external_Bx = 0;
    grid.external_By = 0;
    grid.external_Bz = 0;

    grid.external_Ex = 0;
    grid.external_Ey = 0;
    grid.external_Ez = 0;

    grid.moving_frame = 1;

    %Laser quantities
    %grid.laser1.profile      = Gaussian;          %assumed
    grid.laser1.position    = 9.e-6;               % This point is on the laser plane
    %grid.laser1.direction    = 0. 0. 1.           % The plane normal direction
    %grid.laser1.polarization = 0. 1. 0.           % The main polarization vector
    grid.laser1.E_max        = 16.e12;             % Maximum amplitude of the laser field (in V/m)
    %grid.laser1.profile_waist = 5.e-6             % The waist of the laser (in m)
    grid.laser1.profile_duration = 15.e-15;        % The duration of the laser (in s)
    grid.laser1.profile_t_peak = 30.e-15;          % Time at which the laser reaches its peak (in s)
    %grid.laser1.profile_focal_distance = 100.e-6  % Focal distance from the antenna (in m)
    grid.laser1.wavelength = 0.8e-6;               % The wavelength of the laser (in m)

    %Density
    N = density_func(grid.x1);
    grid.N0 = N(grid.Nx-1);


end


%Case of a plasma wave hitting a beach: JE8
if grid.BC_type == "Propagation into a plasma wave beach"
    %Constants
    grid.c = 299792458.;
    grid.mu_0 = 4*pi*1e-7;
    grid.eps_0 = 8.85418781762039e-12;
    grid.iter = 1;
    grid.m0 = 9.1093837e-31; %Electrons
    grid.e0 = -1.60217663e-19; %Electrons

    %Additional inputs [ SI ]
    %grid.t0 = 10/grid.fd;
    grid.xc = 0.58;
    grid.xmax = 1.0;
    grid.xmin = 0.0;
    %grid.R0 = 5e-3;
    %grid.B0 = 0.536;
    %grid.Te = 1e-2; %Unused [eV]
    grid.J0 = 1.0e-12;

    grid.L = (grid.xmax - grid.xmin);
    grid.dx = grid.L/grid.Nx;
    grid.time = 0;
    grid.cfl = 0.98; %clf = udt/dx <= C_max
    grid.dt = (1/10)*grid.cfl*grid.dx/grid.c;
    grid.deltat = grid.L/(100*grid.c); %grid.cfl*grid.dx/grid.c;
    grid.t_max = 5e-9;
    grid.NT = ceil(grid.t_max/grid.dt);
    grid.omega = pi/(10*grid.deltat);
    %grid.fd = grid.omega/(pi*2);

    %New grids
    grid.x1 = linspace(grid.xmin,grid.xmax,Nx);
    grid.x2 = linspace(grid.xmin+grid.dx/2,grid.xmax-grid.dx/2,Nx-1);

    %Density
    grid.wpdt = 25*((1-grid.x1)/grid.L).^5;
    grid.wp = grid.wpdt/grid.deltat;
    N = (grid.wp.*grid.wp)*grid.eps_0*grid.m0/(grid.e0*grid.e0);
    grid.N0 = N(1);

    %Overwrite for traveling photon case:
    %N = N*0 + 1;

    %TEMP DIAG
%     plot(grid.x1,grid.x1*0. + grid.omega);
%     hold on
%     plot(grid.x1,grid.wp);
%     hold on
%     plot([0.58,0.58],[0,3e10],":black")
%     legend("w","wp")
    

    % External quantities
    grid.external_Bx = 0;
    grid.external_By = 0;
    grid.external_Bz = 0;

    grid.external_Ex = 0;
    grid.external_Ey = 0;
    grid.external_Ez = 0;

    % Build object for plasma beach:
    if grid.BC_type == "Propagation into a plasma wave beach"
        grid.contour_size = 400;
        grid.temp_iter = 1;
        grid.Ey_t_x = zeros(grid.contour_size,grid.Nx);
    end

end


%EC Cuttoff
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
    grid.N0 = N(1);

    %New grids
    grid.x1 = linspace(grid.xmin,grid.xmax,Nx);
    grid.x2 = linspace(grid.xmin+grid.dx/2,grid.xmax-grid.dx/2,Nx-1);

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

    %Constants
    grid.mu_0 = 4*pi*1e-7;
    grid.eps_0 = 8.85418781762039e-12;
    grid.c = sqrt(1/(grid.mu_0*grid.eps_0));
    grid.iter = 1;
    grid.m0 = 9.1093837e-31; %Electrons
    grid.e0 = -1.60217663e-19; %Electrons

    % Frequency (omega = C_frac*wp)
    N = N*0.0 + 1e17; %1e13
    wp = sqrt(grid.e0*grid.e0*mean(N)/(grid.eps_0*grid.m0));
    omega_o_wave = 2.0* wp; % Critical parameter for dampening
    grid.wp = wp;
    grid.N0 = N(1);


    % Phase error
    k_bar = (1/grid.c)*sqrt(omega_o_wave^2 - wp^2);
    phase = atan(imag(k_bar)/real(k_bar));
    k = real(k_bar);
    ik = imag(k_bar);
    K = sqrt(k^2 + ik^2);
    grid.ik = ik;

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
    grid.x2 = linspace(grid.xmin+grid.dx/2,grid.xmax-grid.dx/2,Nx-1);

    % E and B
    E0 = 1.0;
    Ey = E0*sin(K*grid.x1/L);
    Bz = (E0*K/omega_o_wave)*sin(K*grid.x2/L+phase);

    %Initial Vy, Uy for current density with specified N
    if (0)
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
end

%Lastly Print Stability / Stats
fprintf("Stability Requires: (cdt/dx) we have C = %1.3f of max 1.0\n",grid.c*grid.dt/grid.dx);
fprintf("Grid: Nx: %d, NT: %d\n",grid.Nx,grid.NT);
fprintf("Grid-Spacing: dx: %g, dT: %g\n",grid.dx,grid.dt);
wp_mean = sqrt(grid.e0*grid.e0*mean(N)/(grid.eps_0*grid.m0));
fprintf("Average: wp*dt: %f\n",wp_mean*grid.dt);

%Total Energy
grid.Total_Energy_E_field = zeros(1,grid.NT);
grid.Total_Energy_B_field = zeros(1,grid.NT);
grid.Total_Energy_field = zeros(1,grid.NT);
grid.Total_Energy_ptcls = zeros(1,grid.NT);
grid.Total_Momentum_ptcls = zeros(3,grid.NT);
grid.Total_Momentum_fields = zeros(3,grid.NT);
grid.Total_Momentum_Magnitude_ptcls = zeros(1,grid.NT);
grid.Total_Momentum_Magnitude_fields = zeros(1,grid.NT);

end
