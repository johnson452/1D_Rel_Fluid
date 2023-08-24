%Current Deposition
function [Jx,Jy,Jz] = J_deposition(N,Vx,Vy,Vz,grid)

%Interpolate the quantity
if grid.solve_type_field == "FDTD"
    Vx = interp_edge_to_center(Vx,grid);
    N_CC = interp_edge_to_center(N,grid);
else
    %Vx = Vx;
    N_CC = N;
end

% No Gamma, when N evolves according to the original eqs
Jx = N_CC.*Vx.*grid.e0;
Jy = N.*Vy.*grid.e0;
Jz = N.*Vz.*grid.e0;


%Inject laser in the domain
if grid.BC_type == "WFA"

    %Laser quantities
    grid.laser1.position    = 9.e-6;               % This point is on the laser plane
    grid.laser1.E_max        = 16.e12;             % Maximum amplitude of the laser field (in V/m)
    grid.laser1.profile_duration = 15.e-15;        % The duration of the laser (in s)
    grid.laser1.profile_t_peak = 30.e-15;          % Time at which the laser reaches its peak (in s)
    grid.laser1.wavelength = 0.8e-6;               % The wavelength of the laser (in m)

    %Laser amplitude (t + dt/2)
    t = grid.time + grid.dt/2;
    ah_hoc_factor = 16*10*2.4e5;
    Jy_laser = ah_hoc_factor*sin(2*pi*grid.c*t/grid.laser1.wavelength)*(2.0/(grid.mu_0*grid.c))*grid.laser1.E_max * exp(- ((t - grid.laser1.profile_t_peak)^2) / (grid.laser1.profile_duration^2));

    %1th order interpolation to Jy
    reim = mod(grid.Nx*(grid.laser1.position-grid.xmin)/(grid.xmax - grid.xmin),1);
    nearest_i = floor( grid.Nx*(grid.laser1.position-grid.xmin)/(grid.xmax - grid.xmin) );
    if nearest_i >= 1
        Jy(nearest_i) = (1-reim)*Jy_laser + Jy(nearest_i);
        Jy(nearest_i+1) = (reim)*Jy_laser + Jy(nearest_i+1);
    end

end

% Antenna at the boundary:
if grid.BC_type == "Propagation into a plasma wave beach"
    Nx = grid.Nx;
    J0 = grid.J0;
    omega = grid.omega;
    t = grid.time;

    %Current at the boundary
    Jy(Nx-1) = J0*sin(omega*t);

end
if grid.BC_type ==  "Tunneling through an electron-cyclotron cutoff layer"
    Nx = grid.Nx;
    J0 = grid.J0;
    fd = grid.fd;
    t = grid.time;
    t0 = grid.t0;

    %Current at the boundary
    Jy(Nx-1) = J0*sin(2*pi*fd*t)*sin(0.5*pi*min(1,t/t0))^2;

end

end
