%Current Deposition
function [Jx,Jy,Jz] = J_deposition(N,Vx,Vy,Vz,grid)

%Sinusoidal Current
%Grab inital size:
%Field at the center of the grid
%Update the current 
%Nx = grid.Nx;
%Nx_half = ceil(Nx/2);
%Jy(Nx_half) = 3.0* sin(grid.time/30);

%Interpolate the quantity
Vx = interp_edge_to_center(Vx,grid); % NEW
N_CC = interp_edge_to_center(N,grid);

%Fluid J, Recall N0 (original) = Num (1/L0) where length 
%transforms like L(gamma) = L0, so Num(lab)/gamma_trans = N0/gamma_origin
% [gamma_center,gamma_edge] = gamma_on_grids(Vx,Vy,Vz,grid);
% Jx = N_CC.*Vx.*grid.e0.*gamma_center;
% Jy = N.*Vy.*grid.e0.*gamma_edge;
% Jz = N.*Vz.*grid.e0.*gamma_edge;
% >??

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
    Jy_laser = sin(2*pi*grid.c*t/grid.laser1.wavelength)*(2.0/(grid.mu_0*grid.c))*grid.laser1.E_max * exp(- ((t - grid.laser1.profile_t_peak)^2) / (grid.laser1.profile_duration^2));

    %1th order interpolation to Jy
    reim = mod(grid.Nx*(grid.laser1.position+grid.xmin)/(grid.xmax - grid.xmin),1);
    nearest_i = floor( grid.Nx*(grid.laser1.position-grid.xmin)/(grid.xmax - grid.xmin) );
    Jy(nearest_i) = (1-reim)*Jy_laser + Jy(nearest_i);
    Jy(nearest_i+1) = (reim)*Jy_laser + Jy(nearest_i);
    
end

end
