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
end
