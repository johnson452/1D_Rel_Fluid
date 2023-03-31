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
N_CC = interp_edge_to_center(N);

%Fluid J
Jx = N_CC.*Vx.*grid.e0;
Jy = N.*Vy.*grid.e0;
Jz = N.*Vz.*grid.e0;
end
