%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Name: Grant Johnson
%Date: 3/31/2022
%Yee algorithm *1D*
%J from relativistic fluid description.

%Notes:
%-1D
% Fluid Scheme: https://ammar-hakim.org/sj/hancock-muscl.html
% E + (v x B) : H&C 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all


%Setup field arrays (See PPT)
grid = make_grid();
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


%Stability limit
fprintf("Stability Requires: (cdt/dx) we have C = %1.3f of max 1.0\n",grid.c*grid.dt/grid.dx);
%%% Initialize memory end %%%


%%% Main code: %%%
%Initial Conditions (at n - 1)
[Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz,Ux,Uy,Uz] = IC(Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz,Ux,Uy,Uz,grid);

%Half push & Periodic Call (at n - 1/2)
[Bx,By,Bz] = push_B(Bx,By,Bz,Ex,Ey,Ez,grid);

%Make the diagnostic Figure
figure('units','normalized','outerposition',[0 0 1 1])

%%% Time loop %%%
while(grid.time < grid.t_max)
    
    %Advance B field (n-1/2 -> n)
    [Bx,By,Bz] = push_B(Bx,By,Bz,Ex,Ey,Ez,grid);
    
    %Call i/o and diagnostics
    grid = diagnostics(Bx,By,Bz,Ex,Ey,Ez,Jx,Jy,Jz,grid);
    
    %Update the gridtime
    grid.time = grid.time + grid.dt;
    
    %Update the iterator
    grid.iter = grid.iter + 1;

    %Updated the fluid U
    [Ux,Uy,Uz,Vx,Vy,Vz] = Rel_Fluids_1D_u_only(Bx,By,Bz,Ex,Ey,Ez,Ux,Uy,Uz,grid);
    
    %Deposite Current
    [Jx,Jy,Jz] = J_deposition(Vx,Vy,Vz,grid);
    
    %Advance E field (n-1 -> n) & Periodic Boundaries
    [Ex,Ey,Ez] = push_E(Bx,By,Bz,Ex,Ey,Ez,Jx,Jy,Jz,grid);
    [Ex,Ey,Ez] = periodic_E(Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz,grid);
    
    %Advance B field (n -> n+1/2)
    [Bx,By,Bz] = push_B(Bx,By,Bz,Ex,Ey,Ez,grid);
    
end

%%% End Time Loop %%%
%%% End main %%%