%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Name: Grant Johnson
%Date: 3/31/2022
%Yee algorithm *1D*
%J from relativistic fluid description.

%Notes:
%-1D
% Fluid Scheme: https://ammar-hakim.org/sj/hancock-muscl.html
% E + (v x B) : H&C 
% AST 560 Extra notes: https://ast560.readthedocs.io/en/latest/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all


%Setup field arrays (See PPT)
[N,Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz,Ux,Uy,Uz,Vx,Vy,Vz,grid] = make_grid();

%Stability limit
%%% Initialize memory end %%%


%%% Main code: %%%
%Initial Conditions (at n - 1)
[N,Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz,Ux,Uy,Uz,grid] = IC(N,Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz,Ux,Uy,Uz,grid);

%Make the diagnostic Figure
figure('units','normalized','outerposition',[0 0 0.5 0.5])

%%% Time loop %%%
while(grid.time < grid.t_max)
    
    %Call i/o and diagnostics
    grid = diagnostics(Bx,By,Bz,Ex,Ey,Ez,Jx,Jy,Jz,Ux,Uy,Uz,N,grid);
    
    %Update the gridtime
    grid.time = grid.time + grid.dt;
    
    %Update the iterator
    grid.iter = grid.iter + 1;
  
    %Push U (n - 3/2 -> n - 1/2 (Need E and B both at n), U is on the half-grid
    [Ux,Uy,Uz,N,grid] = fluid_grad_U(Ux,Uy,Uz,N,grid);
    [Ux,Uy,Uz,Vx,Vy,Vz] = BC_J(Ux,Uy,Uz,Vx,Vy,Vz,N,grid);
     N = BC_N(N,Vx,Vy,Vz,grid);
    [Ux,Uy,Uz,Vx,Vy,Vz,grid] = fluid_source_U(Bx,By,Bz,Ex,Ey,Ez,Ux,Uy,Uz,grid);
    [Ux,Uy,Uz,Vx,Vy,Vz] = BC_J(Ux,Uy,Uz,Vx,Vy,Vz,N,grid);
     N = BC_N(N,Vx,Vy,Vz,grid);

    %Advance B field (n - 1 -> n - 1/2) using E|n-1
    [Bx,By,Bz] = push_B(Bx,By,Bz,Ex,Ey,Ez,grid);

    %Deposit the Current (time n - 1/2)
    [Jx,Jy,Jz] = J_deposition(N,Vx,Vy,Vz,grid);
    
    %Advance E field (n-1 -> n) & Periodic Boundaries, needs B, J|n - 1/2
    [Ex,Ey,Ez] = push_E(Bx,By,Bz,Ex,Ey,Ez,Jx,Jy,Jz,grid);
    [Ex,Ey,Ez,Bx,By,Bz,Uy,Uz,Jy,Jz] = BC(Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz,Ux,Uy,Uz,grid);
    
    %Advance B field (n - 1/2 -> n)
    [Bx,By,Bz] = push_B(Bx,By,Bz,Ex,Ey,Ez,grid);
    
end

%%% End Time Loop %%%
%%% End main %%%