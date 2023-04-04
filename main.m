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

%Half push & Periodic Call (at n - 1/2)
[Bx,By,Bz] = push_B(Bx,By,Bz,Ex,Ey,Ez,grid);

%Make the diagnostic Figure
figure('units','normalized','outerposition',[0 0 1 1])

%%% Time loop %%%
while(grid.time < grid.t_max)
    
    %Advance B field (n-1/2 -> n) (E & B now both at n)
    [Bx,By,Bz] = push_B(Bx,By,Bz,Ex,Ey,Ez,grid);
    
    %Call i/o and diagnostics
    grid = diagnostics(Bx,By,Bz,Ex,Ey,Ez,Jx,Jy,Jz,Ux,Uy,Uz,N,grid);
    
    %Update the gridtime
    grid.time = grid.time + grid.dt;
    
    %Update the iterator
    grid.iter = grid.iter + 1;

    %Updated the fluid U (Need E and B both at n), U is on the half-grid
    [Ux,Uy,Uz,Vx,Vy,Vz] = fluid_U(Bx,By,Bz,Ex,Ey,Ez,Ux,Uy,Uz,grid);

    %Fix the BC of the current density:
    [Uy,Uz,Vy,Vz] = BC_J(Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz,Ux,Uy,Uz,Vx,Vy,Vz,N,grid);

    %Deposit the Current
    [Jx,Jy,Jz] = J_deposition(N,Vx,Vy,Vz,grid);
    
    %Advance E field (n-1 -> n) & Periodic Boundaries
    [Ex,Ey,Ez] = push_E(Bx,By,Bz,Ex,Ey,Ez,Jx,Jy,Jz,grid);
    [Ey,Ez,Uy,Uz,Jy,Jz] = BC(Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz,Ux,Uy,Uz,grid);
    
    %Advance B field (n -> n+1/2)
    [Bx,By,Bz] = push_B(Bx,By,Bz,Ex,Ey,Ez,grid);
    
end

%%% End Time Loop %%%
%%% End main %%%