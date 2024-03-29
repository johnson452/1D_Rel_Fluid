function [N,Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz,Ux,Uy,Uz,Vx,Vy,Vz,grid] = field_solve(...
    N,Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz,Ux,Uy,Uz,Vx,Vy,Vz,grid)


%Push U (n - 3/2 -> n - 1/2 (Need E and B both at n - 1), U is on the half-time-grid
if (grid.solve_type_field ~= "None")
    [Ux,Uy,Uz,Vx,Vy,Vz,grid] = fluid_source_U(Bx,By,Bz,Ex,Ey,Ez,Ux,Uy,Uz,grid);
    [Ux,Uy,Uz,Vx,Vy,Vz] = BC_J(Ux,Uy,Uz,Vx,Vy,Vz,N,grid);
    N = BC_N(N,Vx,Vy,Vz,grid);
end
%[Ux,Uy,Uz,Vx,Vy,Vz,N,grid] = fluid_grad_U(Ux,Uy,Uz,N,grid);
%[Ux,Uy,Uz,Vx,Vy,Vz,N,grid] = fluid_grad_U_FD(Ux,Uy,Uz,N,grid);
%[Ux,Uy,Uz,Vx,Vy,Vz,N,grid] = fluid_grad_U_FD_Upwind(Ux,Uy,Uz,N,grid);
%[Ux,Uy,Uz,Vx,Vy,Vz,N,grid] = fluid_grad_U_FCT(Ux,Uy,Uz,N,grid);
%[Ux,Uy,Uz,Vx,Vy,Vz,N,grid] = fluid_grad_U_Isothermal(Ux,Uy,Uz,N,grid);
%[Ux,Uy,Uz,Vx,Vy,Vz,N,grid] = fluid_grad_U_fwaves(Ux,Uy,Uz,N,grid);
%[Ux,Uy,Uz,Vx,Vy,Vz,N,grid] = fluid_grad_U_fwaves_Eulderink(Ux,Uy,Uz,N,grid);
%[Ux,Uy,Uz,Vx,Vy,Vz,N,grid] = fluid_grad_U_fwaves_eigen(Ux,Uy,Uz,N,grid);
[Ux,Uy,Uz,Vx,Vy,Vz,N,grid] = fluid_grad_U_prims(Ux,Uy,Uz,N,grid); % BEST
[Ux,Uy,Uz,Vx,Vy,Vz] = BC_J(Ux,Uy,Uz,Vx,Vy,Vz,N,grid);
N = BC_N(N,Vx,Vy,Vz,grid);


if grid.solve_type_field == "FDTD"

    %Advance B field (n - 1 -> n - 1/2) using E|n-1
    [Bx,By,Bz] = push_B(Bx,By,Bz,Ex,Ey,Ez,grid);
    [Ex,Ey,Ez,Bx,By,Bz,Uy,Uz,~,~] = BC(Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz,Ux,Uy,Uz,grid);

    %Deposit the Current (time n - 1/2)
    [Jx,Jy,Jz] = J_deposition(N,Vx,Vy,Vz,grid);

    %Advance E field (n-1 -> n) & Periodic Boundaries, needs B, J|n - 1/2
    [Ex,Ey,Ez] = push_E(Bx,By,Bz,Ex,Ey,Ez,Jx,Jy,Jz,grid);
    [Ex,Ey,Ez,Bx,By,Bz,Uy,Uz,Jy,Jz] = BC(Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz,Ux,Uy,Uz,grid);

    %Advance B field (n - 1/2 -> n)
    [Bx,By,Bz] = push_B(Bx,By,Bz,Ex,Ey,Ez,grid);
    [Ex,Ey,Ez,Bx,By,Bz,Uy,Uz,Jy,Jz] = BC(Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz,Ux,Uy,Uz,grid);

elseif grid.solve_type_field == "Muscl"

    %Deposit the Current (time n - 1/2)
    [Jx,Jy,Jz] = J_deposition(N,Vx,Vy,Vz,grid);

    %Apply E-source via current n - 1/2, requires J|n - 1/2
    %Pushes E n - 1 to E n*
    [Ex,Ey,Ez] = muscl_field_source(Ex,Ey,Ez,Jx,Jy,Jz,grid); 

    %Apply E & B BC at n*, (n - 1)
    [Ex,Ey,Ez,Bx,By,Bz] = muscl_field_BC(Ex,Ey,Ez,Bx,By,Bz,grid); 

    %Advance E & B field (n - 1 -> n) using MUSCL:
    [Ex,Ey,Ez,Bx,By,Bz] = muscl_field_push(Ex,Ey,Ez,Bx,By,Bz,grid);

    %Apply E & B BC at (n - 1) - Corrects bad boundaries
    [Ex,Ey,Ez,Bx,By,Bz] = muscl_field_BC(Ex,Ey,Ez,Bx,By,Bz,grid); 


else
    if grid.iter == 2
        fprintf("No Field Solve selected!");
    end
end

%Moving Frame (at time n):
[N,Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz,Ux,Uy,Uz,grid] = moving_frame(N,Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz,Ux,Uy,Uz,grid);


end