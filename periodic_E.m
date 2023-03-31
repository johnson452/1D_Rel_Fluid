%Boundary conditions E
function [Ex,Ey,Ez] = periodic_E(Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz,grid)
%Retrieve Constants:
dt = grid.dt;
dx = grid.dx;
c_sq_dt = dt*grid.c*grid.c;
mu_0_c_sq_dt = c_sq_dt*grid.mu_0;

%Update Ex, Ey, Ez: at the boundaries (see slides for notes/diagrams)
Ey(1) = Ey(1) - c_sq_dt*(Bz(1) - Bz(grid.Bz_i_end))/dx - mu_0_c_sq_dt*Jy(1);
Ey(grid.Ey_i_end+1) = Ey(1);

Ez(1) = Ez(1) + c_sq_dt*(By(1) - By(grid.By_i_end))/dx - mu_0_c_sq_dt*Jz(1);
Ez(grid.Ez_i_end+1) = Ez(1);

end