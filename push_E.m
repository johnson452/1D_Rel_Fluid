%%% Functions %%%
%Advance E field (n-1 -> n)
function [Ex,Ey,Ez] = push_E(Bx,By,Bz,Ex,Ey,Ez,Jx,Jy,Jz,grid)

%Retrieve Constants:
dt = grid.dt;
dx = grid.dx;
c_sq_dt = dt*grid.c*grid.c;
mu_0_c_sq_dt = c_sq_dt*grid.mu_0;

%Update Ex, Ey, Ez: %Subtract 1 to get to the B domain offset
for i = grid.Ex_i_0:grid.Ex_i_end
    Ex(i) = Ex(i) - mu_0_c_sq_dt*Jx(i);
end
for i = grid.Ey_i_0:grid.Ey_i_end

    Ey(i) = Ey(i) - c_sq_dt*(Bz(i+1 -1) - Bz(i -1))/dx - mu_0_c_sq_dt*Jy(i);
end
for i = grid.Ez_i_0:grid.Ez_i_end
    Ez(i) = Ez(i) + c_sq_dt*(By(i+1 -1) - By(i -1))/dx - mu_0_c_sq_dt*Jz(i);
end
end
