%Advance B field NOTE: halfpush 
function [Bx,By,Bz] = push_B(Bx,By,Bz,Ex,Ey,Ez,grid)

%Retrieve Constants:
dt = grid.dt*0.5;
dx = grid.dx;

%Update Bx, By, Bz:
for i = grid.Bx_i_0:grid.Bx_i_end
    Bx(i) = Bx(i); % No Change in 1D!
end
for i = grid.By_i_0:grid.By_i_end
    By(i) = By(i) + dt*(Ez(i+1) - Ez(i))/dx;
end
for i = grid.Bz_i_0:grid.Bz_i_end
    Bz(i) = Bz(i) - dt*(Ey(i+1) - Ey(i))/dx;
end
end