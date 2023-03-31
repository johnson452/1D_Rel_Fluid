%Intial Conditions
function [Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz,Ux,Uy,Uz] = IC(Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz,Ux,Uy,Uz,grid)


%Grab inital size:
Nx = grid.Nx;

%Intial Current
%Jx = zeros(Nx-1,1);
%Jy = zeros(Nx,1);
%Jz = zeros(Nx,1);


%Create a photon ( E only)
%for i = (grid.Ey_i_0-1):grid.Ey_i_end
%    Ey(i) = 1.0*sin(3*2*(pi)*(i-1)/(grid.Ey_i_end));
%end

end
