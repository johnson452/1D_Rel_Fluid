function [Ey,Ez,Uy,Uz,Jy,Jz] = BC(Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz,Ux,Uy,Uz,grid)

if grid.BC_cond == "Periodic"
[Bx,Ey,Ez,Uy,Uz,Jy,Jz] = periodic(Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz,Ux,Uy,Uz,grid);
end

if grid.BC_cond == "Non_Periodic"
[Ey,Ez,Uy,Uz,Jy,Jz] = non_periodic(Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz,Ux,Uy,Uz,grid);
end

end