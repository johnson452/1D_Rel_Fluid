function N = BC_N(N,Vx,Vy,Vz,grid)

%Periodic BC for N
if grid.BC_cond == "Periodic"

else
    N(grid.Nx) = N(grid.Nx-1);
    N(1) = N(2);
end

end