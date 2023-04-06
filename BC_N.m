function N = BC_N(N,Vx,Vy,Vz,grid)

%Periodic BC for N
if grid.BC_cond == "Periodic"

    %^Grid
    Nx = grid.Nx;

    %Grab the velocity values for the edge calc
    Vx_temp = [Vx(Nx-1),Vx(1)];
    Vy_temp = [Vy(Nx-1),Vy(Nx),Vy(2)];
    Vz_temp = [Vz(Nx-1),Vz(Nx),Vz(2)];

    %update the boundaries of N
    [N_temp,~] = push_N(grid.N_old,Vx_temp,Vy_temp,Vz_temp,grid);
    N(1) = N_temp(2);
    N(Nx) = N_temp(2);
end

end