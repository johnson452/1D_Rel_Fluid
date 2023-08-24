function N = BC_N(N,Vx,Vy,Vz,grid)

%Periodic BC for N
if grid.BC_cond == "Periodic"
    %handled by interpolation
elseif grid.BC_type == "Propagation into a plasma wave beach"
    % All other cases
    N(grid.Nx) = N(grid.Nx-1);
    N(1) = N(2); %(grid.wp(1).*grid.wp(1))*grid.eps_0*grid.m0/(grid.e0*grid.e0);
elseif grid.BC_type == "WFA"
    % Last point is a Ghost cell in WFA (matches WarpX)
    N(grid.Nx) = grid.N0; %N(grid.Nx-2);
    % Only inject density, leave the rest zero:
    N(grid.Nx) = density_func(grid.x1(grid.Nx));
    N(grid.Nx-1) = density_func(grid.x1(grid.Nx-1));
    N(grid.Nx-2) = density_func(grid.x1(grid.Nx-2));
else
    % All other cases
    N(grid.Nx) = N(grid.Nx-1);
    N(1) = N(2);
end

end