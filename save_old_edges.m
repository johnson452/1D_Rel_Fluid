function grid = save_old_edges(Uy,Uz,grid)

%Temp Saves
Nx = grid.Nx;
Uy_temp = [Uy(Nx-1),Uy(Nx),Uy(2)];
Uz_temp = [Uz(Nx-1),Uz(Nx),Uz(2)];

%Remember the prior BCs
grid.prior_bc_Uy = Uy_temp;
grid.prior_bc_Uz = Uz_temp;

end