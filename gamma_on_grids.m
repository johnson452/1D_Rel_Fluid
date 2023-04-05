function [gamma_center,gamma_edge] = gamma_on_grids(Vx,Vy,Vz,grid)
Vx_temp = interp_center_to_edge(Vx);
Vy_temp = interp_edge_to_center(Vy);
Vz_temp = interp_edge_to_center(Vy);
V2_Edge = Vx.*Vx + Vy_temp.*Vy_temp + Vz_temp.*Vz_temp;
V2_Center = Vx_temp.*Vx_temp + Vy.*Vy + Vz.*Vz;
gamma_center = 1./sqrt(1-V2_Edge/(grid.c^2));
gamma_edge = 1./sqrt(1-V2_Center/(grid.c^2));
end