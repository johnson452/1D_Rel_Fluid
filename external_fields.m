function [Bx_e,By_e,Bz_e,Ex_e,Ey_e,Ez_e] = external_fields(Bx,By,Bz,Ex,Ey,Ez,grid)

%Add fields generated outside the system
Bx_e = Bx + grid.external_Bx;
By_e = By + grid.external_By;
Bz_e = Bz + grid.external_Bz;

Ex_e = Ex + grid.external_Ex;
Ey_e = Ey + grid.external_Ey;
Ez_e = Ez + grid.external_Ez;

end