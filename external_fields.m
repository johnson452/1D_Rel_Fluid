function [Bx,By,Bz,Ex,Ey,Ez] = external_fields(Bx,By,Bz,Ex,Ey,Ez,grid)

if grid.problem_name == "JE9"

%Add fields generated outside the system
Bx = Bx + grid.external_Bx;
By = By + grid.external_By;
Bz = Bz + grid.external_Bz;

Ex = Ex + grid.external_Ex;
Ey = Ey + grid.external_Ey;
Ez = Ez + grid.external_Ez;

end

end