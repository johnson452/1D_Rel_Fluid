function [Ex,Ey,Ez] = muscl_field_source(Ex,Ey,Ez,Jx,Jy,Jz,grid)

% MUSCL field source:
option = "CD";

if option == "CD"
    %Retrieve Constants:
    dt = grid.dt;
    c_sq_dt = dt*grid.c*grid.c;
    mu_0_c_sq_dt = c_sq_dt*grid.mu_0;

    %Update Ex, Ey, Ez
    Ex = Ex - mu_0_c_sq_dt*Jx;
    Ey = Ey - mu_0_c_sq_dt*Jy;
    Ez = Ez - mu_0_c_sq_dt*Jz;
end

end