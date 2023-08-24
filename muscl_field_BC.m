function [Ex,Ey,Ez,Bx,By,Bz] = muscl_field_BC(Ex,Ey,Ez,Bx,By,Bz,grid)

if grid.solve_type_field == "Muscl"
    if (grid.BC_type == "WFA")
        N_max = max(size(Bz));

        % PEC:
        % E|| = 0, B_perp = 0

        Ex(1:2) = [0,Ex(3)];
        Ey(1:2) = [0,0];
        Ez(1:2) = [0,0];

        Ex(N_max-1:N_max) = [0,0];
        Ey(N_max-1:N_max) = [0,0];
        Ez(N_max-1:N_max) = [0,0];

        Bx(1:2) = [0,0];
        By(1:2) = [0,By(3)];
        Bz(1:2) = [0,Bz(3)];

        Bx(N_max-1:N_max) = [0,0];
        By(N_max-1:N_max) = [0,0];
        Bz(N_max-1:N_max) = [0,0];

    end
else
    fprintf("muscl_field_BC() - inappropriately called\n");
end

end