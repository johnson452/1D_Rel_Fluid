function [Ex,Ey,Ez,Bx,By,Bz,Uy,Uz,Jy,Jz] = non_periodic(Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz,Ux,Uy,Uz,grid)

% Setup from:
% https://ammar-hakim.org/sj/je/je9/je9-cyclotron-tunneling.html
if (grid.BC_type == "Tunneling through an electron-cyclotron cutoff layer" ||...
    grid.BC_type == "Propagation into a plasma wave beach")
    N_max = max(size(Jy));

    Ex(1) = Ex(2);
    Ey(1) = Ey(2);
    Ez(1) = Ez(2);

    Ex(N_max-1) = Ex(N_max-2);
    Ey(N_max) = Ey(N_max-1);
    Ez(N_max) = Ez(N_max-1);

    Bx(1) = Bx(2);
    By(1) = By(2);
    Bz(1) = Bz(2);

    Bx(N_max) = Bx(N_max-1);
    By(N_max-1) = By(N_max-2);
    Bz(N_max-1) = Bz(N_max-2);

end

end