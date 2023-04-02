function [Ey,Ez,Uy,Uz,Jy,Jz] = non_periodic(Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz,Ux,Uy,Uz,grid)

% Setup from:
% https://ammar-hakim.org/sj/je/je9/je9-cyclotron-tunneling.html
if grid.BC_type == "Tunneling through an electron-cyclotron cutoff layer"
    N_max = max(size(Jy));

    Ey(1) = Ey(2);
    Ez(1) = Ez(2);

    Ey(N_max) = Ey(N_max-1);
    Ez(N_max) = Ez(N_max-1);

end

end