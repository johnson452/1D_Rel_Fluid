function [Ey,Ez,Uy,Uz,Jy,Jz] = non_periodic(Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz,Ux,Uy,Uz,grid)

% Setup from:
% https://ammar-hakim.org/sj/je/je9/je9-cyclotron-tunneling.html
if grid.BC_type == "Tunneling through an electron-cyclotron cutoff layer"
    N_max = max(size(Jy));
    J0 = grid.J0;
    fd = grid.fd;
    t = grid.time;
    t0 = grid.t0;
    Jy(N_max) = J0*sin(2*pi*fd*t)*sin(0.5*pi*min(1,t/t0))^2;

    Ey(1) = 0;
    Ez(1) = 0;
    Uy(1) = 0;
    Uz(1) = 0;
    Jy(1) = 0;
    Jz(1) = 0;

    Ey(N_max) = 0;
    Ez(N_max) = 0;
    Uy(N_max) = 0;
    Uz(N_max) = 0;
    Jz(N_max) = 0;
end

end