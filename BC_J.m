function [Uy,Uz,Vy,Vz] = BC_J(Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz,Ux,Uy,Uz,Vx,Vy,Vz,N,grid)

%Grab the full gridsize
Nx = grid.Nx;

if grid.BC_type == "Periodic"
    %Handled by Interpolations

end

if grid.BC_type == "Tunneling through an electron-cyclotron cutoff layer"
    N_max = max(size(Jy));
    J0 = grid.J0;
    fd = grid.fd;
    t = grid.time;
    t0 = grid.t0;

    Uy(1) = Uy(2);
    Uz(1) = Uz(2);
    Vy(1) = Vy(2);
    Vz(1) = Vz(2);

    Uy(N_max) = 0; %Uy(N_max-1);
    Uz(N_max) = 0; %Uz(N_max-1);
    Vy(N_max) = 0; %Vy(N_max-1);
    Vz(N_max) = 0; %Vz(N_max-1);

    Vy(N_max) = 3*(1/(N(N_max)*grid.e0))*J0*sin(2*pi*fd*t)*sin(0.5*pi*min(1,t/t0))^2;
    gamma = (1/sqrt(1 - ( Vx(N_max-1)*Vx(N_max-1) + Vy(N_max)*Vy(N_max) + Vz(N_max)*Vz(N_max)  )/(grid.c*grid.c))) ;
    Uy(N_max) = Vy(N_max)*gamma;
end

end