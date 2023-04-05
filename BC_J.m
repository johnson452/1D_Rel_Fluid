function [Uy,Uz,Vy,Vz] = BC_J(Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz,Ux,Uy,Uz,Vx,Vy,Vz,N,grid)

%Grab the full gridsize
Nx = grid.Nx;

if grid.BC_type == "Perioidic"

    %Fix the periodic components:
    Bx_temp = [Bx(1),Bx(1)];
    By_temp = [(By(Nx-1) + By(1))/2];
    Bz_temp = [(Bz(Nx-1) + Bz(1))/2];

    Ex_temp = [(Ex(Nx-1) + Ex(1))/2];
    Ey_temp = [Ey(1),Ey(1)];
    Ez_temp = [Ez(1),Ez(1)];

    Ux_temp = [(Ux(Nx-1) + Ux(1))/2];
    Uy_temp = [Uy(1),Uy(1)];
    Uz_temp = [Uz(1),Uz(1)];


    [Ux_temp,Uy_temp,Uz_temp,Vx_temp,Vy_temp,Vz_temp] = ...
        fluid_U(Bx_temp,By_temp,Bz_temp,Ex_temp,Ey_temp,Ez_temp,...
        Ux_temp,Uy_temp,Uz_temp,grid);

    %Correct the boundary conditions
    Uy(1) = Uy_temp(1);
    Uy(Nx) =  Uy_temp(1);
    Uz(1) =  Uz_temp(1);
    Uz(Nx) =  Uz_temp(1);

    Vy(1) =  Vy_temp(1);
    Vy(Nx) =  Vy_temp(1);
    Vz(1) =  Vz_temp(1);
    Vz(Nx) =  Vz_temp(1);
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

    Uy(N_max) = Uy(N_max-1);
    Uz(N_max) = Uz(N_max-1);
    Vy(N_max) = Vy(N_max-1);
    Vz(N_max) = Vz(N_max-1);

    Vy(N_max) = 3*(1/(N(N_max)*grid.e0))*J0*sin(2*pi*fd*t)*sin(0.5*pi*min(1,t/t0))^2;
    gamma = (1/sqrt(1 - ( Vx(N_max-1)*Vx(N_max-1) + Vy(N_max)*Vy(N_max) + Vz(N_max)*Vz(N_max)  )/(grid.c*grid.c))) ;
    Uy(N_max) = Vy(N_max)*gamma;
end

end