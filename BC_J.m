function [Ux,Uy,Uz,Vx,Vy,Vz] = BC_J(Ux,Uy,Uz,Vx,Vy,Vz,N,grid)

if grid.BC_type == "Periodic"
    %Handled by Interpolations

end

if grid.BC_type == "WFA"
    Nx = grid.Nx;

    %Copy Boundaries
    Ux(1) = Ux(3); %%ALL NODAL
    Uy(1) = Uy(3);
    Uz(1) = Uz(3);
    Vx(1) = Vx(3);
    Vy(1) = Vy(3);
    Vz(1) = Vz(3);

    Ux(Nx) = Ux(Nx-2);
    Uy(Nx) = Uy(Nx-2);
    Uz(Nx) = Uz(Nx-2);
    Vx(Nx) = Vx(Nx-2);
    Vy(Nx) = Vy(Nx-2);
    Vz(Nx) = Vz(Nx-2);

end

if grid.BC_type == "Tunneling through an electron-cyclotron cutoff layer"
    Nx = grid.Nx;
    J0 = grid.J0;
    fd = grid.fd;
    t = grid.time;
    t0 = grid.t0;

    %Current at the boundary
    Vy(Nx-1) = 3*(1/(N(Nx-1)*grid.e0))*J0*sin(2*pi*fd*t)*sin(0.5*pi*min(1,t/t0))^2;
    gamma = (1/sqrt(1 - ( Vx(Nx-1)*Vx(Nx-1) + Vy(Nx-1)*Vy(Nx-1) + Vz(Nx-1)*Vz(Nx-1)  )/(grid.c*grid.c))) ;
    Uy(Nx-1) = Vy(Nx-1)*gamma;

    %Copy Boundaries
    Ux(1) = Ux(2);
    Uy(1) = Uy(2);
    Uz(1) = Uz(2);
    Vx(1) = Vx(2);
    Vy(1) = Vy(2);
    Vz(1) = Vz(2);

    Ux(Nx) = Ux(Nx-1);
    Uy(Nx) = Uy(Nx-1);
    Uz(Nx) = Uz(Nx-1);
    Vx(Nx) = Vx(Nx-1);
    Vy(Nx) = Vy(Nx-1);
    Vz(Nx) = Vz(Nx-1);


end


if grid.BC_type == "Propagation into a plasma wave beach"
    Nx = grid.Nx;
    J0 = grid.J0;
    omega = grid.omega;
    t = grid.time;

    %Current at the boundary
    Vx_temp = interp_center_to_edge(Vx,grid);
    Vy(Nx-1) = (1/1e10)*(1/(N(Nx-1)*grid.e0))*J0*sin(omega*t);
    gamma = (1/sqrt(1 - ( Vx_temp(Nx-1)*Vx_temp(Nx-1) + Vy(Nx-1)*Vy(Nx-1) + Vz(Nx-1)*Vz(Nx-1)  )/(grid.c*grid.c))) ;
    Uy(Nx-1) = Vy(Nx-1)*gamma;
    if abs(Vy(Nx-1)/grid.c) > 1
    fprintf("Vy(Nx-1)/c < 1 --> (%f) < 1\n",Vy(Nx-1)/grid.c);
    end

    %Copy Boundaries
    Ux(1) = Ux(2);
    Uy(1) = Uy(2);
    Uz(1) = Uz(2);
    Vx(1) = Vx(2);
    Vy(1) = Vy(2);
    Vz(1) = Vz(2);

    Ux(Nx) = Ux(Nx-1);
    Uy(Nx) = Uy(Nx-1);
    Uz(Nx) = Uz(Nx-1);
    Vx(Nx) = Vx(Nx-1);
    Vy(Nx) = Vy(Nx-1);
    Vz(Nx) = Vz(Nx-1);

end

end