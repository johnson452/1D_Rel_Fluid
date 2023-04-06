function [N,grid] = push_N(N,Vx,Vy,Vz,grid)

%Scheme order
order = "FOE";

%Scheme
if order == "FOE"

    % Time FOE, space (CD - ish? For Vx yes, for N kinda)
    Nx = max(size(N));

    %Old cells
    if grid.Nx == max(size(N))
        grid.N_old = [N(Nx-2),N(Nx),N(2)];
    end

    %Calculate derivatives
    dxVx = (Vx(2:Nx-1) - Vx(1:Nx-2))/grid.dx;
    dxN = ( (N(3:Nx) + N(2:Nx-1))/2 - (N(2:Nx-1) + N(1:Nx-2))/2 )/grid.dx;

    % Update N
    Vx_interp = interp_center_to_edge(Vx);
    N(2:Nx-1) = N(2:Nx-1) - grid.dt*(N(2:Nx-1).*dxVx + Vx_interp(2:Nx-1).*dxN);

end


end