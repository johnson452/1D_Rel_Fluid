function [N,Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz,Ux,Uy,Uz,grid] = moving_frame(N,Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz,Ux,Uy,Uz,grid)
%ASSUMES z+ shift at the speed of light!

%If the moving frame is active
if grid.moving_frame == 1

    % If this time to add in a new frame
    if (grid.c*grid.time - grid.xmax + grid.xmax0)/grid.dx > 1.0

        %Update grid xmax/xmin
        grid.xmin = grid.xmin + grid.dx;
        grid.xmax = grid.xmax + grid.dx;

        % Update grid x1, x2
        grid.x1 = linspace(grid.xmin,grid.xmax,grid.Nx);
        grid.x2 = linspace(grid.xmin+grid.dx/2,grid.xmax-grid.dx/2,grid.Nx-1);

        %Shift quantities
        N = shift_quant(N);
        Ex = shift_quant(Ex);
        Ey = shift_quant(Ey);
        Ez = shift_quant(Ez);
        Bx = shift_quant(Bx);
        By = shift_quant(By);
        Bz = shift_quant(Bz);
        Jx = shift_quant(Jx);
        Jy = shift_quant(Jy);
        Jz = shift_quant(Jz);
        Ux = shift_quant(Ux);
        Uy = shift_quant(Uy);
        Uz = shift_quant(Uz);

        %Inject the boundary values
        if grid.BC_type == "WFA"

            % Only inject density, leave the rest zero:
            N(grid.Nx) = density_func(grid.x1(grid.Nx));
            N(grid.Nx-1) = density_func(grid.x1(grid.Nx-1));
            N(grid.Nx-2) = density_func(grid.x1(grid.Nx-2));
        end

    end

end

end


function [x] = shift_quant(x)

%Compute the size of x
N = max(size(x));
NR = linspace(2,N,N-1);

x = x(NR);
x = [x,0]; %Uninitialized RHS

end