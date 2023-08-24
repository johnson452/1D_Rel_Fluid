%Diagnostics
function grid = diagnostics(Bx,By,Bz,Ex,Ey,Ez,Jx,Jy,Jz,Ux,Uy,Uz,N,grid)

%Total Energy Calcs
if grid.solve_type_field == "FDTD"
    str = "diagnostics.m";
    %grid = energy_momentum_diagnostic(Bx,By,Bz,Ex,Ey,Ez,Jx,Jy,Jz,Ux,Uy,Uz,N,str,grid);
end


% Run only at select iterations:
if (mod ( grid.iter, grid.Output_interval ) == 0 || grid.iter == grid.NT)

    % Clear the figure
    clf()

    %Grab inital size:

    if grid.solve_type_field == "FDTD"
        %Plot the vectorfield for Bx, By, Bz:
        x_bx = grid.x1;
        x_by = grid.x2;
        x_bz = grid.x2;

        %Plot the vectorfield for Ex, Ey, Ez:
        x_ex = grid.x2;
        x_ey = grid.x1;
        x_ez = grid.x1;
    else
        %Plot the vectorfield for Bx, By, Bz:
        x_bx = grid.x1;
        x_by = grid.x1;
        x_bz = grid.x1;

        %Plot the vectorfield for Ex, Ey, Ez:
        x_ex = grid.x1;
        x_ey = grid.x1;
        x_ez = grid.x1;
    end


    %Time spacing:
    t_vec = linspace(0,grid.time,grid.iter);


    % Add any external Fields:
    [Bx,By,Bz,Ex,Ey,Ez] = external_fields(Bx,By,Bz,Ex,Ey,Ez,grid);

    %FILL IN: 1D plotting routines
    if grid.BC_type == "WFA"

        c = grid.c;

        %E-field plot
        subplot(3,3,1)
        plot(x_ex,Ex)
        title("Ex-fields")
        xlabel("x [m]")
        ylabel("Ex [V/m]")

        subplot(3,3,2)
        plot(x_ey,Ey)
        title("Ey-fields")
        xlabel("x [m]")
        ylabel("Ey [V/m]")

        %Current Density
        subplot(3,3,3)
        plot(x_ex,Jx)
        title("Current Density x")
        xlabel("x [m]")
        ylabel("Jx [C/sm^2]")

        subplot(3,3,4)
        plot(x_ey,Jy)
        title("Current Density y")
        xlabel("x [m]")
        ylabel("Jy [C/sm^2]")

        % Number density evolution:
        subplot(3,3,5)
        plot(grid.x1,N*grid.e0)
        title("Charge Density [C/m^3]")
        xlabel("x [m]")
        ylabel("Charge Density [C/m^3]")

        % Compute Vx, Vy, Vz:
        gamma = sqrt(1.0 + (Ux.*Ux + Uy.*Uy + Uz.*Uz)/(c^2));
        Vx = Ux./gamma;
        Vy = Uy./gamma;

        % Number density evolution:
        subplot(3,3,6)
        plot(grid.x1,Vx/c)
        title("Vx [v/c]")
        xlabel("x [m]")
        ylabel("Vx [v/c]")

        % Number density evolution:
        subplot(3,3,7)
        plot(grid.x1,Vy/c)
        title("Vy [v/c]")
        xlabel("x [m]")
        ylabel("Vy [v/c]")

        % Number density evolution:
        subplot(3,3,8)
        plot(grid.x1,Ux/c)
        title("Ux [u/c]")
        xlabel("x [m]")
        ylabel("Ux [u/c]")

        % Number density evolution:
        subplot(3,3,9)
        plot(grid.x1,Uy/c)
        title("Uy [u/c]")
        xlabel("x [m]")
        ylabel("Uy [u/c]")

    else

        %E-field plot
        subplot(2,3,1)
        plot(x_ex,Ex)
        hold on
        plot(x_ey,Ey)
        hold on
        plot(x_ez,Ez)
        title("E-fields")
        xlabel("x [m]")
        ylabel("E [V/m]")
        legend("Ex","Ey","Ez")

        %J-field plot
        subplot(2,3,2)
        plot(x_bx,Bx)
        hold on
        plot(x_by,By)
        hold on
        plot(x_bz,Bz)
        title("B-fields")
        xlabel("x [m]")
        ylabel("B [V/s]")
        legend("Bx","By","Bz")

        %B-field plot
        subplot(2,3,3)
        plot(x_ex,Jx)
        hold on
        plot(x_ey,Jy)
        hold on
        plot(x_ez,Jz)
        title("Current Density")
        xlabel("x [m]")
        ylabel("J [C/sm^2]")
        legend("Jx","Jy","Jz")

        % Total Energy in the Fields
        subplot(2,3,4)
        plot(t_vec,grid.Total_Energy_B_field(1:grid.iter))
        hold on
        plot(t_vec,grid.Total_Energy_E_field(1:grid.iter))
        hold on
        plot(t_vec,grid.Total_Energy_field(1:grid.iter))
        title("Field Energy")
        xlabel("t [s]")
        ylabel("Energy [J]")
        legend("B-Field","E-Field","Total Field")

        % Total ptcl/field energy partition
        subplot(2,3,5)
        plot(t_vec,grid.Total_Energy_field(1:grid.iter))
        hold on
        plot(t_vec,grid.Total_Energy_ptcls(1:grid.iter))
        hold on
        plot(t_vec,grid.Total_Energy_field(1:grid.iter)+grid.Total_Energy_ptcls(1:grid.iter))
        title("Ptcl and Field Energy")
        xlabel("t [s]")
        ylabel("Energy [J]")
        legend("Fields", "Particles","Total (Ptcls+Fields)")

        % Number density evolution:
        subplot(2,3,6)
        plot(grid.x1,N*grid.e0)
        title("Charge Density [C/m^3]")
        xlabel("x [m]")
        ylabel("Charge Density [C/m^3]")
    end

    %Print that it runs:
    fprintf("Output for: iteration %d\n",grid.iter);

    %Pause and then clear figure
    pause(0.01)
end

if grid.BC_type == "Propagation into a plasma wave beach"
    %Build extra plot
    output_interval = floor(linspace(0,grid.NT,grid.contour_size)) + 1;
    if ((grid.iter == grid.NT || max(output_interval == grid.iter)) ...
            && (grid.BC_type == "Propagation into a plasma wave beach"))

        % Ey vs time, x profile
        grid.Ey_t_x(grid.temp_iter,:) = Ey;
        grid.temp_iter = grid.temp_iter + 1;

        % Plot only at the last time:
        if grid.iter == grid.NT
            figure(2)
            time = linspace(0,grid.time,grid.contour_size);
            x = grid.x1;
            [~,h] = contourf(  x, time, grid.Ey_t_x, 50);
            ylabel("t")
            xlabel("x")
            title("Ey")
            colorbar()
            set(h,'LineColor','none')

            fprintf("Done!\n");

        end

    end
end



end