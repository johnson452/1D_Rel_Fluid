%Diagnostics
function grid = diagnostics(Bx,By,Bz,Ex,Ey,Ez,Jx,Jy,Jz,Ux,Uy,Uz,N,grid)

%Total Energy Calcs
str = "diagnostics.m";
grid = energy_momentum_diagnostic(Bx,By,Bz,Ex,Ey,Ez,Jx,Jy,Jz,Ux,Uy,Uz,N,str,grid);


% Run only at select iterations:
if (mod ( grid.iter, grid.Output_interval ) == 0)

    % Clear the figure
    clf()

    %Grab inital size:

    %Plot the vectorfield for Bx, By, Bz:
    x_bx = grid.x1;
    x_by = grid.x2;
    x_bz = grid.x2;

    %Plot the vectorfield for Ex, Ey, Ez:
    x_ex = grid.x2;
    x_ey = grid.x1;
    x_ez = grid.x1;

    %Time spacing:
    t_vec = linspace(0,grid.time,grid.iter);


    % Add any external Fields:
    [Bx,By,Bz,Ex,Ey,Ez] = external_fields(Bx,By,Bz,Ex,Ey,Ez,grid);

    %FILL IN: 1D plotting routines

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

    %Print that it runs:
    fprintf("Output for: iteration %d\n",grid.iter);

    %Pause and then clear figure
    pause(0.01)
end
end