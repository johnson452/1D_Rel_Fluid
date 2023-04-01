%Diagnostics
function grid = diagnostics(Bx,By,Bz,Ex,Ey,Ez,Jx,Jy,Jz,grid)

% Extra calculations:
% Eliminate double dounting edge cells

grid.Total_Energy_E_field(grid.iter) = (grid.eps_0/2) * ( sum(Ex.*Ex) +...
    sum(Ey(grid.Ey_i_0:grid.Ey_i_end).*Ey(grid.Ey_i_0:grid.Ey_i_end)) +...
    sum(Ez(grid.Ez_i_0:grid.Ez_i_end).*Ez(grid.Ez_i_0:grid.Ez_i_end)) );
grid.Total_Energy_B_field(grid.iter) = (1/(grid.mu_0*2)) * ( sum(Bx.*Bx) + sum(By.*By) + sum(Bz.*Bz) );
grid.Total_Energy_field(grid.iter) = grid.Total_Energy_E_field(grid.iter)  + grid.Total_Energy_B_field(grid.iter);


% Run only at select iterations:
if (mod ( grid.iter, 10) == 0)
    
    % Clear the figure
    clf()
    
    %Grab inital size:
    Nx = grid.Nx;
    
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
    title("J-fields")
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
    
    %Print that it runs:
    fprintf("Output for: iteration %d\n",grid.iter);
    
    %Pause and then clear figure
    pause(0.001)
end
end