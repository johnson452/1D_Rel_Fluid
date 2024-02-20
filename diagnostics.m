%Diagnostics
function grid_sim = diagnostics(Bx,By,Bz,Ex,Ey,Ez,Jx,Jy,Jz,Ux,Uy,Uz,N,grid_sim)

%Total Energy Calcs
%if grid_sim.solve_type_field == "FDTD" 
    str = "diagnostics.m";
    grid_sim = energy_momentum_diagnostic(Bx,By,Bz,Ex,Ey,Ez,Jx,Jy,Jz,Ux,Uy,Uz,N,str,grid_sim);
%end

%Save the state of the simulation
if (mod(grid_sim.iter,5000)==0)
    file_name = sprintf("wfa_save_state_%d.mat",grid_sim.iter);
    save(file_name);
end





% Run only at select iterations:
if grid_sim.iter == 1 || (mod ( grid_sim.iter, grid_sim.Output_interval ) == 0 || grid_sim.iter == grid_sim.NT)

    % Define a unique tag for your diagnostic figure
    diagnosticPlotTag = 'diagnostic';
    
    % Find the figure with this tag. If it doesn't exist, create a new figure.
    figHandle = findobj('Type', 'figure', 'Tag', diagnosticPlotTag);
    if isempty(figHandle)
        figHandle = figure('units','normalized','outerposition',[0 0 0.5 0.45]); % Create a new figure with specified shape
        set(figHandle, 'Tag', diagnosticPlotTag); % Assign the unique tag to the figure
    else
        figure(figHandle); % Bring the existing figure to the front
        clf; % Clear the figure for new plots
    end

    %Grab inital size:

    if grid_sim.solve_type_field == "FDTD"
        %Plot the vectorfield for Bx, By, Bz:
        x_bx = grid_sim.x1;
        x_by = grid_sim.x2;
        x_bz = grid_sim.x2;

        %Plot the vectorfield for Ex, Ey, Ez:
        x_ex = grid_sim.x2;
        x_ey = grid_sim.x1;
        x_ez = grid_sim.x1;
    else
        %Plot the vectorfield for Bx, By, Bz:
        x_bx = grid_sim.x1;
        x_by = grid_sim.x1;
        x_bz = grid_sim.x1;

        %Plot the vectorfield for Ex, Ey, Ez:
        x_ex = grid_sim.x1;
        x_ey = grid_sim.x1;
        x_ez = grid_sim.x1;
    end


    %Time spacing:
    t_vec = linspace(0,grid_sim.time,grid_sim.iter);


    % Add any external Fields:
    [Bx,By,Bz,Ex,Ey,Ez] = external_fields(Bx,By,Bz,Ex,Ey,Ez,grid_sim);

    %FILL IN: 1D plotting routines
    if grid_sim.BC_type == "WFA" 

        c = grid_sim.c;

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
        plot(grid_sim.x1,N*grid_sim.e0)
        title("Charge Density [C/m^3]")
        xlabel("x [m]")
        ylabel("Charge Density [C/m^3]")

        % Compute Vx, Vy, Vz:
        gamma = sqrt(1.0 + (Ux.*Ux + Uy.*Uy + Uz.*Uz)/(c^2));
        Vx = Ux./gamma;
        Vy = Uy./gamma;

        % Number density evolution:
        subplot(3,3,6)
        plot(grid_sim.x1,Vx/c)
        title("Vx [v/c]")
        xlabel("x [m]")
        ylabel("Vx [v/c]")

        % Number density evolution:
        subplot(3,3,7)
        plot(grid_sim.x1,Vy/c)
        title("Vy [v/c]")
        xlabel("x [m]")
        ylabel("Vy [v/c]")

        % Number density evolution:
        subplot(3,3,8)
        plot(grid_sim.x1,Ux/c)
        title("Ux [u/c]")
        xlabel("x [m]")
        ylabel("Ux [u/c]")

        % Number density evolution:
        subplot(3,3,9)
        %         plot(grid_sim.x1,Uy/c)
        %         title("Uy [u/c]")
        %         xlabel("x [m]")
        %         ylabel("Uy [u/c]")

        if grid_sim.BC_type == "WFA"
            rho = grid_sim.e0*N;
            x1 = grid_sim.x1;

            T = grid_sim.dx;
            Fs = 1/T;
            t = x1;
            L = length(t);
            X = rho;
            n = 2^nextpow2(L);
            Y = fft(X,n);
            f = 1./(Fs*(0:(n/2))/n);
            P = abs(Y/n).^2;

            loglog(f/grid_sim.dx,P(1:n/2+1)) 
            title("FFT(\rho)")
            xlabel("\lambda/dx")
            ylabel("|P(f)|^2")
        end

    elseif  grid_sim.BC_type == "fluid_only_diagnostic"
        
        % Number density evolution:
        subplot(2,2,1)
        plot(grid_sim.x1,N)
        title("N")
        xlabel("x [m]")
        ylabel("N [1/m^3]")
        if grid_sim.iter == 3000
            grid on
        end

        % Compute Vx, Vy, Vz:
        c = grid_sim.c;
        gamma = sqrt(1.0 + (Ux.*Ux + Uy.*Uy + Uz.*Uz)/(c^2));
        Vx = Ux./gamma;
        Vy = Uy./gamma;
        NUx = N.*Ux;
        NUy = N.*Uy;

        %Current Density
        subplot(2,2,2)
        plot(grid_sim.x1,NUx)
        title("NUx")
        xlabel("x [m]")
        ylabel("NUx [1/m^2s]")
        if grid_sim.iter == 3000
            grid on
        end

        subplot(2,2,3)
        plot(grid_sim.x1,NUy)
        title("NUy")
        xlabel("x [m]")
        ylabel("NUy [1/m^2s]")
        if grid_sim.iter == 3000
            grid on
        end

        subplot(2,2,4)
        plot(grid_sim.x1,Ux)
        title("Ux")
        xlabel("x [m]")
        ylabel("Ux [1/m^2s]")
        if grid_sim.iter == 3000
            grid on
        end
        ylim([-0.5,0.5]*1e9)


    else

        if grid_sim.solve_type_field ~= "None"

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
        plot(t_vec,grid_sim.Total_Energy_B_field(1:grid_sim.iter))
        hold on
        plot(t_vec,grid_sim.Total_Energy_E_field(1:grid_sim.iter))
        hold on
        plot(t_vec,grid_sim.Total_Energy_field(1:grid_sim.iter))
        title("Field Energy")
        xlabel("t [s]")
        ylabel("Energy [J]")
        legend("B-Field","E-Field","Total Field")

        % Total ptcl/field energy partition
        subplot(2,3,5)
        plot(t_vec,grid_sim.Total_Energy_field(1:grid_sim.iter))
        hold on
        plot(t_vec,grid_sim.Total_Energy_ptcls(1:grid_sim.iter))
        hold on
        plot(t_vec,grid_sim.Total_Energy_field(1:grid_sim.iter)+grid_sim.Total_Energy_ptcls(1:grid_sim.iter))
        title("Ptcl and Field Energy")
        xlabel("t [s]")
        ylabel("Energy [J]")
        legend("Fields", "Particles","Total (Ptcls+Fields)")

        % Number density evolution:
        subplot(2,3,6)
        plot(grid_sim.x1,N*grid_sim.e0)
        title("Charge Density [C/m^3]")
        xlabel("x [m]")
        ylabel("Charge Density [C/m^3]")
        end
    end

    %Print that it runs:
    fprintf("Output for: iteration %d\n",grid_sim.iter);

    %Pause and then clear figure
    pause(0.01)
end

if grid_sim.BC_type == "Propagation into a plasma wave beach"
    %Build extra plot
    output_interval = floor(linspace(0,grid_sim.NT,grid_sim.contour_size)) + 1;
    if ((grid_sim.iter == grid_sim.NT || max(output_interval == grid_sim.iter)) ...
            && (grid_sim.BC_type == "Propagation into a plasma wave beach"))

        % Ey vs time, x profile
        grid_sim.Ey_t_x(grid_sim.temp_iter,:) = Ey;
        grid_sim.temp_iter = grid_sim.temp_iter + 1;

        % Plot only at the last time:
        if grid_sim.iter == grid_sim.NT
            figure(2)
            time = linspace(0,grid_sim.time,grid_sim.contour_size);
            x = grid_sim.x1;
            [~,h] = contourf(  x, time, grid_sim.Ey_t_x, 50);
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