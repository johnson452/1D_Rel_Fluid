function grid = energy_momentum_diagnostic(Bx,By,Bz,Ex,Ey,Ez,Jx,Jy,Jz,Ux,Uy,Uz,N,str,grid)
% Diagnostics for the total electromagnetic energy-momentum

%Externally applied fields
Nx = grid.Nx;
[Bx,By,Bz,Ex,Ey,Ez] = external_fields(Bx,By,Bz,Ex,Ey,Ez,grid);

% Compute derivaitives before interpolation:
dx_NUx = grid.c*grid.c*grid.m0*(N(2:Nx).*Ux(2:Nx) - N(1:Nx-1).*Ux(1:Nx-1))/grid.dx;

%Interp all Values to be cell-centered
N = interp_edge_to_center_diag(N,grid);
Bx = interp_edge_to_center_diag(Bx,grid);
Ey = interp_edge_to_center_diag(Ey,grid);
Ez = interp_edge_to_center_diag(Ez,grid);
Jy = interp_edge_to_center_diag(Jy,grid);
Jz = interp_edge_to_center_diag(Jz,grid);
Ux = interp_edge_to_center_diag(Ux,grid);
Uy = interp_edge_to_center_diag(Uy,grid);
Uz = interp_edge_to_center_diag(Uz,grid);


%Correctly place E and B in time half a step (so they match fluid times)
if grid.iter >1
    Ex = (Ex + grid.Ex_old)/2;
    Ey = (Ey + grid.Ey_old)/2;
    Ez = (Ez + grid.Ez_old)/2;

    Bx = (Bx + grid.Bx_old)/2;
    By = (By + grid.By_old)/2;
    Bz = (Bz + grid.Bz_old)/2;
end

%Compute diff for derivative terms
sz_by = max(size(By));
sz_bz = max(size(Bz));
dx_By = (By(2:sz_by) - By(1:sz_by-1))/grid.dx;
dx_Bz = (Bz(2:sz_bz) - Bz(1:sz_bz-1))/grid.dx;

% Compute the extra energy terms of Maxwell:
field_energy_term_Ey_dx_Bz = -(1/grid.mu_0)*interp_edge_to_center_diag(Ey,grid).*dx_Bz;
field_energy_term_Ez_dx_By =  (1/grid.mu_0)*interp_edge_to_center_diag(Ez,grid).*dx_By;
x3 = interp_edge_to_center_diag(grid.x2,grid);


% Replace J with Fluid quantities (eliminates laser profile)
gamma = sqrt(1.0 + (Ux.*Ux + Uy.*Uy + Uz.*Uz)/(grid.c*grid.c));
Jx_fiuld = grid.e0*N.*Ux./gamma;
Jy_fiuld = grid.e0*N.*Uy./gamma;
Jz_fiuld = grid.e0*N.*Uz./gamma;
%Jx = Jx_fiuld;
%Jy = Jy_fiuld;
%Jz = Jz_fiuld;

% J dot E - Ohmic Heating P = integral x (J dot E)
grid.JdotE(grid.iter) = grid.dx*(sum(Jx.*Ex) + sum(Jy.*Ey) + sum(Jz.*Ez));
grid.JdotE_x = (Jx.*Ex + Jy.*Ey + Jz.*Ez);


%Total Energy Calcs ( Eliminate double dounting edge cell )
grid.Total_Energy_E_field(grid.iter) = grid.dx*(grid.eps_0/2) * ( sum(Ex.*Ex) + sum(Ey.*Ey) + sum(Ez.*Ez) );
grid.Total_Energy_B_field(grid.iter) = grid.dx*(1/(grid.mu_0*2)) * ( sum(Bx.*Bx) + sum(By.*By) + sum(Bz.*Bz) );
grid.Total_Energy_field(grid.iter) = grid.Total_Energy_E_field(grid.iter)  + grid.Total_Energy_B_field(grid.iter);

% Spatial grid of energy:
grid.Energy_E_field_x = (grid.eps_0/2) * ( Ex.*Ex + Ey.*Ey + Ez.*Ez );
grid.Energy_B_field_x = (1/(grid.mu_0*2)) * ( Bx.*Bx + By.*By + Bz.*Bz);
grid.Energy_fluid_x =  grid.m0*grid.c*grid.c*...
    N.*(sqrt(1 + (Ux.*Ux + Uy.*Uy + Uz.*Uz)/(grid.c*grid.c)));

%Relativistic energy is difficult to calculate when U << c
grid.Total_Energy_ptcls(grid.iter) = 0;
for i = 1:Nx-1
    if (Ux(i).*Ux(i) + Uy(i).*Uy(i) + Uz(i).*Uz(i))/(grid.c*grid.c) < 1e-12
        grid.Total_Energy_ptcls(grid.iter) = grid.Total_Energy_ptcls(grid.iter) +...
            grid.dx*grid.m0*...
            (N(i).* (0.5*(Ux(i).*Ux(i) + Uy(i).*Uy(i) + Uz(i).*Uz(i)))) + ...
        grid.dx*N(i)*grid.m0*grid.c*grid.c;
    else
        grid.Total_Energy_ptcls(grid.iter) = grid.Total_Energy_ptcls(grid.iter) +...
            grid.dx*grid.m0*grid.c*grid.c*...
            (N(i).* ...
            (sqrt(1 + (Ux(i).*Ux(i) + Uy(i).*Uy(i) + Uz(i).*Uz(i))/(grid.c*grid.c))));
    end
end

t = linspace(0,grid.time,grid.iter);


% Compute the discrete differences of each energy component
dE_field = diff(grid.Total_Energy_E_field(1:grid.iter))/grid.dt;
dB_field = diff(grid.Total_Energy_B_field(1:grid.iter))/grid.dt;
dTotal_field = diff(grid.Total_Energy_field(1:grid.iter))/grid.dt;
dPtcls = diff(grid.Total_Energy_ptcls(1:grid.iter))/grid.dt;

% Adjust time array for difference plots (using midpoints)
t_diff = t(1:end-1) + grid.dt/2;

if grid.iter == 1 || (mod ( grid.iter, grid.Output_interval ) == 0 || grid.iter == grid.NT)

    % Define a unique tag for your figure
    energyPlotTag = 'EnergyMomentumPlot';

    % Find the figure with this tag. If it doesn't exist, create a new figure.
    figHandle = findobj('Type', 'figure', 'Tag', energyPlotTag);
    if isempty(figHandle)
        figHandle = figure('units','normalized','outerposition',[0.0 0.5 0.5 0.55]); % Create a new figure
        set(figHandle, 'Tag', energyPlotTag); % Assign the unique tag to the figure
    else
        figure(figHandle); % Bring the existing figure to the front
        clf; % Clear the figure for new plots
    end

    subplot(2,3,1)
    hold on; % Hold on to plot multiple lines

    % Plot each energy component
    plot(t, grid.Total_Energy_E_field(1:grid.iter), 'b', 'LineWidth', 2);
    plot(t, grid.Total_Energy_B_field(1:grid.iter), 'r', 'LineWidth', 2);
    plot(t, grid.Total_Energy_field(1:grid.iter), 'black', 'LineWidth', 2);
    plot(t, grid.Total_Energy_ptcls(1:grid.iter), "color",[1, 0.843, 0], 'LineWidth', 2);

    % Adding labels and legend
    xlabel('Time');
    ylabel('Energy');
    title('Energy vs Time');
    legend('Total Energy E-field', 'Total Energy B-field', 'Total Field Energy', 'Total Energy Particles',"Location","best");

    hold off; % Release the plot

    if grid.iter >1
        subplot(2,3,2)
        hold on;

        % Plot the discrete differences
        plot(t_diff, dE_field, 'b', 'LineWidth', 2);
        plot(t_diff, dB_field, 'r', 'LineWidth', 2);
        plot(t_diff, dTotal_field, 'black', 'LineWidth', 2);
        plot(t_diff, dPtcls,"color",[1, 0.843, 0], 'LineWidth', 2);

        % Adding labels and legend for difference plots
        xlabel('Time');
        ylabel('Change in Energy');
        title('d/dt (Energy) v time');
        legend('Change in E-field Energy', 'Change in B-field Energy', 'Change in Total Field Energy', 'Change in Particle Energy', "Location", "best");

        hold off;

        subplot(2,3,3)
        hold on;

        % Plot the discrete differences
        plot(t_diff, dE_field, 'b', 'LineWidth', 2);
        plot(t_diff, dB_field, 'r', 'LineWidth', 2);
        plot(t_diff, dTotal_field, 'black', 'LineWidth', 2);
        plot(t_diff, dPtcls,"color",[1, 0.843, 0], 'LineWidth', 2);

        % Adding labels and legend for difference plots
        xlabel('Time');
        ylabel('Change in Energy');
        title('d/dt (Energy) v time (zoomed in)');
        legend('Change in E-field Energy', 'Change in B-field Energy', 'Change in Total Field Energy', 'Change in Particle Energy', "Location", "best");

        dL = max(dPtcls)-min(dPtcls) +1e-13;
        ylim([min(dPtcls)-dL,max(dPtcls)+dL])

        hold off;


        subplot(2,3,4)
        hold on;

        % Plot the discrete differences
        plot(t_diff, dE_field, 'b', 'LineWidth', 2);
        plot(t_diff, dB_field, 'r', 'LineWidth', 2);
        plot(t_diff, dTotal_field, 'black', 'LineWidth', 2);
        plot(t_diff, dPtcls,"color",[0.5, 0, 0.5], 'LineWidth', 2);
        plot(t, grid.JdotE(1:grid.iter),"color",[1, 0.843, 0], 'LineWidth', 2);
        %fprintf("max(dE/dt): %1.3e, max(dU/dt): %1.3e, max(JdotE): %1.3e\n", ...
        %max(abs(dTotal_field)),max(abs(dPtcls)),max(abs(grid.JdotE(1:grid.iter))));

        % Adding labels and legend for difference plots
        xlabel('Time');
        ylabel('Change in Energy');
        title('d/dt (Energy) v time (with J \cdot E)');
        legend("E Field","B Field",'Total Field ', 'Particle','J \cdot E', "Location", "best");


        hold off;


        subplot(2,3,5)
        hold on;

        % Plot the discrete differences
        dx_NUx_half  = (dx_NUx + grid.dx_NUx_old)/2;
        JdotE_x_half  = (grid.JdotE_x + grid.JdotE_x_old)/2;

        Total_energy_eq = ((grid.Energy_fluid_x) - (grid.Energy_fluid_x_old))/grid.dt  + dx_NUx_half - JdotE_x_half;

        % Plotting
        plot(grid.x2, (grid.Energy_E_field_x - grid.Energy_E_field_x_old)/grid.dt, 'b', 'LineWidth', 2);
        plot(grid.x2, Total_energy_eq,"color","black" ,'LineWidth', 2);
        plot(grid.x2, ((grid.Energy_fluid_x) - (grid.Energy_fluid_x_old))/grid.dt,"color","r", 'LineWidth', 2);
        plot(grid.x2, JdotE_x_half,"color",[1, 0.843, 0], 'LineWidth', 2);
        plot(grid.x2, dx_NUx_half,"color",[0.5, 0, 0.5], 'LineWidth', 2);


        % Adding labels and legend for difference plots
        xlabel('x [m]');
        ylabel('d/dt (Energy)  [J/s]');
        title('d/dt (Energy) vs x [Fluids]');
        legend("\epsilon_0\partial_t(E^2)/2","Fluid Energy Eqn.",'Fluid Energy','J \cdot E', '\nabla(NUx)',"Location", "best");

        dL = max( grid.JdotE_x)-min( grid.JdotE_x) +1e-13;
        ylim([min(grid.JdotE_x)-dL,max(grid.JdotE_x)+dL])


        hold off;


        subplot(2,3,6)
        hold on;

        % Compute half-timestep terms:
        field_energy_term_Ez_dx_By_half = (grid.field_energy_term_Ez_dx_By_old + field_energy_term_Ez_dx_By)/2;
        field_energy_term_Ey_dx_Bz_half = (grid.field_energy_term_Ey_dx_Bz_old + field_energy_term_Ey_dx_Bz)/2;

        % Plotting
        plot(grid.x2, -(grid.Energy_E_field_x - grid.Energy_E_field_x_old)/grid.dt, 'b', 'LineWidth', 2);
        %plot(x3, field_energy_term_Ez_dx_By_half,"color","black" ,'LineWidth', 2);
        plot(x3, field_energy_term_Ey_dx_Bz_half,"color",[0.5, 0, 0.5], 'LineWidth', 2);
        plot(grid.x2, JdotE_x_half,"color",[1, 0.843, 0], 'LineWidth', 2);
        plot(x3, -interp_edge_to_center(grid.Energy_E_field_x - grid.Energy_E_field_x_old,grid)/grid.dt + field_energy_term_Ey_dx_Bz_half,'--','color', 'black', 'LineWidth', 2);


        % Adding labels and legend for difference plots
        xlabel('x [m]');
        ylabel('d/dt (Energy)  [J/s]');
        title('d/dt (Energy) vs x [Fields]');
        %legend("d/dt (E-field energy) -E_y(\partial_x B_z)/\epsilon_0",'J \cdot E',"Location", "best");
        legend("\epsilon_0\partial_t(E^2)/2",'-E_y(\partial_x B_z)/\epsilon_0','J \cdot E',"\epsilon_0\partial_t(E^2)/2-E_y(\partial_x B_z)/\epsilon_0","Location", "best");

        dL = max( grid.JdotE_x)-min( grid.JdotE_x) +1e-13;
        ylim([min(grid.JdotE_x)-dL,max(grid.JdotE_x)+dL])


        hold off;




        if grid.BC_type == "fluid_only_diagnostic"
            clf();
            subplot(1,2,1)
            hold on;


            % Plot the discrete differences
            dx_NUx_half  = (dx_NUx + grid.dx_NUx_old)/2;

            Total_energy_eq = ((grid.Energy_fluid_x) - (grid.Energy_fluid_x_old))/grid.dt  + dx_NUx_half - JdotE_x_half;

            % Plotting
            plot(grid.x2, Total_energy_eq,"color","black" ,'LineWidth', 2);
            plot(grid.x2, ((grid.Energy_fluid_x) - (grid.Energy_fluid_x_old))/grid.dt,"color","r", 'LineWidth', 2);
            plot(grid.x2, dx_NUx_half,"color",[0.5, 0, 0.5], 'LineWidth', 2);


            % Adding labels and legend for difference plots
            xlabel('x [m]');
            ylabel('\partial_t(Energy)  [J/s]');
            title('\partial_t(Energy) vs x [Fluids]');
            legend("Fluid Energy Eqn.",'Change in Fluid Energy', '\nabla(NUx)',"Location", "best");

            hold off;

            subplot(1,2,2)
            hold on;

            % Plotting
            plot(grid.x2, Total_energy_eq,"color","black" ,'LineWidth', 2);


            % Adding labels and legend for difference plots
            xlabel('x [m]');
            ylabel('\partial_t(Energy)  [J/s]');
            title('\partial_t(Energy) vs x [Fluids]');
            legend("Fluid Energy Eqn.","Location", "best");

            hold off;
        end


    end
end


%Save old timesteps:
grid.Energy_E_field_x_old = grid.Energy_E_field_x;
grid.Energy_B_field_x_old = grid.Energy_B_field_x;
grid.Energy_fluid_x_old = grid.Energy_fluid_x;
grid.JdotE_x_old = grid.JdotE_x;
grid.dx_NUx_old = dx_NUx;
grid.Ex_old = Ex;
grid.Ey_old = Ey;
grid.Ez_old = Ez;
grid.Bx_old = Bx;
grid.By_old = By;
grid.Bz_old = Bz;
grid.field_energy_term_Ey_dx_Bz_old = field_energy_term_Ey_dx_Bz;
grid.field_energy_term_Ez_dx_By_old = field_energy_term_Ez_dx_By;


% % Total Momentum Calcs
% %Ptcls
% total_U = grid.dx*grid.m0*[sum(N.*Ux),sum(N.*Uy),sum(N.*Uz)];
% magnitude_momentum = sqrt( total_U(1)*total_U(1) + total_U(2)*total_U(2) +...
%     total_U(3)*total_U(3) );
% grid.Total_Momentum_ptcls(:,grid.iter) = total_U;
% grid.Total_Momentum_Magnitude_ptcls(grid.iter) = magnitude_momentum;
%
% %Fields:
% S = zeros(3,Nx-1);
% for i = 1:Nx-1
%     E = [Ex(i), Ey(i), Ez(i)];
%     B = [Bx(i), By(i), Bz(i)];
%     S(:,i) = (1/grid.mu_0)*cross(E, B);
% end
% total_U_fields = grid.dx*grid.mu_0*grid.eps_0*[sum(S(1,:)),sum(S(2,:)),sum(S(3,:))];
% magnitude_momentum_fields = sqrt( total_U_fields(1)*total_U_fields(1) + ...
%     total_U_fields(2)*total_U_fields(2) +...
%     total_U_fields(3)*total_U_fields(3) );
% grid.Total_Momentum_fields(1:3,grid.iter) = total_U_fields;
% grid.Total_Momentum_Magnitude_fields(grid.iter) = magnitude_momentum_fields;
%
%
% % Output (CML)
% if (mod ( grid.iter, grid.Output_interval ) == 0)
%     fileID = fopen(grid.filename,'a');
%     fprintf(fileID,"\n*** (START) Energy-Momentum Diagnostic Output ***\n");
%     fprintf(fileID,"Printed at iteration: %d\n",grid.iter);
%     fprintf(fileID," -  Called from function: %s -\n",str);
%     fprintf(fileID,"E-Field Energy: %e\n",grid.Total_Energy_E_field(grid.iter));
%     fprintf(fileID,"B-Field Energy: %e\n",grid.Total_Energy_B_field(grid.iter));
%     fprintf(fileID,"Total Field Energy: %e\n",grid.Total_Energy_field(grid.iter));
%     fprintf(fileID,"Total Particle Energy: %e\n",grid.Total_Energy_ptcls(grid.iter));
%
%     fprintf(fileID,"Total Momentum Ptcls: x: %e, y: %e, z: %e\n",...
%         grid.Total_Momentum_ptcls(1,grid.iter),grid.Total_Momentum_ptcls(2,grid.iter), ...
%         grid.Total_Momentum_ptcls(3,grid.iter));
%     fprintf(fileID,"Total Momentum Ptcls: %e\n",grid.Total_Momentum_Magnitude_ptcls(grid.iter));
%     fprintf(fileID,"Total Momentum Fields: x: %e, y: %e, z: %e\n",...
%         grid.Total_Momentum_fields(1,grid.iter),grid.Total_Momentum_fields(2,grid.iter), ...
%         grid.Total_Momentum_fields(3,grid.iter));
%     fprintf(fileID,"Total Momentum Fields: %e\n",grid.Total_Momentum_Magnitude_fields(grid.iter));
%
%     fprintf(fileID,"J dot E Heating Rate (Total): %e\n",grid.JdotE(grid.iter));
%
%     % Print the differences in field/energy
%     if grid.iter > 1
%         %dE, dB, dU^2
%         fprintf(fileID,"dE-Field Energy: %e\n",grid.Total_Energy_E_field(grid.iter)-grid.Total_Energy_E_field(grid.iter-1));
%         fprintf(fileID,"dB-Field Energy: %e\n",grid.Total_Energy_B_field(grid.iter)-grid.Total_Energy_B_field(grid.iter-1));
%         fprintf(fileID,"dTotal Field Energy: %e\n",grid.Total_Energy_field(grid.iter)-grid.Total_Energy_field(grid.iter-1));
%         fprintf(fileID,"dTotal Particle Energy: %e\n",grid.Total_Energy_ptcls(grid.iter)-grid.Total_Energy_ptcls(grid.iter-1));
%     end
%     fprintf(fileID,"*** (END) Energy-Momentum Diagnostic Output ***\n");
%     fclose(fileID);
% end
end