function grid = energy_momentum_diagnostic(Bx,By,Bz,Ex,Ey,Ez,Jx,Jy,Jz,Ux,Uy,Uz,N,str,grid)
% Diagnostics for the total electromagnetic energy-momentum

%Externally applied fields
Nx = grid.Nx;
[Bx,By,Bz,Ex,Ey,Ez] = external_fields(Bx,By,Bz,Ex,Ey,Ez,grid);

% J dot E - Ohmic Heating P = integral x (J dot E)
grid.JdotE(grid.iter) = sum(Jx.*Ex) + sum(Jy.*Ey) + sum(Jz.*Ez);

%Interp all values to cell-centered quantities?
Bx = interp_edge_to_center(Bx);
Ey = interp_edge_to_center(Ey);
Ez = interp_edge_to_center(Ez);
Uy = interp_edge_to_center(Uy);
Uz = interp_edge_to_center(Uz);
N = interp_edge_to_center(N);

%Total Energy Calcs ( Eliminate double dounting edge cell )
grid.Total_Energy_E_field(grid.iter) = (grid.eps_0/2) * ( sum(Ex.*Ex) +...
    sum(Ey.*Ey) +...
    sum(Ez.*Ez) );
grid.Total_Energy_B_field(grid.iter) = (1/(grid.mu_0*2)) * ( sum(Bx.*Bx) + sum(By.*By) + sum(Bz.*Bz) );
grid.Total_Energy_field(grid.iter) = grid.Total_Energy_E_field(grid.iter)  + grid.Total_Energy_B_field(grid.iter);
grid.Total_Energy_ptcls(grid.iter) = grid.m0*grid.c*grid.c*...
    sum(N.* ...
    sqrt(1 + (Ux.*Ux + Uy.*Uy + Uz.*Uz)/(grid.c*grid.c)) - 1);

% Total Momentum Calcs
%Ptcls
total_U = grid.m0*[sum(N.*Ux),sum(N.*Uy),sum(N.*Uz)];
magnitude_momentum = sqrt( total_U(1)*total_U(1) + total_U(2)*total_U(2) +...
    total_U(3)*total_U(3) );
grid.Total_Momentum_ptcls(:,grid.iter) = total_U;
grid.Total_Momentum_Magnitude_ptcls(grid.iter) = magnitude_momentum;

%Fields:
S = zeros(3,Nx-1);
for i = 1:Nx-1
    E = [Ex(i), Ey(i), Ez(i)];
    B = [Bx(i), By(i), Bz(i)];
    S(:,i) = (1/grid.mu_0)*cross(E, B);
end
total_U_fields = [sum(S(1,:)),sum(S(2,:)),sum(S(3,:))];
magnitude_momentum_fields = sqrt( total_U_fields(1)*total_U_fields(1) + ...
    total_U_fields(2)*total_U_fields(2) +...
    total_U_fields(3)*total_U_fields(3) );
grid.Total_Momentum_fields(1:3,grid.iter) = total_U_fields;
grid.Total_Momentum_Magnitude_fields(grid.iter) = magnitude_momentum_fields;


% Output (CML)
if (mod ( grid.iter, grid.Output_interval ) == 0)
fileID = fopen(grid.filename,'a');
fprintf(fileID,"\n*** (START) Energy-Momentum Diagnostic Ouptput ***\n");
fprintf(fileID,"Printed at iteration: %d\n",grid.iter);
fprintf(fileID," -  Called from function: %s -\n",str);
fprintf(fileID,"E-Field Energy: %e\n",grid.Total_Energy_E_field(grid.iter));
fprintf(fileID,"B-Field Energy: %e\n",grid.Total_Energy_B_field(grid.iter));
fprintf(fileID,"Total Field Energy: %e\n",grid.Total_Energy_field(grid.iter));
fprintf(fileID,"Total Particle Energy: %e\n",grid.Total_Energy_ptcls(grid.iter));

fprintf(fileID,"Total Momentum Ptcls: x: %e, y: %e, z: %e\n",...
    grid.Total_Momentum_ptcls(1,grid.iter),grid.Total_Momentum_ptcls(2,grid.iter), ...
    grid.Total_Momentum_ptcls(3,grid.iter));
fprintf(fileID,"Total Momentum Ptcls: %e\n",grid.Total_Momentum_Magnitude_ptcls(grid.iter));
fprintf(fileID,"Total Momentum Fields: x: %e, y: %e, z: %e\n",...
    grid.Total_Momentum_fields(1,grid.iter),grid.Total_Momentum_fields(2,grid.iter), ...
    grid.Total_Momentum_fields(3,grid.iter));
fprintf(fileID,"Total Momentum Fields: %e\n",grid.Total_Momentum_Magnitude_fields(grid.iter));

fprintf(fileID,"J dot E Heating Rate (Total): %e\n",grid.JdotE(grid.iter));
fprintf(fileID,"*** (END) Energy-Momentum Diagnostic Ouptput ***\n");
fclose(fileID);
end
end