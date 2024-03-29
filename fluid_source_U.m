%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Name: Grant Johnson
%Date: 3/31/2022
%Fluid algroithm (ties into Yee_1D_fluid_u.m

%Notes:
%-1D
% Fluid Scheme: https://ammar-hakim.org/sj/hancock-muscl.html
% E + (v x B) : H&C
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Update the quanitites Ux, Uy, Uz (t -> t + dt)
function [Ux,Uy,Uz,Vx,Vy,Vz,grid] = fluid_source_U(Bx,By,Bz,Ex,Ey,Ez,Ux,Uy,Uz,grid)

% **** IMPORTANT ****

%External Fields if applied
[Bx,By,Bz,Ex,Ey,Ez] = external_fields(Bx,By,Bz,Ex,Ey,Ez,grid);
if grid.BC_type == "WFA" && grid.laser_envelope_model == "true"
    %Envelope model of the laser
    [Ex,Ey,Ez,Bx,By,Bz] = external_field_fun(Ex,Ey,Ez,Bx,By,Bz,grid);
end

%Interp all values to cell-centered quantities if FDTD
if grid.solve_type_field == "FDTD"
    Ex = interp_center_to_edge(Ex,grid);
    By = interp_center_to_edge(By,grid);
    Bz = interp_center_to_edge(Bz,grid);
end
% **** IMPORTANT ****

%%% H&C Algo %%%
%Grab the vector values (at original time)
h = grid.dt;
m0 = grid.m0;
qmdt = grid.e0*0.5*h./(m0);
c = grid.c;
E_0 = qmdt.*Ex;
E_1 = qmdt.*Ey;
E_2 = qmdt.*Ez;
B_0 = qmdt.*Bx;
B_1 = qmdt.*By;
B_2 = qmdt.*Bz;

u_0_minus = Ux + E_0;
u_1_minus = Uy + E_1;
u_2_minus = Uz + E_2;

u_star = u_0_minus.*(B_0/c) + u_1_minus.*(B_1/c) + u_2_minus.*(B_2/c); 
gamma_minus = sqrt(1 + (u_0_minus.*u_0_minus + u_1_minus.*u_1_minus + u_2_minus.*u_2_minus)/(c*c)); 
dot_tau_tau = (B_0.*B_0 + B_1.*B_1 + B_2.*B_2); 
sigma = gamma_minus.*gamma_minus - dot_tau_tau; 
gamma_new = sqrt(  0.5*(   sigma + sqrt( sigma.*sigma + 4*(dot_tau_tau + u_star.*u_star ) ) )  ); 

t_0 = B_0./gamma_new; 
t_1 = B_1./gamma_new; 
t_2 = B_2./gamma_new; 
s = 1./(1+(t_0.*t_0 + t_1.*t_1 + t_2.*t_2)); 

% FROM WarpX:
umt = u_0_minus.*t_0 + u_1_minus.*t_1 + u_2_minus.*t_2;
u_0_plus = s.*( u_0_minus + umt.*t_0 + (u_1_minus.*t_2 - u_2_minus.*t_1));
u_1_plus = s.*( u_1_minus + umt.*t_1 + (u_2_minus.*t_0 - u_0_minus.*t_2));
u_2_plus = s.*( u_2_minus + umt.*t_2 + (u_0_minus.*t_1 - u_1_minus.*t_0));
Ux = u_0_plus + E_0 + (u_1_plus.*t_2 - u_2_plus.*t_1);
Uy = u_1_plus + E_1 + (u_2_plus.*t_0 - u_0_plus.*t_2);
Uz = u_2_plus + E_2 + (u_0_plus.*t_1 - u_1_plus.*t_0);

% EQUIV:
% u_0_prime = u_0_minus + (u_1_minus.*t_2 - u_2_minus.*t_1);
% u_1_prime = u_1_minus + (u_2_minus.*t_0 - u_0_minus.*t_2);
% u_2_prime = u_2_minus + (u_0_minus.*t_1 - u_1_minus.*t_0);
% dot_u_prime_t = u_0_prime.*t_0 + u_1_prime.*t_1 + u_2_prime.*t_2;
% u_0_plus = s.*( u_0_prime + (dot_u_prime_t).*t_0 + (u_1_prime.*t_2 - u_2_prime.*t_1));
% u_1_plus = s.*( u_1_prime + (dot_u_prime_t).*t_1 + (u_2_prime.*t_0 - u_0_prime.*t_2));
% u_2_plus = s.*( u_2_prime + (dot_u_prime_t).*t_2 + (u_0_prime.*t_1 - u_1_prime.*t_0));
%
% %Now for the updated quantities
% Ux = u_0_plus + E_0;
% Uy = u_1_plus + E_1;
% Uz = u_2_plus + E_2;

%Save the output
gamma_n_plus_1 = sqrt(1+(Ux.*Ux + Uy.*Uy + Uz.*Uz)/(c*c));
Vx = Ux./gamma_n_plus_1;
Vy = Uy./gamma_n_plus_1;
Vz = Uz./gamma_n_plus_1;

%Interpolate Back
% interp_center_to_edge() Procedure
%[  A   B   C   D  ] - Centered Values
%    \ / \ / \ /
%[P   X   Y   Z   P] - Edge Values
% P = 0 -> fix then calculated in BC_J
%Uy = interp_center_to_edge(Uy,grid);
%Uz = interp_center_to_edge(Uz,grid);
%Vy = interp_center_to_edge(Vy,grid);
%Vz = interp_center_to_edge(Vz,grid);



end