%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Name: Grant Johnson
%Date: 5/17/2023
%Fluid algroithm, V (divergence U)

%Notes:
%-1D
% https://gkeyll.readthedocs.io/en/latest/dev/ssp-rk.html
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Update the quanitites Ux, Uy, Uz (t -> t + dt)
function [Ux,Uy,Uz,Vx,Vy,Vz,N,grid] = fluid_grad_U_FD(Ux,Uy,Uz,N,grid)

% Right and left 
Nx = grid.Nx;
grid.R = mod( linspace(1,Nx,Nx), Nx-1) + 1; %mod( linspace(1,Nx,Nx), Nx) + 1; %Good
grid.L = mod( linspace(-1,Nx-2,Nx), Nx-1) + 1; %mod( linspace(-1,Nx-2,Nx), Nx) + 1; %Good

%Update N & U with SSP, RK3:
N0 = N;
Ux0 = Ux;
Uy0 = Uy;
Uz0 = Uz;

%Stage 1:
[N, Ux, Uy, Uz] = stage1(N, Ux, Uy, Uz, grid);
%Stage 2:
[N, Ux, Uy, Uz] = stage2(N, Ux, Uy, Uz, N0, Ux0, Uy0, Uz0, grid);
%Stage 3:
[N, Ux, Uy, Uz] = stage3(N, Ux, Uy, Uz, N0, Ux0, Uy0, Uz0, grid);

%Save the output
gamma = sqrt(1+(Ux.*Ux + Uy.*Uy + Uz.*Uz)/(grid.c*grid.c));
Vx = Ux./gamma;
Vy = Uy./gamma;
Vz = Uz./gamma;

end


%SSP-RK3 (stage 3)
function [N, Ux, Uy, Uz] = stage3(N, Ux, Uy, Uz, N0, Ux0, Uy0, Uz0, grid)

%First stage calculation
[N_star, Ux_star, Uy_star, Uz_star] = euler(N, Ux, Uy, Uz, grid);
N = (1/3)*N0 + (2/3)*N_star;
Ux = (1/3)*Ux0 + (2/3)*Ux_star;
Uy = (1/3)*Uy0 + (2/3)*Uy_star;
Uz = (1/3)*Uz0 + (2/3)*Uz_star;

end

%SSP-RK3 (stage 2)
function [N, Ux, Uy, Uz] = stage2(N, Ux, Uy, Uz, N0, Ux0, Uy0, Uz0, grid)

%First stage calculation
[N_star, Ux_star, Uy_star, Uz_star] = euler(N, Ux, Uy, Uz, grid);
N = (3/4)*N0 + (1/4)*N_star;
Ux = (3/4)*Ux0 + (1/4)*Ux_star;
Uy = (3/4)*Uy0 + (1/4)*Uy_star;
Uz = (3/4)*Uz0 + (1/4)*Uz_star;

end

%SSP-RK3 (stage 1)
function [N, Ux, Uy, Uz] = stage1(N, Ux, Uy, Uz, grid)

%First stage calculation
[N_star, Ux_star, Uy_star, Uz_star] = euler(N, Ux, Uy, Uz, grid);
N = N_star;
Ux = Ux_star;
Uy = Uy_star;
Uz = Uz_star;

end

function [N_star, Ux_star, Uy_star, Uz_star] = euler(N, Ux, Uy, Uz, grid)

%Simplify the equations;
c = grid.dt/grid.dx;
R = grid.R;
L = grid.L;
RR = grid.R(grid.R);
LL = grid.L(grid.L);

% Compute the velocity
Vx = Ux./sqrt(1.0 + (Ux.*Ux + Uy.*Uy + Uz.*Uz)/(grid.c*grid.c) );

% Compute a
Zero = zeros(size(Vx));
a_p = max(Vx,Zero);
a_m = min(Vx,Zero);

% Compute u_+ and u_- (first order)
f_m1 = N - N(L);
f_p1 = N(R) - N;
f_m2 = Ux - Ux(L);
f_p2 = Ux(R) - Ux;
f_m3 = Uy - Uy(L);
f_p3 = Uy(R) - Uy;
f_m4 = Uz - Uz(L);
f_p4 = Uz(R) - Uz;

% Compute u_+ and u_- (2nd order)
% f_m1 = 0.5*(3.0*N - 4.0*N(L) + N(LL));
% f_p1 = 0.5*(-N(RR) + 4.0*N(R) - 3.0*N);
% f_m2 = 0.5*(3.0*Ux - 4.0*Ux(L) + Ux(LL));
% f_p2 = 0.5*(-Ux(RR) + 4.0*Ux(R) - 3.0*Ux);
% f_m3 = 0.5*(3.0*Uy - 4.0*Uy(L) + Uy(LL));
% f_p3 = 0.5*(-Uy(RR) + 4.0*Uy(R) - 3.0*Uy);
% f_m4 = 0.5*(3.0*Uz - 4.0*Uz(L) + Uz(LL));
% f_p4 = 0.5*(-Uz(RR) + 4.0*Uz(R) - 3.0*Uz);

% Update to get the star values
N_star  = N  - c*( a_p.*f_m1 + f_p1.*a_m );
Ux_star = Ux - c*( a_p.*f_m2 + f_p2.*a_m );
Uy_star = Uy - c*( a_p.*f_m3 + f_p3.*a_m );
Uz_star = Uz - c*( a_p.*f_m4 + f_p4.*a_m );

end

