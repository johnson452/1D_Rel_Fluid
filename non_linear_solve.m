% Nonlinear solve the momentum conservation equations
function [u_hat_L_vector] = non_linear_solve(N,Ux,Uy,Uz,L,Nx,c_light)

% Nonlinear solve:
%Constants
c = c_light;
rho_max = max(N(L),N);
rhoL = N(L)./rho_max;
rhoR = N./rho_max;

% SCALE
uLx = Ux(L)/c;
uRx = Ux/c;
uLy = Uy(L)/c;
uRy = Uy/c;
uLz = Uz(L)/c;
uRz = Uz/c;
c = 1.0;
gammaL = sqrt(1 + (uLx.*uLx + uLy.*uLy + uLz.*uLz )/(c*c));
gammaR = sqrt(1 + (uRx.*uRx + uRy.*uRy + uRz.*uRz )/(c*c));
vLx = uLx./gammaL;
vRx = uRx./gammaR;


%%%% (RELATIVSITIC 3V) %%%%
Ux_hat = zeros(1,Nx);
Uy_hat = zeros(1,Nx);
Uz_hat = zeros(1,Nx);
for i = 1:Nx
    if (uLx(i) < 0 && 0 < uRx(i)) || (uLx(i) == 0 && 0 == uRx(i))
        Ux_hat(i) = 0;
        Uy_hat(i) = 0;
        Uz_hat(i) = 0;
    else
        [Ux_hat(i), Uy_hat(i), Uz_hat(i) ] = NL_solve(c,rhoR(i),rhoL(i),uLx(i),uRx(i),uLy(i),uRy(i),uLz(i),uRz(i),vLx(i),vRx(i),gammaL(i),gammaR(i));
    end
end

%Return the vector of momentums
u_hat_L_vector = [Ux_hat; Uy_hat; Uz_hat ];

% Rescale
u_hat_L_vector = u_hat_L_vector*c_light;


%Check u0
if check(u_hat_L_vector) == 0
    fprintf("IC Failed!\n");
end

end


% Conservation of momentum, returns U
function F = momentum_cons_eqs(x,c,rhoR,rhoL,uLx,uRx,uLy,uRy,uLz,uRz,vLx,vRx)

% Left and right masses accreted
% mL = abs(rhoL*(vLx - x(1)));
% mR = abs(rhoR*(x(1) - vRx));
% F(1) = (mL + mR)*x(1)/sqrt(1.0 - (x(1)*x(1) + x(2)*x(2) + x(3)*x(3))/(c*c)) - mL*uLx - mR*uRx;
% F(2) = (mL + mR)*x(2)/sqrt(1.0 - (x(1)*x(1) + x(2)*x(2) + x(3)*x(3))/(c*c)) - mL*uLy - mR*uRy;
% F(3) = (mL + mR)*x(3)/sqrt(1.0 - (x(1)*x(1) + x(2)*x(2) + x(3)*x(3))/(c*c)) - mL*uLz - mR*uRz;

Usx = x(1);
Usy = x(2);
Usz = x(3);
Vsx = Usx/sqrt(1.0 + (Usx*Usx + Usy*Usy + Usz*Usz)/(c*c));
f1 = Usx*( (rhoR-rhoL)*Vsx + rhoL*vLx - rhoR*vRx ) + rhoL*uLx*(Vsx - vLx) + rhoR*uRx*(vRx - Vsx);
f2 = Usy*( (rhoR-rhoL)*Vsx + rhoL*vLx - rhoR*vRx ) + rhoL*uLy*(Vsx - vLx) + rhoR*uRy*(vRx - Vsx);
f3 = Usz*( (rhoR-rhoL)*Vsx + rhoL*vLx - rhoR*vRx ) + rhoL*uLz*(Vsx - vLx) + rhoR*uRz*(vRx - Vsx);
F = [f1,f2,f3];
end

%Condition test
function [cond] = conditions(v,~,~,u,uL,uR)

%Check that the velocity is logical
cond = 0;
sz = max(size(v));
for i = 1:sz
    %if min(vL(i),vR(i)) <= v(i) && v(i) <= max(vL(i),vR(i))
    %else
    %fprintf("Spooky Velocity (Not nessesarily a true cond.)\n")
    %end

    if min(uL(i),uR(i)) <= u(i) && u(i) <= max(uL(i),uR(i))
    else
        %fprintf("Violation of the momentum\n")
        cond = 1; % Failed
    end
end
end


function [ux,uy,uz] = NL_solve(c,rhoR,rhoL,uLx,uRx,uLy,uRy,uLz,uRz,vLx,vRx,gammaL,gammaR,options)

%Guesses
vLy = uLy./gammaL;
vRy = uRy./gammaR;
vLz = uLz./gammaL;
vRz = uRz./gammaR;
vR = [vRx,vRy,vRz];
vL = [vLx,vLy,vLz];
uR = [uRx,uRy,uRz];
uL = [uLx,uLy,uLz];
u0 = (sqrt(rhoL)*uL + sqrt(rhoR)*uR )/(sqrt(rhoL) + sqrt(rhoR));
%u0 = 0.5*[uLx+uRx,uLy+uRy,uLz+uRz]; % Might need a u0 guess first


%Check u0
if check(u0) == 0
    fprintf("IC Failed!\n");
end

global_tol = max(rhoL,rhoR)*max(abs(u0));
Final_tols = [];
newton_tol = global_tol*1e-20;

[v,cond4] = direct_ridders(u0, c,rhoR,rhoL,uLx,uRx,uLy,uRy,uLz,uRz,vLx,vRx,newton_tol);
u = v./sqrt(1-(v(1)^2 + v(2)^2 + v(3)^2)/c^2);
F_x = momentum_cons_eqs(v,c,rhoR,rhoL,uLx,uRx,uLy,uRy,uLz,uRz,vLx,vRx);
%fprintf("Ridders F(x) = (0) is actually : %e ",F_x)
%fprintf("\n")
cond4 = conditions(v,vL,vR,u,uL,uR);

%Stop if cond 4 breaks
if cond4 == 1 && max(abs(F_x)) > 1e-6
    fprintf("Ridders' fails\n")
elseif isnan(u)
    u = [0,0,0];
else
    cond4 = 0; 
end

cond = 0;
cond3 = 0;
if cond4 == 1
    [u,cond3] = direct_iterations(u0,c,rhoR,rhoL,uLx,uRx,uLy,uRy,uLz,uRz,vLx,vRx,newton_tol);
    v = u./sqrt(1+(u(1)^2 + u(2)^2 + u(3)^2)/c^2);
    cond = conditions(v,vL,vR,u,uL,uR);
    F_x = momentum_cons_eqs(u,c,rhoR,rhoL,uLx,uRx,uLy,uRy,uLz,uRz,vLx,vRx);
    %fprintf("Direct max(abs(F_x)): %e of tol: %e, ratio: %e\n",max(abs(F_x)),global_tol,max(abs(F_x))/global_tol);
    Final_tols = [Final_tols,max(abs(F_x))/global_tol];
    u_vec = u;
end


% Attemp newton solve
cond2 = 0;
if cond == 1  || cond3 == 1
    [u,cond2] = newton_iterations(u0, c,rhoR,rhoL,uLx,uRx,uLy,uRy,uLz,uRz,vLx,vRx,newton_tol);
    v = u./sqrt(1+(u(1)^2 + u(2)^2 + u(3)^2)/c^2);
    cond = conditions(v,vL,vR,u,uL,uR);
    F_x = momentum_cons_eqs(u,c,rhoR,rhoL,uLx,uRx,uLy,uRy,uLz,uRz,vLx,vRx);
    %fprintf("Newton max(abs(F_x)): %e of tol: %e, ratio: %e\n",max(abs(F_x)),global_tol,max(abs(F_x))/global_tol);
    Final_tols = [Final_tols,max(abs(F_x))/global_tol];
    %u_vec = u;
    u_vec = [u_vec;u];
end
if cond == 1  || cond2 == 1
    u0 = [uLx,uLy,uLz];
    [u,cond2] = newton_iterations(u0, c,rhoR,rhoL,uLx,uRx,uLy,uRy,uLz,uRz,vLx,vRx,newton_tol);
    v = u./sqrt(1+(u(1)^2 + u(2)^2 + u(3)^2)/c^2);
    cond = conditions(v,vL,vR,u,uL,uR);
    F_x = momentum_cons_eqs(u,c,rhoR,rhoL,uLx,uRx,uLy,uRy,uLz,uRz,vLx,vRx);
    %fprintf("Newton max(abs(F_x)): %e of tol: %e, ratio: %e\n",max(abs(F_x)),global_tol,max(abs(F_x))/global_tol);
    Final_tols = [Final_tols,max(abs(F_x))/global_tol];
    u_vec = [u_vec;u];
end

if cond == 1  || cond2 == 1
    u0 = [uRx,uRy,uRz];
    [u,cond2] = newton_iterations(u0, c,rhoR,rhoL,uLx,uRx,uLy,uRy,uLz,uRz,vLx,vRx,newton_tol);
    v = u./sqrt(1+(u(1)^2 + u(2)^2 + u(3)^2)/c^2);
    cond = conditions(v,vL,vR,u,uL,uR);
    F_x = momentum_cons_eqs(u,c,rhoR,rhoL,uLx,uRx,uLy,uRy,uLz,uRz,vLx,vRx);
    %fprintf("Newton max(abs(F_x)): %e of tol: %e, ratio: %e\n",max(abs(F_x)),global_tol,max(abs(F_x))/global_tol);
    Final_tols = [Final_tols,max(abs(F_x))/global_tol];
    u_vec = [u_vec;u];
end


if cond == 1  || cond2 == 1 % || max(abs(F_x))/global_tol > 1e-15
    %[S,parameters,conditions]

    %Options for fsolve
    %options = optimset('Display','off');
    options = optimoptions('fsolve', 'TolX', 1e-20); % Set a smaller tolerance (e.g., 1e-6)
    options.OptimalityTolerance = 1.000000e-20;
    options.FunctionTolerance = 1.000000e-20;
    options.Display = "off";

    u0 = (sqrt(rhoL)*uL + sqrt(rhoR)*uR )/(sqrt(rhoL) + sqrt(rhoR));
    F0_x = momentum_cons_eqs(u0,c,rhoR,rhoL,uLx,uRx,uLy,uRy,uLz,uRz,vLx,vRx);
    if check(F0_x)
        u = fsolve(@(x) momentum_cons_eqs(x,c,rhoR,rhoL,uLx,uRx,uLy,uRy,uLz,uRz,vLx,vRx), u0, options);
    end
    F_x = momentum_cons_eqs(u,c,rhoR,rhoL,uLx,uRx,uLy,uRy,uLz,uRz,vLx,vRx);
    %fprintf("Fsolve max(abs(F_x)): %e of tol: %e, ratio: %e\n",max(abs(F_x)),global_tol,max(abs(F_x))/global_tol);
    v = u./sqrt(1+(u(1)^2 + u(2)^2 + u(3)^2)/c^2);
    cond = conditions(v,vL,vR,u,uL,uR);
    Final_tols = [Final_tols,max(abs(F_x))/global_tol];
    u_vec = [u_vec;u];
end
if cond == 1
    % New inital guess
    u0 = [uLx,uLy,uLz];
    %v0 = u0./sqrt(1+(u0(1)^2 + u0(2)^2 + u0(3)^2)/c^2);
    F0_x = momentum_cons_eqs(u0,c,rhoR,rhoL,uLx,uRx,uLy,uRy,uLz,uRz,vLx,vRx);
    if check(F0_x)
        u = fsolve(@(x) momentum_cons_eqs(x,c,rhoR,rhoL,uLx,uRx,uLy,uRy,uLz,uRz,vLx,vRx), u0, options);
    end
    F_x = momentum_cons_eqs(u,c,rhoR,rhoL,uLx,uRx,uLy,uRy,uLz,uRz,vLx,vRx);
    %fprintf("Fsolve max(abs(F_x)): %e of tol: %e, ratio: %e\n",max(abs(F_x)),global_tol,max(abs(F_x))/global_tol);
    v = u./sqrt(1+(u(1)^2 + u(2)^2 + u(3)^2)/c^2);
    cond = conditions(v,vL,vR,u,uL,uR);
    Final_tols = [Final_tols,max(abs(F_x))/global_tol];
    u_vec = [u_vec;u];
end
if cond == 1
    % New inital guess
    u0 = [uRx,uRy,uRz];
    %v0 = u0./sqrt(1+(u0(1)^2 + u0(2)^2 + u0(3)^2)/c^2);
    F0_x = momentum_cons_eqs(u0,c,rhoR,rhoL,uLx,uRx,uLy,uRy,uLz,uRz,vLx,vRx);
    if check(F0_x)
        u = fsolve(@(x) momentum_cons_eqs(x,c,rhoR,rhoL,uLx,uRx,uLy,uRy,uLz,uRz,vLx,vRx), u0, options);
    end
    F_x = momentum_cons_eqs(u,c,rhoR,rhoL,uLx,uRx,uLy,uRy,uLz,uRz,vLx,vRx);
    %fprintf("Fsolve max(abs(F_x)): %e of tol: %e, ratio: %e\n",max(abs(F_x)),global_tol,max(abs(F_x))/global_tol);
    v = u./sqrt(1+(u(1)^2 + u(2)^2 + u(3)^2)/c^2);
    cond = conditions(v,vL,vR,u,uL,uR);
    Final_tols = [Final_tols,max(abs(F_x))/global_tol];
    u_vec = [u_vec;u];
end
if cond == 1
    %fprintf("Issue on NL Solve: Using Best\n");
    indx = find_indx(Final_tols);
    u = u_vec(indx,:);
    for it_tols = 1:max(size(Final_tols))
        %    fprintf("Tol %d: %e ",it_tols,Final_tols(it_tols))
    end
    %fprintf("\n")
    %fprintf("Selected index %d: %e\n",indx,Final_tols(indx))

    if Final_tols(indx) > 1e-5
        fprintf("True Failure of the root solving algo.")
        %u0 = 0.5*[uLx+uRx,uLy+uRy,uLz+uRz];
    end
    %diags(v,u,rhoR,rhoL,uLx,uRx,uLy,uRy,uLz,uRz,vLx,vRx,vLy,vRy,vLz,vRz,gammaL,gammaR);

    %Total failure diagnostic:
    %total_failure(v,u,rhoR,rhoL,uLx,uRx,uLy,uRy,uLz,uRz,vLx,vRx,vLy,vRy,vLz,vRz,gammaL,gammaR)

    % Assume defaults:
    %u = 0.5*[uLx+uRx,uLy+uRy,uLz+uRz];
end

%Split the output
ux = u(1);
uy = u(2);
uz = u(3);

end


function diags(v,u,rhoR,rhoL,uLx,uRx,uLy,uRy,uLz,uRz,vLx,vRx,vLy,vRy,vLz,vRz,gammaL,gammaR)

clf()

% Left and right
N = 100;
L = zeros(1,N);
R = zeros(1,N);
for i = 1:N
    if i <= N/2
        L(i) = 1;
    elseif i > N/2
        R(i) = 1;
    end
end
x = linspace(-1,1,N);

%Build the vectors
rho = L.*rhoL + R.*rhoR;
Ux = L.*uLx + R.*uRx;
Uy = L.*uLy + R.*uRy;
Uz = L.*uLz + R.*uRz;
Vx = L.*vLx + R.*vRx;
Vy = L.*vLy + R.*vRy;
Vz = L.*vLz + R.*vRz;
gamma = L.*gammaL + R.*gammaR;

%Plot the results
subplot(2,2,1)
plot(x,rho,"red")
title("rho")

subplot(2,2,2)
plot(x,Ux,"red")
hold on
plot(x,Uy,"blue")
hold on
plot(x,Uz,"black")
hold on
plot(0.0,u(1),"*red")
hold on
plot(0.0,u(2),"*blue")
hold on
plot(0.0,u(3),"*black")
title("U")
legend("Ux","Uy","Uz")


subplot(2,2,3)
plot(x,Vx,"red")
hold on
plot(x,Vy,"blue")
hold on
plot(x,Vz,"black")
hold on
plot(0.0,v(1),"*red")
hold on
plot(0.0,v(2),"*blue")
hold on
plot(0.0,v(3),"*black")
title("V")
legend("Vx","Vy","Vz")

subplot(2,2,4)
plot(x,gamma,"red")
title("Gamma")

%Check on conditions we expect
% vR = [vRx,vRy,vRz];
% vL = [vLx,vLy,vLz];
% uR = [uRx,uRy,uRz];
% uL = [uLx,uLy,uLz];
% conditions(v,vL,vR,u,uL,uR)

pause(0.0)

end


% Newton functions:
function [u0,cond] = newton_iterations(u0,c,rhoR,rhoL,uLx,uRx,uLy,uRy,uLz,uRz,vLx,vRx,newton_tol)

%Call momentum conservation equation with v0
u0 = u0';
tol = newton_tol*2;
iter = 0;
cond = 0;
while (tol > newton_tol)
    J_inverse = J_inv(u0,c,rhoR,rhoL,uLx,uRx,uLy,uRy,uLz,uRz,vLx,vRx);
    F_x = momentum_cons_eqs(u0,c,rhoR,rhoL,uLx,uRx,uLy,uRy,uLz,uRz,vLx,vRx);
    u0_new = u0 - J_inverse*F_x';
    %tol = max(abs(u0_new-u0));
    tol = max(abs(F_x));
    u0 = u0_new;
    iter = iter + 1;
    v0 = u0./sqrt(1 + (u0(1)^2 + u0(2)^2 + u0(3)^2)/(c^2));
    %fprintf("Error in iter %d: %e (max F(x))\n",iter,tol)
    if max(v0) > c || iter > 10
        %fprintf("Sol. diverged! (max v0/c: %f) \n",max(v0)/c)
        tol = 0.0;
        cond = 1;
    end
end

u0 = u0';
%fprintf("Converged in %d iterations\n",iter)

end


% Conservation of momentum,
function [J_in] = J_inv(x,c,rhoR,rhoL,uLx,uRx,uLy,uRy,uLz,uRz,vLx,vRx)

% Left and right masses accreted
Usx = x(1);
Usy = x(2);
Usz = x(3);
gamma = sqrt(1 + (Usx^2 + Usy^2 + Usz^2)/c^2);
J = [[-(Usx^3*rhoL - Usx^3*rhoR - c^2*rhoL*uLx + c^2*rhoR*uRx + 2*Usx*Usy^2*rhoL - 2*Usx*Usy^2*rhoR + 2*Usx*Usz^2*rhoL - 2*Usx*Usz^2*rhoR + 2*Usx*c^2*rhoL - 2*Usx*c^2*rhoR - Usy^2*rhoL*uLx - Usz^2*rhoL*uLx + Usy^2*rhoR*uRx + Usz^2*rhoR*uRx - c^2*gamma^3*rhoL*vLx + c^2*gamma^3*rhoR*vRx)/(c^2*gamma^3),                                                                                                                          (Usx*Usy*(Usx*rhoL - Usx*rhoR - rhoL*uLx + rhoR*uRx))/(c^2*gamma^3),                                                                                                                          (Usx*Usz*(Usx*rhoL - Usx*rhoR - rhoL*uLx + rhoR*uRx))/(c^2*gamma^3)];...
    [                                                                                                                                                                                                                    -((Usy^2 + Usz^2 + c^2)*(Usy*rhoL - Usy*rhoR - rhoL*uLy + rhoR*uRy))/(c^2*gamma^3), -(Usx^3*rhoL - Usx^3*rhoR + Usx*Usz^2*rhoL - Usx*Usz^2*rhoR + Usx*c^2*rhoL - Usx*c^2*rhoR - c^2*gamma^3*rhoL*vLx + c^2*gamma^3*rhoR*vRx + Usx*Usy*rhoL*uLy - Usx*Usy*rhoR*uRy)/(c^2*gamma^3),                                                                                                                          (Usx*Usz*(Usy*rhoL - Usy*rhoR - rhoL*uLy + rhoR*uRy))/(c^2*gamma^3)];...
    [                                                                                                                                                                                                                    -((Usy^2 + Usz^2 + c^2)*(Usz*rhoL - Usz*rhoR - rhoL*uLz + rhoR*uRz))/(c^2*gamma^3),                                                                                                                          (Usx*Usy*(Usz*rhoL - Usz*rhoR - rhoL*uLz + rhoR*uRz))/(c^2*gamma^3), -(Usx^3*rhoL - Usx^3*rhoR + Usx*Usy^2*rhoL - Usx*Usy^2*rhoR + Usx*c^2*rhoL - Usx*c^2*rhoR - c^2*gamma^3*rhoL*vLx + c^2*gamma^3*rhoR*vRx + Usx*Usz*rhoL*uLz - Usx*Usz*rhoR*uRz)/(c^2*gamma^3)]];
J_in = inv(J);
end


function [minIndex] = find_indx(myArray)

% Initialize variables to store the minimum value and its index
minValue = myArray(1);  % Assume the first element is the minimum
minIndex = 1;          % Index of the first element

% Iterate through the array to find the minimum value and its index
for i = 2:length(myArray)
    if myArray(i) < minValue
        minValue = myArray(i);
        minIndex = i;
    end
end
end


% Newton functions:
function [u0,cond] = direct_iterations(u0,c,rhoR,rhoL,uLx,uRx,uLy,uRy,uLz,uRz,vLx,vRx,newton_tol)

%Check that we don't divide by zero: (Vx-hat is defined)

%Call momentum conservation equation with v0
Ux0 = u0(1);
Uy0 = u0(2);
Uz0 = u0(3);
tol = newton_tol*2;
iter = 0;
cond = 0;

%Compute the constants:
dp = rhoR - rhoL;
a = rhoL*vLx - rhoR*vRx;
Ax = - rhoL*vLx*uLx + rhoR*vRx*uRx;
Ay = - rhoL*vLx*uLy + rhoR*vRx*uRy;
Az = - rhoL*vLx*uLz + rhoR*vRx*uRz;
Bx = rhoL*uLx - rhoR*uRx;
By = rhoL*uLy - rhoR*uRy;
Bz = rhoL*uLz - rhoR*uRz;

% Checks to make sure we don't div by zero:
if uLz == 0.0 && uRz == 0.0 && uLy == 0.0 && uRy == 0.0 && a ~= 0.0

    % Print which case:
    %fprintf("Case 1: (Uy == 0) and (Uz == 0)\n")

    %1D case:
    while (tol > newton_tol) || iter < 5

        % Compute values with old u0
        Vx = (-Ux0*a - Ax)/(Ux0*dp + Bx);

        % Newton Iterate: (From CAS/F_x_prime_Uy_for_3D_rel_solve.m)
        F_Ux_prime = F_Ux_non_zero_prime(a, c, dp, Ax, Bx, Ux0);
        F_Ux = Vx -  Ux0 / sqrt(1 + (Ux0^2)/c^2);
        Ux0 = Ux0 - F_Ux/F_Ux_prime;

        % Check the tolerance
        tol = max(abs(F_Ux));
        iter = iter + 1;
        if iter > 10 || isnan(Ux0)
            %fprintf("Sol. diverged!\n")
            tol = 0.0;
            cond = 1;
        end
    end

    %Convert Back to Ux Uy Uz:
    Ux = Ux0;
    Uy = 0;
    Uz = 0;

elseif uLy == 0.0 && uRy == 0.0  && a ~= 0.0

    % Print which case:
    %fprintf("Case 2: (Uy == 0)\n")

    %2D case (uy == 0):
    while (tol > newton_tol) || iter < 5

        % Compute values with old u0
        Vx = (-Uz0*a - Ay)/(Uz0*dp + Bz);
        Ux = (Ax + Bx*Vx)/(-dp*Vx - a);


        % Newton Iterate: (From CAS/F_x_prime_Uy_for_3D_rel_solve.m)
        F_Ux_prime =  F_Ux_Uz_non_zero_prime(a, c, dp, Ax, Ay, Bx, Bz, Uz0);
        F_Ux = Vx -  Ux / sqrt(1 + (Ux^2 + Uz0^2)/c^2);
        Uz0 = Uz0 - F_Ux/F_Ux_prime;

        % Check the tolerance
        tol = max(abs(F_Ux));
        iter = iter + 1;
        if iter > 10 || isnan(Uz0)
            %fprintf("Sol. diverged!\n")
            tol = 0.0;
            cond = 1;
        end
    end

    %Convert Back to Ux Uy Uz:
    Vx = (-Uz0*a - Ay)/(Uz0*dp + Bz);
    Ux = (Ax + Bx*Vx)/(-dp*Vx - a);
    Uy = 0;
    Uz = Uz0;

elseif uLz == 0.0 && uRz == 0.0  && a ~= 0.0

    % Print which case:
    %fprintf("Case 3: (Uz == 0)\n")

    %2D case (uz == 0):
    while (tol > newton_tol) || iter < 5

        % Compute values with old u0
        Vx = (-Uy0*a - Ay)/(Uy0*dp + By);
        Ux = (Ax + Bx*Vx)/(-dp*Vx - a);


        % Newton Iterate: (From CAS/F_x_prime_Uy_for_3D_rel_solve.m)
        F_Ux_prime =  F_Ux_Uy_non_zero_prime(a, c, dp, Ax, Ay, Bx, By, Uy0);
        F_Ux = Vx -  Ux / sqrt(1 + (Ux^2 + Uy0^2)/c^2);
        Uy0 = Uy0 - F_Ux/F_Ux_prime;

        % Check the tolerance
        tol = max(abs(F_Ux));
        iter = iter + 1;
        if iter > 10 || isnan(Uy0)
            %fprintf("Sol. diverged!\n")
            tol = 0.0;
            cond = 1;
        end
    end

    %Convert Back to Ux Uy Uz:
    Vx = (-Uy0*a - Ay)/(Uy0*dp + By);
    Ux = (Ax + Bx*Vx)/(-dp*Vx - a);
    Uy = Uy0;
    Uz = 0;

elseif a ~= 0.0 && By ~= 0.0

    % Print which case:
    %fprintf("Case 4: No Zeros \n")

    %Full 3D case
    while (tol > newton_tol) || iter < 5

        % Compute values with old u0
        Vx = (-Uy0*a - Ay)/(Uy0*dp + By);
        Ux = (Ax + Bx*Vx)/(-dp*Vx - a);
        Uz = (Az + Bz*Vx)/(-dp*Vx - a);

        % Newton Iterate: (From CAS/F_x_prime_Uy_for_3D_rel_solve.m)
        F_Uy_prime = F_Ux_Uy_Uz_non_zero_prime(a, c, dp, Ax, Ay, Az, Bx, By, Bz, Uy0);
        F_Uy = Vx -  Ux / sqrt(1 + (Ux^2 + Uy0^2 + Uz^2)/c^2);
        Uy0 = Uy0 - F_Uy/F_Uy_prime;

        % Check the tolerance
        tol = max(abs(F_Uy));
        iter = iter + 1;
        if iter > 10 || isnan(Uy0)
            %fprintf("Sol. diverged!\n")
            tol = 0.0;
            cond = 1;
        end
    end

    %Convert Back to Ux Uy Uz:
    Vx = (-Uy0*a - Ay)/(Uy0*dp + By);
    Ux = (Ax + Bx*Vx)/(-dp*Vx - a);
    Uy = Uy0;
    Uz = (Az + Bz*Vx)/(-dp*Vx - a);

else
    % If not the right criterion, use simpler solves:
    %fprintf("Undefined Cases!\n")
    Ux = 0;
    Uy = 0;
    Uz = 0;
end


%Compute the exact values to hand back
gamma = sqrt(1 + (Ux^2 + Uy^2 + Uz^2)/c^2);
Vx = Ux/gamma;
Vy = Uy/gamma;
Vz = Uz/gamma;
v0 = [Vx,Vy,Vz];
u0 = [Ux,Uy,Uz];

%If nan's etc then set to zero
if check(u0) == 0
    u0 = [0,0,0];
end

%fprintf("(Direct) Converged in %d iterations\n",iter)

end

% Compute the derivative of F(Uy)
function [val] = F_Ux_Uy_Uz_non_zero_prime(a, c, dp, Ax, Ay, Az, Bx, By, Bz, Uy)
val = (dp*(Ay + Uy*a))/(By + Uy*dp)^2 - Bx/((((Ay*Bx - Ax*By + Bx*Uy*a - Ax*Uy*dp)^2/(By*a - Ay*dp)^2 + (Ay*Bz - Az*By + Bz*Uy*a - Az*Uy*dp)^2/(By*a - Ay*dp)^2 + Uy^2)/c^2 + 1)^(1/2)*(By + Uy*dp)) - a/(By + Uy*dp) + ((Ay*Bx - Ax*By + Bx*Uy*a - Ax*Uy*dp)*(Uy*Ay^2*dp^2 - Ay*Az*Bz*dp + Ay*Bx^2*a - Ax*Ay*Bx*dp - 2*Uy*Ay*By*a*dp + Ay*Bz^2*a + Az^2*By*dp + Uy*Az^2*dp^2 - Az*By*Bz*a - 2*Uy*Az*Bz*a*dp + Uy*Bx^2*a^2 - Ax*Bx*By*a - 2*Uy*Ax*Bx*a*dp + Uy*By^2*a^2 + Ax^2*By*dp + Uy*Bz^2*a^2 + Uy*Ax^2*dp^2))/(c^2*(((Ay*Bx - Ax*By + Bx*Uy*a - Ax*Uy*dp)^2/(By*a - Ay*dp)^2 + (Ay*Bz - Az*By + Bz*Uy*a - Az*Uy*dp)^2/(By*a - Ay*dp)^2 + Uy^2)/c^2 + 1)^(3/2)*(By*a - Ay*dp)^3) - (dp*(Ay*Bx - Ax*By + Bx*Uy*a - Ax*Uy*dp))/((((Ay*Bx - Ax*By + Bx*Uy*a - Ax*Uy*dp)^2/(By*a - Ay*dp)^2 + (Ay*Bz - Az*By + Bz*Uy*a - Az*Uy*dp)^2/(By*a - Ay*dp)^2 + Uy^2)/c^2 + 1)^(1/2)*(By*a - Ay*dp)*(By + Uy*dp));
end

function [val] = F_Ux_non_zero_prime(a, c, dp, Ax, Bx, Ux)
val = Ux^2/(c^2*(Ux^2/c^2 + 1)^(3/2)) - a/(Bx + Ux*dp) - 1/(Ux^2/c^2 + 1)^(1/2) + (dp*(Ax + Ux*a))/(Bx + Ux*dp)^2;
end

function [val] = F_Ux_Uy_non_zero_prime(a, c, dp, Ax, Ay, Bx, By, Uy)
val = (dp*(Ay + Uy*a))/(By + Uy*dp)^2 - Bx/((By + Uy*dp)*(((Ay*Bx - Ax*By + Bx*Uy*a - Ax*Uy*dp)^2/(By*a - Ay*dp)^2 + Uy^2)/c^2 + 1)^(1/2)) - a/(By + Uy*dp) + ((Ay*Bx - Ax*By + Bx*Uy*a - Ax*Uy*dp)*(Uy*Ay^2*dp^2 + Ay*Bx^2*a - Ax*Ay*Bx*dp - 2*Uy*Ay*By*a*dp + Uy*Bx^2*a^2 - Ax*Bx*By*a - 2*Uy*Ax*Bx*a*dp + Uy*By^2*a^2 + Ax^2*By*dp + Uy*Ax^2*dp^2))/(c^2*(By*a - Ay*dp)^3*(((Ay*Bx - Ax*By + Bx*Uy*a - Ax*Uy*dp)^2/(By*a - Ay*dp)^2 + Uy^2)/c^2 + 1)^(3/2)) - (dp*(Ay*Bx - Ax*By + Bx*Uy*a - Ax*Uy*dp))/((By*a - Ay*dp)*(By + Uy*dp)*(((Ay*Bx - Ax*By + Bx*Uy*a - Ax*Uy*dp)^2/(By*a - Ay*dp)^2 + Uy^2)/c^2 + 1)^(1/2));
end

function [val] = F_Ux_Uz_non_zero_prime(a, c, dp, Ax, Ay, Bx, Bz, Uz)
val = (dp*(Ay + Uz*a))/(Bz + Uz*dp)^2 - Bx/((Bz + Uz*dp)*(((Ay*Bx - Ax*Bz + Bx*Uz*a - Ax*Uz*dp)^2/(Bz*a - Ay*dp)^2 + Uz^2)/c^2 + 1)^(1/2)) - a/(Bz + Uz*dp) + ((Ay*Bx - Ax*Bz + Bx*Uz*a - Ax*Uz*dp)*(Uz*Ay^2*dp^2 + Ay*Bx^2*a - Ax*Ay*Bx*dp - 2*Uz*Ay*Bz*a*dp + Uz*Bx^2*a^2 - Ax*Bx*Bz*a - 2*Uz*Ax*Bx*a*dp + Uz*Bz^2*a^2 + Ax^2*Bz*dp + Uz*Ax^2*dp^2))/(c^2*(Bz*a - Ay*dp)^3*(((Ay*Bx - Ax*Bz + Bx*Uz*a - Ax*Uz*dp)^2/(Bz*a - Ay*dp)^2 + Uz^2)/c^2 + 1)^(3/2)) - (dp*(Ay*Bx - Ax*Bz + Bx*Uz*a - Ax*Uz*dp))/((Bz*a - Ay*dp)*(Bz + Uz*dp)*(((Ay*Bx - Ax*Bz + Bx*Uz*a - Ax*Uz*dp)^2/(Bz*a - Ay*dp)^2 + Uz^2)/c^2 + 1)^(1/2));
end


function [bool_ret] = check(x)

% 1 is okay to continue
bool_ret = 1;

% If it fails conditions then return bool-ret = 0
if max(isnan(x),[],"all") == 1
    bool_ret = 0;
end
if min(isfinite(x),[],"all") == 0
    bool_ret = 0;
end
if ~isreal(x)
    bool_ret = 0;
end

if bool_ret == 0
    %fprintf("F_solve failed I.C.\n")
end

end







%%%% Ridder functions %%%%

% Newton functions:
function [v0,cond] = direct_ridders(~,c,rhoR,rhoL,uLx,uRx,uLy,uRy,uLz,uRz,vLx,vRx,newton_tol)

%Check that we don't divide by zero: (Vx-hat is defined)

%Call momentum conservation equation with v0
tol = newton_tol*2;

%Compute the constants:
dp = rhoR - rhoL;
a = rhoL*vLx - rhoR*vRx;
Ax = - rhoL*vLx*uLx + rhoR*vRx*uRx;
Ay = - rhoL*vLx*uLy + rhoR*vRx*uRy;
Az = - rhoL*vLx*uLz + rhoR*vRx*uRz;
Bx = rhoL*uLx - rhoR*uRx;
By = rhoL*uLy - rhoR*uRy;
Bz = rhoL*uLz - rhoR*uRz;

% Checks to make sure we don't div by zero:
if uLz == 0.0 && uRz == 0.0 && uLy == 0.0 && uRy == 0.0 && a ~= 0.0

    % Print which case:
    %fprintf("Case 1: (Uy == 0) and (Uz == 0)\n")

    %1D case:
    [Ux0,cond,iter] = ridders(@f_case1,uRx, uLx, tol, a, Ay, dp, By, Ax, Bx, Az, Bz, c);

    %Convert Back to Ux Uy Uz:
    Ux = Ux0;
    Uy = 0;
    Uz = 0;

elseif uLy == 0.0 && uRy == 0.0  && a ~= 0.0

    % Print which case:
    %fprintf("Case 2: (Uy == 0)\n")

    %2D case (uy == 0):
    [Uz0,cond,iter] = ridders(@f_case2,uRz, uLz, tol, a, Ay, dp, By, Ax, Bx, Az, Bz, c);

    %Convert Back to Ux Uy Uz:
    Vx = (-Uz0*a - Az)/(Uz0*dp + Bz);
    Ux = (Ax + Bx*Vx)/(-dp*Vx - a);
    Uy = 0;
    Uz = Uz0;

elseif uLz == 0.0 && uRz == 0.0  && a ~= 0.0

    % Print which case:
    %fprintf("Case 3: (Uz == 0)\n")

    %Full 2D case, Uz = 0;
    [Uy0,cond,iter] = ridders(@f_case3,uRy, uLy, tol, a, Ay, dp, By, Ax, Bx, Az, Bz, c);

    %Convert Back to Ux Uy Uz:
    Vx = (-Uy0*a - Ay)/(Uy0*dp + By);
    Ux = (Ax + Bx*Vx)/(-dp*Vx - a);
    Uy = Uy0;
    Uz = 0;

elseif a ~= 0.0 && By ~= 0.0

    % Print which case:
    %fprintf("Case 4: No Zeros \n")
    %diagnostic_on_direct_solve(a, c, dp, Ax, Ay, Az, Bx, By, Bz, Uy0, uRy, uLy)

    %Full 3D case
    [Uy0,cond,iter] = ridders(@f_case4,uRy, uLy, tol, a, Ay, dp, By, Ax, Bx, Az, Bz, c);

    %Convert Back to Ux Uy Uz:
    Vx = (-Uy0*a - Ay)/(Uy0*dp + By);
    Ux = (Ax + Bx*Vx)/(-dp*Vx - a);
    Uy = Uy0;
    Uz = (Az + Bz*Vx)/(-dp*Vx - a);

else
    % If not the right criterion, use simpler solves:
    fprintf("Undefined Cases!\n")
    Ux = 0;
    Uy = 0;
    Uz = 0;
    cond = 0;
end


%Compute the exact values to hand back
gamma = sqrt(1 + (Ux^2 + Uy^2 + Uz^2)/c^2);
Vx = Ux/gamma;
Vy = Uy/gamma;
Vz = Uz/gamma;
v0 = [Vx,Vy,Vz];
%fprintf("(Ridders) Converged in %d iterations\n",iter)

end


% Find the root
function [sol,cond,iter] = ridders(f,bound1, bound2, tol, a, Ay, dp, By, Ax, Bx, Az, Bz, c)

%Select upper and lower
if bound1 > bound2
    upper = bound1;
    lower = bound2;
else
    upper = bound2;
    lower = bound1;
end

%Setup the variables
x0 = lower;
x1 = 0.5*(lower + upper);
x2 = upper;
d = x2 - x1;
F0 = f(x0,a,Ay,dp,By,Ax,Bx,Az,Bz,c);
F2 = f(x2,a,Ay,dp,By,Ax,Bx,Az,Bz,c);
F3 = 2*tol;
iter = 0;

% condition 0 -> Okay, 1 -> Fail
cond = 0;

while( cond == 0 && iter < 15)  %abs(F3) > tol &&

    %Fails, multiple roots/no root etc, bad domain
    if F0*F2 >= 0
        cond = 1;
        %fprintf("Ridders' algo. fails to have two signed ends!\n")
        sol = 0;

        if d > tol && tol > 1e-40
        % Create a figure in this case:
        %figure
        N = 100;
        x = linspace(lower,upper,N);
        f_x = zeros(1,N);
        for i = 1:N
            f_x(i) = f(x(i),a,Ay,dp,By,Ax,Bx,Az,Bz,c);
        end
        plot(x,f_x)
        title("Ridders Diag, x, f_x")
        xlabel("U")
        ylabel("f(U)")
        ylim([-max(abs(f_x)),max(abs(f_x))])

        else % fails but with condition 2, small domain might not be valid
            cond = 2;
        end


        % Otherwise continue
    else

        %Compute more variables
        F1 = f(x1,a,Ay,dp,By,Ax,Bx,Az,Bz,c);
        %W = F1^2 - F0*F2;
        %m = (1/d)*ln( (F1 - sign(F0)*sqrt(W))/F2 );
        %H1 = H(x1);
        %H2 = H(x2);
        diff = d*(F1/F0)/sqrt( (F1/F0)^2 - F2/F0 );
        x3 = x1 + diff;

        % Setup the next interval
        F3 = f(x3,a,Ay,dp,By,Ax,Bx,Az,Bz,c);

        % Solution
        if x2 >= x3 && x3 >= x1
            if F1*F3 < 0
                upper = x3;
                lower = x1;
            elseif F2*F3 < 0
                upper = x2;
                lower = x3;
            else
                done = 1;
            end
        elseif x0 <= x3 && x3 <= x1
            if F1*F3 < 0
                upper = x1;
                lower = x3;
            elseif F0*F3 < 0
                upper = x3;
                lower = x0;
            else
                done = 1;
            end
        else
            done = 1;
        end

        % Save the new bounds, and new solution
        x2 = upper;
        x0 = lower;
        x1 = 0.5*(lower + upper);
        d = x2 - x1;
        F0 = f(x0,a,Ay,dp,By,Ax,Bx,Az,Bz,c);
        F2 = f(x2,a,Ay,dp,By,Ax,Bx,Az,Bz,c);
        sol = x3;
        iter = iter + 1;

        %fprintf("Iteration: %d, returns val: %1.16f, f(x3): %1.16f\n",iter,sol,F3)

    end

end


end

function [eval] = f_case4(Uy0,a,Ay,dp,By,Ax,Bx,Az,Bz,c)

% Compute values with old u0
Vx = (-Uy0*a - Ay)/(Uy0*dp + By);
Ux = (Ax + Bx*Vx)/(-dp*Vx - a);
Uz = (Az + Bz*Vx)/(-dp*Vx - a);

% Eval of the function
eval = Vx -  Ux / sqrt(1 + (Ux^2 + Uy0^2 + Uz^2)/c^2);
end


function [eval] = f_case3(Uy0,a,Ay,dp,By,Ax,Bx,~,~,c)
% Compute values with old u0
Vx = (-Uy0*a - Ay)/(Uy0*dp + By);
Ux = (Ax + Bx*Vx)/(-dp*Vx - a);

% Newton Iterate: (From CAS/F_x_prime_Uy_for_3D_rel_solve.m)
eval = Vx -  Ux / sqrt(1 + (Ux^2 + Uy0^2)/c^2);
end

function [eval] = f_case2(Uz0,a,~,dp,~,Ax,Bx,Az,Bz,c)
% Compute values with old u0
Vx = (-Uz0*a - Az)/(Uz0*dp + Bz);
Ux = (Ax + Bx*Vx)/(-dp*Vx - a);

% Newton Iterate: (From CAS/F_x_prime_Uy_for_3D_rel_solve.m)
eval = Vx -  Ux / sqrt(1 + (Ux^2 + Uz0^2)/c^2);
end

function [eval] = f_case1(Ux0,a,~,dp,~,Ax,Bx,~,~,c)
% Compute values with old u0
Vx = (-Ux0*a - Ax)/(Ux0*dp + Bx);

% Newton Iterate: (From CAS/F_x_prime_Uy_for_3D_rel_solve.m)
eval = Vx -  Ux0 / sqrt(1 + (Ux0^2)/c^2);
end