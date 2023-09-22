%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grant Johnson
% Solves for the shock speed for various equaitons
%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all


%%%% (RELATIVSITIC 3D) %%%%
syms  vsx vsy vsz
c = 1.0;
rhoL = 1.0;
rhoR = 2.0;
uLx = 2.0;
uRx = 1.0;
uLy = 2.0;
uRy = 1.0;
uLz = 2.0;
uRz = 1.0;
vLx = uLx/(sqrt(1 + (uLx*uLx + uLy*uLy + uLz*uLz )/(c*c)));
vRx = uRx/(sqrt(1 + (uRx*uRx + uRy*uRy + uRz*uRz )/(c*c)));


% Left and right masses accreted
mL = rhoL*(vLx - vsx);
mR = rhoR*(vsx - vRx);

% Conservation of momentum
eqn1 = (mL + mR)*vsx/sqrt(1.0 - (vsx*vsx + vsy*vsy + vsz*vsz)/(c*c)) == mL*uLx + mR*uRx;
eqn2 = (mL + mR)*vsy/sqrt(1.0 - (vsx*vsx + vsy*vsy + vsz*vsz)/(c*c)) == mL*uLy + mR*uRy;
eqn3 = (mL + mR)*vsz/sqrt(1.0 - (vsx*vsx + vsy*vsy + vsz*vsz)/(c*c)) == mL*uLz + mR*uRz;

% Solve
eqns = [eqn1, eqn2, eqn3];

%[S,parameters,conditions] = 
S = solve(eqns,[vsx, vsy, vsz],'ReturnConditions',true);

% Select the right function
disp(simplify(S.vsx))
