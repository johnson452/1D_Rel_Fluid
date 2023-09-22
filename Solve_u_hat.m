%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grant Johnson
% Solves for the shock speed for various equaitons
%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all


%%%% (NON-RELATIVSITIC) %%%%
% Make the symbolics 
syms rhoL rhoR vL vR t vs
assume(rhoR >= 0 )
assume(rhoL >= 0 )
assume(t >= 0 )
assume(vL > vR )
assume(vL > vs)
assume(vs > vR)

% Left and right masses accreted
mL = rhoL*(vL - vs)*t;
mR = rhoR*(vs - vR)*t;

% Conservation of momentum
eqn1 = (mL + mR)*vs == mL*vL + mR*vR;

% Solve
eqns = eqn1;

[S,parameters,conditions] = solve(eqns,vs,'ReturnConditions',true);

% Select the right function
disp(S)
