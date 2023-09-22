%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grant Johnson
% Solves for the shock speed for various equaitons
%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all


%%%% (RELATIVSITIC 1D) %%%%
% Make the symbolics 
syms rhoL rhoR vL vR t vs c
assume(rhoR >= 0 )
assume(rhoL >= 0 )
assume(t > 0 )
assume(vL > vR )
assume(vL > vs)
assume(vs > vR)
assume(c > 0)
assume(vs < c)
assume(vL < c)
assume(vR < c)
assume(-c < vs)
assume(-c < vL)
assume(-c < vR)
% uL = vL/sqrt(1.0 - (vL*vL)/(c*c));
% uR = vR/sqrt(1.0 - (vR*vR)/(c*c));
% us = vs/sqrt(1.0 - (vs*vs)/(c*c));
% assume(us,'real')
% assume(uL,'real')
% assume(uR,'real')

% Left and right masses accreted
mL = rhoL*(vL - vs)*t;
mR = rhoR*(vs - vR)*t;

% Conservation of momentum
eqn1 = (mL + mR)*vs/sqrt(1.0 - (vs*vs)/(c*c)) == mL*vL/sqrt(1.0 - (vL*vL)/(c*c)) + mR*vR/sqrt(1.0 - (vR*vR)/(c*c));

% Solve
eqns = eqn1;

%[S,parameters,conditions] = 
S = solve(eqns,vs,'ReturnConditions',true);

% Select the right function
disp(S)
