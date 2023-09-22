%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Name: Grant Johnson
%Date: 5/17/2023
%Fluid algroithm, V (divergence U)

%Notes:
%-1D
% Fluid Scheme: https://ammar-hakim.org/sj/hancock-muscl.html
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Update the quanitites Ux, Uy, Uz (t -> t + dt)
function [Ux,Uy,Uz,Vx,Vy,Vz,N,grid] = fluid_grad_U_prims(Ux,Uy,Uz,N,grid)

% Build Q
Q_prim = construct(N, Ux, Uy, Uz);
Nx = grid.Nx;

% TVD Diagnostic
%[grid] = TVD_diagnostic("start",Q(1,:),Q(2,:),Q(3,:),Q(4,:),grid);

%Iterate over the domain
% I = linspace(1,Nx-1,Nx-1); %DEFAULT
R = mod( linspace(1,Nx,Nx), Nx-1) + 1; %mod( linspace(1,Nx,Nx), Nx) + 1; %Good
L = mod( linspace(-1,Nx-2,Nx), Nx-1) + 1; %mod( linspace(-1,Nx-2,Nx), Nx) + 1; %Good

%Reconstruct U(x) (soln) within one cell, also slope dU
dQ = reconstruct(Q_prim,Q_prim(:,L),Q_prim(:,R));

%Update solution, primative variables
%Q_tilde = Q - grid.dt/(2*grid.dx)*AQ(Q,grid).*dQ;
A = AQ(Q_prim,grid);
Q_tilde = zeros(4,grid.Nx);
for i = 1:Nx
    Q_tilde(:,i) = Q_prim(:,i) - grid.dt/(2*grid.dx)*A(:,:,i)*dQ(:,i);
end

%Get edge values (edges_linear)
% Output: (cell edge values) in cell i
% [ i - 1/2 (+), i + 1/2 (-) ]
% Positivity Limiter:
[Q_plus_I, Q_minus_I] = edges_linear(Q_tilde,dQ);
limited_cases = 0;
% Limit for positivity
for i = 1:grid.Nx
    if (Q_plus_I(1,i)<0) || (Q_minus_I(1,i)<0)
        dQ(1,i) = 0.0;
        dQ(2,i) = 0.0;
        dQ(3,i) = 0.0;
        dQ(4,i) = 0.0;
        limited_cases = limited_cases + 1;
    end
end
%fprintf("Limited Cases per iteration: %d\n",limited_cases);
% Recompute Q-tilde
for i = 1:Nx
    Q_tilde(:,i) = Q_prim(:,i) - grid.dt/(2*grid.dx)*A(:,:,i)*dQ(:,i);
end
[Q_plus_I, Q_minus_I] = edges_linear(Q_tilde,dQ);
[Q_plus_R, ~] = edges_linear(Q_tilde(:,R),dQ(:,R));
[~, Q_minus_L] = edges_linear(Q_tilde(:,L),dQ(:,L));

%Update the fluxes
check(grid,A,"A_main");
check(grid,dQ,"dQ_main");
check(grid,Q_prim,"Q_main");
check(grid,Q_tilde,"Q_tilde");
F_R =  Flux(Q_minus_I,Q_plus_R,grid);
F_L = Flux(Q_minus_L,Q_plus_I,grid);

%Stop if isnan fails
if max(max((isnan(F_R))))
    str = "F_R";
    fprintf("FAILS: Nan Soln.\n");
    fprintf("STR: %s\n",str);
    pause(1000)
end

%Compute the updated Q
Q = construct_cons(N, Ux, Uy, Uz);
Q = Q - grid.dt/(grid.dx)*(F_R - F_L);

%Stop if isnan fails
if max(max((isnan(Q_prim))))
    str = "Q";
    fprintf("FAILS: Nan Soln.\n");
    fprintf("STR: %s\n",str);
    pause(1000)
end

check(grid,Q_prim,"Q");
check(grid,F_R,"F_R");
check(grid,F_L,"F_L");

% TVD Diagnostic
%[grid] = TVD_diagnostic("end",Q(1,:),Q(2,:),Q(3,:),Q(4,:),grid);

% Destruct Q into it's components
[N, Ux, Uy, Uz] = destruct_cons(grid,Q);

%Save the output
gamma = sqrt(1+(Ux.*Ux + Uy.*Uy + Uz.*Uz)/(grid.c*grid.c));
Vx = Ux./gamma;
Vy = Uy./gamma;
Vz = Uz./gamma;

check(grid,Vx,"Vx");
check(grid,Vy,"Vy");
check(grid,Vz,"Vz");

end


%Locally defined functions for MUSCL-Handcock:


%Construct Q
function [Q] = construct(N, Ux, Uy, Uz)
Q = [N; Ux; Uy; Uz];
end


%Destruct Q into N, Ux, Uy, Uz
function [N, Ux, Uy, Uz] = destruct(Q)
N = Q(1,:);
Ux = Q(2,:);
Uy = Q(3,:);
Uz = Q(4,:);
end

%Construct Q
function [Q] = construct_cons(N, Ux, Uy, Uz)
Q = [N', (N.*Ux)', (N.*Uy)', (N.*Uz)'];
Q = Q';
end

%Destruct Q into N, Ux, Uy, Uz
function [N, Ux, Uy, Uz] = destruct_cons(grid,Q)
N = Q(1,:);
Ux = Q(2,:)./N;
Uy = Q(3,:)./N;
Uz = Q(4,:)./N;
Ux(N==0) = 0;
Uy(N==0) = 0;
Uz(N==0) = 0;
check(grid,N,"N");
check(grid,Ux,"Ux");
check(grid,Uy,"Uy");
check(grid,Uz,"Uz");
end

%Fluxes N
function [Flux_vec] = Flux(QL,QR,grid)

c = grid.c;
[NR, UxR, UyR, UzR] = destruct(QR);
[NL, UxL, UyL, UzL] = destruct(QL);

% compute c for the lax flux
UR2 = ( UxR.^2 + UyR.^2 + UzR.^2 )/(c^2);
UL2 = ( UxL.^2 + UyL.^2 + UzL.^2 )/(c^2);
vRx = UxR./(sqrt(1 + UR2));
vLx = UxL./(sqrt(1 + UL2));
c = max( abs(vRx), abs(vLx) ); %%% FIX ABS | lamdba^p|

%Rusanov Flux F = vxQ
FR = [ NR.*vRx  ; NR.*vRx.*UxR ; NR.*vRx.*UyR  ; NR.*vRx.*UzR ];
FL = [ NL.*vLx  ; NL.*vLx.*UxL ; NL.*vLx.*UyL  ; NL.*vLx.*UzL ];
QR_vec = [ NR  ; NR.*UxR ; NR.*UyR  ; NR.*UzR ];
QL_vec = [ NL  ; NL.*UxL ; NL.*UyL  ; NL.*UzL ];
Flux_vec = (1/2) * ( FR + FL  - c.*( QR_vec - QL_vec ) );


%Stop if isfinite fails
if min(min((isfinite(NR.*UxR)))) == 0
    str = "NR.*UxR";
    fprintf("FAILS: Inf Soln.\n");
    fprintf("STR: %s\n",str);
    pause(1000)
end

check(grid,NR,"NR")
check(grid,UxR,"UxR")
check(grid,UyR,"UyR")
check(grid,UzR,"UzR")

check(grid,NR,"NR")
check(grid,NR.*UxR,"NUxR")
check(grid,NR.*UyR,"NUyR")
check(grid,NR.*UzR,"NUzR")

check(grid,FR,"FR")
check(grid,FL,"FL")
check(grid,QR_vec,"QR_vec")
check(grid,QL_vec,"QL_vec")
check(grid,Flux_vec,"Flux_vec")

%Stop if isnan fails
if max(max((isnan(FR))))
    str = "FR";
    fprintf("FAILS: Nan Soln.\n");
    fprintf("STR: %s\n",str);
    pause(1000)
end

%Stop if isnan fails
if max(max((isnan(QR))))
    str = "QR";
    fprintf("FAILS: Nan Soln.\n");
    fprintf("STR: %s\n",str);
    pause(1000)
end

%Stop if isnan fails
if max(max((isnan(FL))))
    str = "FL";
    fprintf("FAILS: Nan Soln.\n");
    fprintf("STR: %s\n",str);
    pause(1000)
end

%Stop if isnan fails
if max(max((isnan(QL))))
    str = "QL";
    fprintf("FAILS: Nan Soln.\n");
    fprintf("STR: %s\n",str);
    pause(1000)
end

%Stop if isnan fails
if max(max((isnan(c))))
    str = "c";
    fprintf("FAILS: Nan Soln.\n");
    fprintf("STR: %s\n",str);
    pause(1000)
end

%Stop if isnan fails
if max(max((isnan(Flux_vec))))
    str = "Flux_vec";
    fprintf("FAILS: Nan Soln.\n");
    fprintf("STR: %s\n",str);
    pause(1000)
end

end

%Reconstruction
function [dW] = reconstruct(Wi,Wm,Wp)

%Average Dw
dW = ave( Wi - Wm, Wp - Wi );

end

% Averaging
function [dW] = ave( Wm, Wp )

avg_type = "minmod"; %"zero"; %"minmod_more_diff"; %"minmod"; %"Supebee"; %"standard";
a = Wm; b = Wp;
ab = a.*b;
sz_a = size(a);
dW = zeros(sz_a);

% Standard Averaging
if avg_type == "standard"
    dW = (Wm + Wp)/2;
elseif avg_type == "zero"
    dW = zeros(sz_a);
elseif avg_type == "minmod_more_diff"
    for i = 1:sz_a(1)
        for j = 1:sz_a(2)
            if ab(i,j) > 0
                dW(i,j) = minmod( [(a(i,j) + b(i,j))/2 ,a(i,j), b(i,j)] ); %minmod( [(a(i,j) + b(i,j))/2 , a(i,j), b(i,j)] );
            else
                dW(i,j) = 0;
            end
        end
    end
elseif avg_type == "minmod"
    for i = 1:sz_a(1)
        for j = 1:sz_a(2)
            if ab(i,j) > 0
                dW(i,j) = minmod( [(a(i,j) + b(i,j))/2 , 2*a(i,j), 2*b(i,j)] );
            else
                dW(i,j) = 0;
            end
        end
    end
elseif avg_type == "Supebee"
    for i = 1:sz_a(1)
        for j = 1:sz_a(2)
            if ab(i,j) > 0
                max_v =  maxmod(  [ a(i,j)  , b(i,j)  ] );
                min_v = minmod (  [ 2*a(i,j), 2*b(i,j)] );
                dW(i,j) = minmod( [ max_v   , min_v   ] );
            else
                dW(i,j) = 0;
            end
        end
    end
end


end

% Minmod used for averaging
function [val] = minmod(a)
if (max(a) > 0) && (min(a) >  0)
    val = min(a);
elseif (max(a) < 0) && (min(a) <  0)
    val = max(a);
else
    val = 0;
end
end

%Maxmod used for averaging
function [val] = maxmod(a)
if (max(a) > 0) && (min(a) >  0)
    val = max(a);
elseif (max(a) < 0) && (min(a) <  0)
    val = min(a);
else
    val = 0;
end
end


% Function A(Q)
function [A] = AQ(Q,grid)

% Build the 4x4 flux Jacobian... (See Maxima generated values)
% [ A11, A12, A13, A14 ]
% [ A21, A22, A23, A24 ]
% [ A31, A32, A33, A34 ]
% [ A41, A42, A43, A44 ]

%Recover primitive variables
[N, Ux, Uy, Uz] = destruct(Q);

%Definte helper functions
c = grid.c;
c_sq = c^2;
gamma = sqrt(1.0 + (Ux.*Ux + Uy.*Uy + Uz.*Uz)/(c_sq) );
Vx = Ux./gamma;

%Assemble A
A = zeros(4,4,grid.Nx);
for i = 1:grid.Nx
    A01 = (N(i)/gamma(i))*(1 - Vx(i)*Vx(i)/(c*c));
    A02 = -N(i)*Ux(i)*Uy(i)/(c*c*gamma(i)*gamma(i)*gamma(i));
    A03 = -N(i)*Ux(i)*Uz(i)/(c*c*gamma(i)*gamma(i)*gamma(i));
    A(:,:,i) = [[ Vx(i), A01, A02, A03 ];...
                [ 0, Vx(i), 0, 0 ];...
                [ 0, 0, Vx(i), 0 ];...
                [ 0, 0, 0, Vx(i) ] ];
end

end


function [W_plus, W_minus] = edges_linear(W_tilde,dW)
W_plus = W_tilde - dW/2;
W_minus = W_tilde + dW/2;
end


%Check that X is real and isnotnan
function check(grid,X,str)

%Stop if isreal fails
if (~isreal(X))
    fprintf("FAILS: Not-Real Soln.\n");
    fprintf("STR: %s\n",str);
    pause(1000)
end

%Stop if isnan fails
if max(max((isnan(X))))
    fprintf("FAILS: Nan Soln.\n");
    fprintf("STR: %s\n",str);
    pause(1000)
end

%Stop if isfinite fails
if min(min((isfinite(X)))) == 0
    fprintf("FAILS: Inf Soln.\n");
    fprintf("STR: %s\n",str);
    pause(1000)
end

end
