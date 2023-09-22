%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Name: Grant Johnson
%Date: 5/17/2023
%Fluid algroithm, V (divergence U)

%Notes:
%-1D
% Fluid Scheme: https://ammar-hakim.org/sj/hancock-muscl.html
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Update the quanitites Ux, Uy, Uz (t -> t + dt)
function [F_L,F_R] = FH(Ux,Uy,Uz,N,grid)

% Build Q
Q = construct(N, Ux, Uy, Uz);
Nx = grid.Nx;

% TVD Diagnostic
%[grid] = TVD_diagnostic("start",Q(1,:),Q(2,:),Q(3,:),Q(4,:),grid);

%Iterate over the domain
% I = linspace(1,Nx-1,Nx-1); %DEFAULT
R = mod( linspace(1,Nx,Nx), Nx-1) + 1; %mod( linspace(1,Nx,Nx), Nx) + 1; %Good
L = mod( linspace(-1,Nx-2,Nx), Nx-1) + 1; %mod( linspace(-1,Nx-2,Nx), Nx) + 1; %Good

%Reconstruct U(x) (soln) within one cell, also slope dU
dQ = reconstruct(Q,Q(:,L),Q(:,R));

%Update solution, primative variables
%Q_tilde = Q - grid.dt/(2*grid.dx)*AQ(Q,grid).*dQ;
A = AQ(Q,grid);
Q_tilde = zeros(4,grid.Nx);
for i = 1:Nx
    Q_tilde(:,i) = Q(:,i) - grid.dt/(2*grid.dx)*A(:,:,i)*dQ(:,i);
end

%Get edge values (edges_linear)
% Output: (cell edge values) in cell i
% [ i - 1/2 (+), i + 1/2 (-) ]
% Positivity Limiter:
[Q_plus_I, Q_minus_I] = edges_linear(Q_tilde,dQ);
limited_cases = 0;
% Limit for positivity
for i = 1:grid.Nx
    if (Q_plus_I(1,i)<0) || (Q_minus_I(1,i)<0) || ...
        sign(Q_tilde(2,i)) ~= sign(Q_plus_I(2,i)) || sign(Q_tilde(2,i)) ~= sign(Q_minus_I(2,i)) ||...
        sign(Q_tilde(3,i)) ~= sign(Q_plus_I(3,i)) || sign(Q_tilde(3,i)) ~= sign(Q_minus_I(3,i)) ||...
        sign(Q_tilde(4,i)) ~= sign(Q_plus_I(4,i)) || sign(Q_tilde(4,i)) ~= sign(Q_minus_I(4,i))
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
    Q_tilde(:,i) = Q(:,i) - grid.dt/(2*grid.dx)*A(:,:,i)*dQ(:,i);
end
[Q_plus_I, Q_minus_I] = edges_linear(Q_tilde,dQ);
[Q_plus_R, ~] = edges_linear(Q_tilde(:,R),dQ(:,R));
[~, Q_minus_L] = edges_linear(Q_tilde(:,L),dQ(:,L));

%Update the fluxes
F_R =  Flux(Q_minus_I,Q_plus_R,grid);
F_L = Flux(Q_minus_L,Q_plus_I,grid);

end


%Locally defined functions for MUSCL-Handcock:


%Construct Q
function [Q] = construct(N, Ux, Uy, Uz)
Q = [N', (N.*Ux)', (N.*Uy)', (N.*Uz)'];
Q = Q';
end


%Destruct Q into N, Ux, Uy, Uz
function [N, Ux, Uy, Uz] = destruct(Q)
N = Q(1,:);
Ux = Q(2,:)./N;
Uy = Q(3,:)./N;
Uz = Q(4,:)./N;
end

%Fluxes N
function [Fl] = Flux(QL,QR,grid)

c = grid.c;
[~, UxR, UyR, UzR] = destruct(QR);
[~, UxL, UyL, UzL] = destruct(QL);

% compute c for the lax flux
UR2 = ( UxR.^2 + UyR.^2 + UzR.^2 )/(c^2);
UL2 = ( UxL.^2 + UyL.^2 + UzL.^2 )/(c^2);
vRx = UxR./(sqrt(1 + UR2));
vLx = UxL./(sqrt(1 + UL2));
c = max( abs(vRx), abs(vLx) ); %%% FIX ABS | lamdba^p|

%Rusanov Flux F = vxQ
FL = vLx.*QL;
FR = vRx.*QR;
Fl = (1/2) * ( FR + FL  - c.*( QR - QL ) );

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
[~, Ux, Uy, Uz] = destruct(Q);

%Definte helper functions
c = grid.c;
c_sq = c^2;
gamma = sqrt(1.0 + (Ux.*Ux + Uy.*Uy + Uz.*Uz)/(c_sq) );
a = c_sq*gamma.*gamma.*gamma;

%Build the components of A
A00 = (Ux.*(Uz.*Uz)+Ux.*(Uy.*Uy)+(Ux.*Ux.*Ux))./a;
A01 = ((c_sq)+(Uz.*Uz)+(Uy.*Uy))./a;
A02 = -(Ux.*Uy)./a;
A03 = -(Ux.*Uz)./a;
A10 = -(Ux.*Ux)./(gamma.*gamma.*gamma);
A11 = (2.0*Ux*(c_sq)+2.0*Ux.*(Uz.*Uz)+2.0*Ux.*(Uy.*Uy)+(Ux.*Ux.*Ux))./a;
A12 = -((Ux.*Ux).*Uy)./a;
A13 = -((Ux.*Ux).*Uz)./a;
A20 = -(Ux.*Uy)./(gamma.*gamma.*gamma);
A21 = (Uy*(c_sq)+Uy.*(Uz.*Uz)+(Uy.*Uy.*Uy))./a;
A22 = (Ux*(c_sq)+Ux.*(Uz.*Uz)+(Ux.*Ux.*Ux))./a;
A23 = -(Ux.*Uy.*Uz)./a;
A30 = -(Ux.*Uz)./(gamma.*gamma.*gamma);
A31 = (Uz.*(c_sq)+(Uz.*Uz.*Uz)+(Uy.*Uy).*Uz)./a;
A32 = -(Ux.*Uy.*Uz)./a;
A33 = (Ux*(c_sq)+Ux.*(Uy.*Uy)+(Ux.*Ux.*Ux))./a;

%Assemble A
A = zeros(4,4,grid.Nx);
for i = 1:grid.Nx
A(:,:,i) = [ [ A00(i), A01(i), A02(i), A03(i) ];...
             [ A10(i), A11(i), A12(i), A13(i) ];...
             [ A20(i), A21(i), A22(i), A23(i) ];...
             [ A30(i), A31(i), A32(i), A33(i) ] ];
end

end


function [W_plus, W_minus] = edges_linear(W_tilde,dW)
W_plus = W_tilde - dW/2;
W_minus = W_tilde + dW/2;
end

