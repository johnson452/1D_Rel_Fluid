%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Name: Grant Johnson
%Date: 5/17/2023
%Fluid algroithm, V (divergence U)

%Notes:
%-1D
% Fluid Scheme: https://ammar-hakim.org/sj/hancock-muscl.html
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Update the quanitites Ux, Uy, Uz (t -> t + dt)
function [Ux,Uy,Uz,N,grid] = fluid_grad_U(Ux,Uy,Uz,N,grid)

%Locally interpolate U, N to cell centers for update
%Uy = interp_edge_to_center(Uy,grid);
%Uz = interp_edge_to_center(Uz,grid);
%N = interp_edge_to_center(N,grid);

% Build Q
Q = construct(N, Ux, Uy, Uz);
Nx = grid.Nx;


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
[Q_plus_I, Q_minus_I] = edges_linear(Q_tilde,dQ);
[Q_plus_R, ~] = edges_linear(Q_tilde(:,R),dQ(:,R));
[~, Q_minus_L] = edges_linear(Q_tilde(:,L),dQ(:,L));

%Get the edge values (cell edges, in cell i)
% [Q_plus_I, Q_minus_I] = edges(Q_tilde);
% [Q_plus_R, ~] = edges(Q_tilde(:,R));
% [~, Q_minus_L] = edges(Q_tilde(:,L));

%Update the fluxes
F_R =  Flux(Q_minus_I,Q_plus_R);
F_L = Flux(Q_minus_L,Q_plus_I);

%Compute the updated Q
Q = Q - grid.dt/(grid.dx)*(F_R - F_L);

% Destruct Q into it's components
[N, Ux, Uy, Uz] = destruct(Q);

% Interpolate back
%Uy = interp_center_to_edge(Uy,grid);
%Uz = interp_center_to_edge(Uz,grid);
%N = interp_center_to_edge(N,grid);

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
function [Fl] = Flux(QL,QR)

[~, UxR, UyR, UzR] = destruct(QR);
[~, UxL, UyL, UzL] = destruct(QL);

% compute c for the lax flux
UR2 = ( UxR.^2 + UyR.^2 + UzR.^2 );
UL2 = ( UxL.^2 + UyL.^2 + UzL.^2 );
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

avg_type = "minmod"; %"Supebee"; %"standard";
a = Wm; b = Wp;
ab = a.*b;
sz_a = size(a);
dW = zeros(sz_a);

% Standard Averaging
if avg_type == "standard"
    dW = (Wm + Wp)/2;
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
a = Uz.^4+(2*Uy.^2+2*Ux.^2+2).*Uz.^2+Uy.^4+(2*Ux.^2+2).*Uy.^2+Ux.^4+2*Ux.^2+1;
gamma = sqrt(1 + Ux.*Ux + Uy.*Uy + Uz.*Uz );

%Build the components of A
A11 = ((Ux.*Uz.^2+Ux.*Uy.^2+Ux.^3).*gamma)./a;
A12 = ((Uz.^2+Uy.^2+1).*gamma)./a;
A13 = -(Ux.*Uy)./gamma.^3;
A14 = -(Ux.*Uz)./gamma.^3;
A21 = -Ux.^2./gamma.^3;
A22 = ((2.*Ux.*Uz.^2+2.*Ux.*Uy.^2+Ux.^3+2.*Ux).*gamma)./a;
A23 = -(Ux.^2.*Uy)./gamma.^3;
A24 = -(Ux.^2.*Uz)./gamma.^3;
A31 = -(Ux.*Uy)./gamma.^3;
A32 = ((Uy.*Uz.^2+Uy.^3+Uy).*gamma)./a;
A33 = ((Ux.*Uz.^2+Ux.^3+Ux).*gamma)./a;
A34 = -(Ux.*Uy.*Uz)./gamma.^3;
A41 = -(Ux.*Uz)./gamma.^3;
A42 = ((Uz.^3+(Uy.^2+1).*Uz).*gamma)./a;
A43 = -(Ux.*Uy.*Uz)./gamma.^3;
A44 = ((Ux.*Uy.^2+Ux.^3+Ux).*gamma)./a;

%Assemble A
A = zeros(4,4,grid.Nx);
for i = 1:grid.Nx
A(:,:,i) = [ [ A11(i), A12(i), A13(i), A14(i) ];...
             [ A21(i), A22(i), A23(i), A24(i) ];...
             [ A31(i), A32(i), A33(i), A34(i) ];...
             [ A41(i), A42(i), A43(i), A44(i) ] ];
end

end


function [W_plus, W_minus] = edges_linear(W_tilde,dW)
W_plus = W_tilde - dW/2;
W_minus = W_tilde + dW/2;
end


% Compute the edge values
function [W_plus, W_minus] = edges(W_tilde)
% |+  x  -|
%Linear assumption

%Compute W_plus and W_minus
% (Linear)
% W_plus = W_tilde - dW/2;
% W_minus = W_tilde + dW/2;

%Interpolate 1:
%W_tilde = W_tilde';
sz = size(W_tilde);
sz2 = [sz(1),sz(2)+1];
W_tilde_interp = zeros(sz2);
W_plus = zeros(sz);
W_minus = zeros(sz);
for i = 1:sz2(1)
    W_tilde_interp(i,:) = interp_center_to_edge_local(W_tilde(i,:));
    W_plus(i,:) = W_tilde_interp(i,1:end-1);
    W_minus(i,:) = W_tilde_interp(i,2:end);
end


end

% Local center to edge, with unique fluxes
function [y_interp] = interp_center_to_edge_local(y)
Nx = max(size(y));
x = linspace(0,1,Nx);
dx = x(2)-x(1);
x2 = linspace(0+dx/2,1-dx/2,Nx-1);
y_interp = interp1(x,y,x2,'spline');
y_interp = [y_interp(1),y_interp,y_interp(end)]; % Signal issue if not delt with!
end

%Sup function
function [sup_val] = sup(values)
[maxB, index] = max(abs(values));
sup_val = maxB * sign(values(index));
end