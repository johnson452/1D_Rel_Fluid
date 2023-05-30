%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Name: Grant Johnson
%Date: 5/17/2022
%Fluid algroithm, V (divergence U)

%Notes:
%-1D
% Fluid Scheme: https://ammar-hakim.org/sj/hancock-muscl.html
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Update the quanitites Ux, Uy, Uz (t -> t + dt)
function [Ux,Uy,Uz,N,grid] = fluid_grad_U(Ux,Uy,Uz,N,grid)

%Locally interpolate U, N to cell centers for update
Uy = interp_edge_to_center(Uy,grid);
Uz = interp_edge_to_center(Uz,grid);
N = interp_edge_to_center(N,grid);


% Build U
U = [Ux', Uy', Uz'];
N = N';
Nx = grid.Nx;

%Iterate over the domain
% I = linspace(1,Nx-1,Nx-1); %DEFAULT
R = mod( linspace(1,Nx-1,Nx-1), Nx-1) + 1;
L = mod( linspace(-1,Nx-3,Nx-1), Nx-1) + 1;
% R2 = mod( linspace(1,Nx,Nx), Nx) + 1;
% L2 = mod( linspace(-1,Nx-2,Nx), Nx) + 1;

%Reconstruct U(x) (soln) within one cell, also slope dU
[~,dU] = reconstruct(U,U(L,:),U(R,:),0,grid.x2,grid.dx);
[~,dN] = reconstruct(N,N(L),N(R),0,grid.x2,grid.dx);

%Update solution, primative variables
U_tilde = U - grid.dt/(2*grid.dx)*Aprim(U,N,"fluid").*dU;
N_tilde = N - grid.dt/(2*grid.dx)*AQ(U,N,"density").*dN;


%Build Q: (it's a vector) [NUx, NUy, NUz]
Q_NU_tilde = construct_UN(N,U_tilde);
Q_NU = construct_UN(N,U);
%No corresponding N updates since Q = N
Q_N_tilde = N_tilde;
Q_N = N;

%Get the edge values (cell edges, in cell i)
[Q_NU_plus_I, Q_NU_minus_I] = edges(Q_NU_tilde,grid);
[~, Q_NU_minus_R] = edges(Q_NU_tilde(R,:),grid);
[Q_NU_plus_L, ~] = edges(Q_NU_tilde(L,:),grid);

[Q_N_plus_I, Q_N_minus_I] = edges(Q_N_tilde,grid);
[~, Q_N_minus_R] = edges(Q_N_tilde(R),grid);
[Q_N_plus_L, ~] = edges(Q_N_tilde(L),grid);

%Update the fluxes
F_NU_plus =  F_NU(Q_NU_plus_I,Q_NU_minus_R, Q_N_plus_I, Q_N_minus_R);
F_NU_minus = F_NU(Q_NU_plus_L,Q_NU_minus_I, Q_N_plus_L, Q_N_minus_I);
F_N_plus =  F_N(Q_NU_plus_I,Q_NU_minus_R, Q_N_plus_I, Q_N_minus_R);
F_N_minus = F_N(Q_NU_plus_L,Q_NU_minus_I, Q_N_plus_L, Q_N_minus_I);

%Compute the updated Q
Q_NU = Q_NU - grid.dt/(grid.dx)*(F_NU_plus - F_NU_minus);
Q_N = Q_N - grid.dt/(grid.dx)*(F_N_plus - F_N_minus);

%Destruct Q into N, Ux, Uy, Uz
[Ux,Uy,Uz] = destruct_UN(Q_NU,N);
Ux = Ux';
Uy = Uy';
Uz = Uz';
N = Q_N';

% Interpolate back
Uy = interp_center_to_edge(Uy,grid);
Uz = interp_center_to_edge(Uz,grid);
N = interp_center_to_edge(N,grid);

end


%Locally defined functions for MUSCL-Handcock:

%Fluxes N
function [Flux] = F_N(Q_NU_R, Q_NU_L, Q_NR, Q_NL)

% compute c for the lax flux
UR2 = ( Q_NU_R(:,1).^2 + Q_NU_R(:,2).^2 + Q_NU_R(:,3).^2 )./(Q_NR.^2);
UL2 = ( Q_NU_L(:,1).^2 + Q_NU_L(:,2).^2 + Q_NU_L(:,3).^2 )./(Q_NL.^2);
vRx = (Q_NU_R(:,1)./Q_NR)./(sqrt(1 + UR2));
vLx = (Q_NU_L(:,1)./Q_NL)./(sqrt(1 + UL2));
c = max( vRx, vLx ); %Must be a vector!

%Rusanov Flux F = vxQ
F_Q_minus = vLx.*Q_NL;
F_Q_plus = vRx.*Q_NR;
Flux = (1/2) * ( F_Q_minus + F_Q_plus  - c.*( Q_NR - Q_NL ) );

end



%Fluxes NU
function [Flux] = F_NU(QR, QL, NR, NL)

% compute c for the lax flux
UR2 = ( QR(:,1).^2 + QR(:,2).^2 + QR(:,3).^2 )./(NR.^2);
UL2 = ( QL(:,1).^2 + QL(:,2).^2 + QL(:,3).^2 )./(NL.^2);
vRx = (QR(:,1)./NR)./(sqrt(1 + UR2));
vLx = (QL(:,1)./NL)./(sqrt(1 + UL2));
c = max( vRx, vLx ); %Must be a vector!

%Rusanov Flux F = vxQ
F_Q_minus = vLx.*QL;
F_Q_plus = vRx.*QR;
Flux = (1/2) * ( F_Q_minus + F_Q_plus  - c.*( QR - QL ) );

end

%Construct Q
function [Q] = construct_UN(N,U)

% Build Q
Q = N.*U;

end


%Destruct Q
function [Ux,Uy,Uz] = destruct_UN(Q,N)

%destruct Q
Ux = Q(:,1)./N;
Uy = Q(:,2)./N;
Uz = Q(:,3)./N;

end

%Reconstruction
function [W_rec,dW] = reconstruct(Wi,Wm,Wp,x,xi,dx)

%Average Dw
dW = ave( Wi - Wm, Wp - Wi );

%Rec with linear:
W_rec = Wi + ((x-xi)/dx)'.*dW;

end

% Averaging
function [dW] = ave( Wm, Wp )

avg_type = "Supebee"; %"standard";
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
function [AQ] = AQ(Q,N,problem)

% Compute A(Q^n), can be vectors, for dt(UN), (matrix?)
if problem == "fluid"
    AQ = 2*(Q(1)/N)*(1/sqrt(1+(Q(1)^2 + Q(2)^2 + Q(3)^2)/(N^2))) + ...
        - ((Q(1)/N)^3)*(1/sqrt(1+(Q(1)^2 + Q(2)^2 + Q(3)^2)/(N^2)));
    % Compute A(Q^n), can be vectors, for dt(N)
elseif problem == "density"
    % Note Qi = U here:
    AQ =  (Q(:,1)).*(1./sqrt(1+(Q(:,1).^2 + Q(:,2).^2 + Q(:,3).^2)));
else
    fprintf("A(Q) Not Found! \n");
    exit();
end
end


% Function Ap(Q)
function [AQ] = Aprim(U,N,problem)
N = N';

% Compute A(Q^n), This is  the velocity (x) in each cell which NU advects with
if problem == "fluid"
    AQ = (U(:,1)).*(1./sqrt(1+(U(:,1).^2 + U(:,2).^2 + U(:,3).^2)));
    % Compute A(Q^n), can be vectors, for dt(N)
elseif problem == "particles"
    fprintf("Particles incomplete\n");
    exit();
else
    fprintf("A(Q) Not Found! \n");
    exit();
end

end

% Compute the edge values
function [W_plus, W_minus] = edges(W_tilde,grid)
%Linear assumption

%Compute W_plus and W_minus
% (Linear)
% W_plus = W_tilde - dW/2;
% W_minus = W_tilde + dW/2;

%Interpolate 1:
W_tilde = W_tilde';
sz = size(W_tilde);
sz2 = [sz(1),sz(2)+1];
W_tilde_interp = zeros(sz2);
W_plus = zeros(sz);
W_minus = zeros(sz);
for i = 1:sz2(1)
    W_tilde_interp(i,:) = interp_center_to_edge_local(W_tilde(i,:));
    W_plus(i,:) = W_tilde_interp(i,2:end);
    W_minus(i,:) = W_tilde_interp(i,1:end-1);
end
W_plus = W_plus';
W_minus = W_minus';


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