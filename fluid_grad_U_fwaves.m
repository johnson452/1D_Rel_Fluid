%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Name: Grant Johnson
%Date: 5/17/2023
%Fluid algroithm, V (divergence U)

%Notes:
%-1D
% F_wave Scheme: Leveque Dust Preprint (2003)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Update the quanitites Ux, Uy, Uz (t -> t + dt)
function [Ux,Uy,Uz,Vx,Vy,Vz,N,grid] = fluid_grad_U_fwaves(Ux,Uy,Uz,N,grid)

% Build Q
U_prim = construct_prim(N, Ux, Uy, Uz);
Q = construct_cons(N, Ux, Uy, Uz);
Nx = grid.Nx;

%Iterate over the domain
% I = linspace(1,Nx-1,Nx-1); %DEFAULT
if grid.iter <= 2
    grid.R = mod( linspace(1,Nx,Nx), Nx-1) + 1; %mod( linspace(1,Nx,Nx), Nx) + 1; %Good
    grid.L = mod( linspace(-1,Nx-2,Nx), Nx-1) + 1; %mod( linspace(-1,Nx-2,Nx), Nx) + 1; %Good
end
R = grid.R;
L = grid.L;

%Compute u_hat at the interfaces (L: i-1/2, R: i+1/2)
u_hat_L_vector = non_linear_solve(N,Ux,Uy,Uz,L,Nx,grid.c);

%Calulcate fluxes:
FG_L = Flux_Godunov(U_prim,u_hat_L_vector,grid);
FG_R =  FG_L(:,R);

%Calculate the corrector Fluxes
[FCorr_L,grid] = Flux_Corr(U_prim,u_hat_L_vector,grid);
FCorr_R = FCorr_L(:,R);

%Compute the updated U
FR = FG_R + FCorr_R;
FL = FG_L + FCorr_L;
Q_new = Q - grid.dt/(grid.dx)*(FR - FL);

% Check if the density fails the positivity test
% [N_new, ~, ~, ~] = destruct_cons(grid,Q_new);
% if min(N_new) < 0
%     X = (N_new<0);
%     fprintf("Second order failed in iteration: %d\n",grid.iter)
%     for i = 1:Nx
%         if X(i) == 1
%             FCorr_R(:,i) = zeros(4,1); %FG_R(X(i));
%             FCorr_L(:,i) = zeros(4,1); %FG_L(X(i));
%             FCorr_R(:,i-1) = zeros(4,1); %FG_R(X_L(i));
%             FCorr_L(:,i+1) = zeros(4,1); %FG_L(X_R(i));
%         end
%     end
%     FR = FG_R + FCorr_R;
%     FL = FG_L + FCorr_L;
% end

%Push with the fluxes
Q = Q - grid.dt/(grid.dx)*(FR - FL);

% Destruct Q into it's components
[N, Ux, Uy, Uz] = destruct_cons(grid,Q);

%Save the output
gamma = sqrt(1+(Ux.*Ux + Uy.*Uy + Uz.*Uz)/(grid.c*grid.c));
Vx = Ux./gamma;
Vy = Uy./gamma;
Vz = Uz./gamma;
end

%Locally defined functions for f-wave:
%Construct Q
function [Q] = construct_prim(N, Ux, Uy, Uz)
Q = [N; Ux; Uy; Uz];
end

%Destruct Q into N, Ux, Uy, Uz
function [N, Ux, Uy, Uz] = destruct_prim(Q)
N = Q(1,:);
Ux = Q(2,:);
Uy = Q(3,:);
Uz = Q(4,:);
end

%Construct Q
function [Q] = construct_cons(N, Ux, Uy, Uz)
Q = [N; (N.*Ux); (N.*Uy); (N.*Uz)];
end

%Destruct Q into N, Ux, Uy, Uz
function [N, Ux, Uy, Uz] = destruct_cons(grid,Q)
N = Q(1,:);
Ux = Q(2,:)./N;
Uy = Q(3,:)./N;
Uz = Q(4,:)./N;

% Positivity Limiter
Ux(N<=0) = 0;
Uy(N<=0) = 0;
Uz(N<=0) = 0;
N(N<=0) = 1.e10;
end


%Fluxes Godunov ( u_hat(i) | Q(i) )
function [Fl] = Flux_Godunov(Q,u_hat,grid)

%Grab the elements of Q
[N, Ux, Uy, Uz] = destruct_prim(Q);

%Grab individual ui, ui-1
gamma = sqrt(1+(Ux.*Ux + Uy.*Uy + Uz.*Uz)/(grid.c*grid.c));
Vx = Ux./gamma;
L = grid.L;

%Build Flux Array
Fl = zeros(min(size(Q)),grid.Nx);

% Loop over the grid, making the fluxes
for i = 1:grid.Nx
    if Ux(L(i)) < 0 && 0 < Ux(i)
        Fl(:,i) = [0;0;0;0];
    else
        if u_hat(1,i) > 0
            % f(Q_{i-1})
            F1 = N(L(i))*Vx(L(i));
            F2 = N(L(i))*Vx(L(i))*Ux(L(i));
            F3 = N(L(i))*Vx(L(i))*Uy(L(i));
            F4 = N(L(i))*Vx(L(i))*Uz(L(i));
        elseif u_hat(1,i) == 0
            % 0.5* ( f(Q_{i-1}) + f(Q_{i}) )
            F1i = N(i)*Vx(i);
            F2i = N(i)*Vx(i)*Ux(i);
            F3i = N(i)*Vx(i)*Uy(i);
            F4i = N(i)*Vx(i)*Uz(i);
            F1m = N(L(i))*Vx(L(i));
            F2m = N(L(i))*Vx(L(i))*Ux(L(i));
            F3m = N(L(i))*Vx(L(i))*Uy(L(i));
            F4m = N(L(i))*Vx(L(i))*Uz(L(i));
            F1 = 0.5*(F1i + F1m);
            F2 = 0.5*(F2i + F2m);
            F3 = 0.5*(F3i + F3m);
            F4 = 0.5*(F4i + F4m);
        elseif u_hat(1,i) < 0
            % f(Q_{i})
            F1 = N(i)*Vx(i);
            F2 = N(i)*Vx(i)*Ux(i);
            F3 = N(i)*Vx(i)*Uy(i);
            F4 = N(i)*Vx(i)*Uz(i);
        end
        Fl(:,i) = [F1;F2;F3;F4];
    end
end
end

%Fluxes Godunov ( Q(i-1) | u_hat(i) | Q(i) )
function [Fl,grid] = Flux_Corr(Q,u_hat,grid)

% mw number of waves

%Grab the elements of Q
[N, Ux, Uy, Uz] = destruct_prim(Q);

%Grab individual ui, ui-1
gamma = sqrt(1+(Ux.*Ux + Uy.*Uy + Uz.*Uz)/(grid.c*grid.c));
Vx = Ux./gamma;
gamma_hat = sqrt(1+(u_hat(1,:).*u_hat(1,:) + u_hat(2,:).*u_hat(2,:) + u_hat(3,:).*u_hat(3,:))/(grid.c*grid.c));
Vx_hat = u_hat(1,:)./gamma_hat;
L = grid.L;
R = grid.R;

%Build Flux Array
Fl = zeros(min(size(Q)),grid.Nx);
Z1 = zeros(min(size(Q)),grid.Nx);
Z2 = zeros(min(size(Q)),grid.Nx);
s1 = zeros(1,grid.Nx);
s2 = zeros(1,grid.Nx);

% Loop over the grid, making the fluxes
for i = 1:grid.Nx
    if Ux(L(i)) < 0 && 0 < Ux(i)
        Z1_1 = - N(L(i))*Vx(L(i));
        Z1_2 = - N(L(i))*Vx(L(i))*Ux(L(i));
        Z1_3 = - N(L(i))*Vx(L(i))*Uy(L(i));
        Z1_4 = - N(L(i))*Vx(L(i))*Uz(L(i));
        s1(i) = Vx(L(i));
        Z1(:,i) = [Z1_1;Z1_2;Z1_3;Z1_4];
        Z2_1 = N(i)*Vx(i);
        Z2_2 = N(i)*Vx(i)*Ux(i);
        Z2_3 = N(i)*Vx(i)*Uy(i);
        Z2_4 = N(i)*Vx(i)*Uz(i);
        Z2(:,i) = [Z2_1;Z2_2;Z2_3;Z2_4];
        s2(i) = Vx(i);
    else
        s1(i) = Vx_hat(i);
        s2(i) = Vx_hat(i);
        if u_hat(1,i) < 0
            Z1_1 = N(i)*Vx(i) - N(L(i))*Vx(L(i));
            Z1_2 = N(i)*Vx(i)*Ux(i) - N(L(i))*Vx(L(i))*Ux(L(i));
            Z1_3 = N(i)*Vx(i)*Uy(i) - N(L(i))*Vx(L(i))*Uy(L(i));
            Z1_4 = N(i)*Vx(i)*Uz(i) - N(L(i))*Vx(L(i))*Uz(L(i));
            Z1(:,i) = [Z1_1;Z1_2;Z1_3;Z1_4];
            Z2_1 = 0;
            Z2_2 = 0;
            Z2_3 = 0;
            Z2_4 = 0;
            Z2(:,i) = [Z2_1;Z2_2;Z2_3;Z2_4];
        elseif u_hat(1,i) >= 0
            Z1_1 = 0;
            Z1_2 = 0;
            Z1_3 = 0;
            Z1_4 = 0;
            Z1(:,i) = [Z1_1;Z1_2;Z1_3;Z1_4];
            Z2_1 = N(i)*Vx(i) - N(L(i))*Vx(L(i));
            Z2_2 = N(i)*Vx(i)*Ux(i) - N(L(i))*Vx(L(i))*Ux(L(i));
            Z2_3 = N(i)*Vx(i)*Uy(i) - N(L(i))*Vx(L(i))*Uy(L(i));
            Z2_4 = N(i)*Vx(i)*Uz(i) - N(L(i))*Vx(L(i))*Uz(L(i));
            Z2(:,i) = [Z2_1;Z2_2;Z2_3;Z2_4];
        end
    end
end

Z1_tilde = Z1*0.0;
Z2_tilde = Z2*0.0;
for i = 1:grid.Nx
    %Minmod Limit:
    for j = 1:4
        phi_Z1 = limiter(Z1(j,L(i)),Z1(j,i),Z1(j,R(i)),s1(i));
        phi_Z2 = limiter(Z2(j,L(i)),Z2(j,i),Z2(j,R(i)),s2(i));
        Z1_tilde(j,i) = phi_Z1*Z1(j,i);
        Z2_tilde(j,i) = phi_Z2*Z2(j,i);
   end
%          phi_Z1 = limiter(Z1(:,L(i)),Z1(:,i),Z1(:,R(i)),s1(i));
%          phi_Z2 = limiter(Z2(:,L(i)),Z2(:,i),Z2(:,R(i)),s2(i));
%          Z1_tilde(:,i) = phi_Z1*Z1(:,i);
%          Z2_tilde(:,i) = phi_Z2*Z2(:,i);
end

% Compute the fluxes
for i = 1:4
    Fl(i,:) = 0.5*sign(s1).*(1-(grid.dt/grid.dx)*abs(s1)).*Z1_tilde(i,:) +...
              0.5*sign(s2).*(1-(grid.dt/grid.dx)*abs(s2)).*Z2_tilde(i,:);
end

end

%Limit the waves
function phi = limiter(wL,wI,wR,s)

%Iterate over the waves:
dotl = dot_manual(wL,wI);
dotr = dot_manual(wI,wR);
wnorm2 = dot_manual(wI,wI);

% Check for zero norm:
if wnorm2 > 0
    if s > 0
        r = dotl/wnorm2;
    else
        r = dotr/wnorm2;
    end
    phi = minmod(r);
else
    phi = 0.0;
end

end


%Compute the manual dotproduct
function W2 = dot_manual(a,b)

sz = max(size(a));
W2 = 0.0;
for i = 1:sz
    W2 = W2 + a(i)*b(i);
end

end

% Minmod used for averaging
function [phi] = minmod(r)
theta = 2.0;
phi = max(0, min([theta*r, (1+r)/2, theta]));
%theta = max(0+zeros(size(r)), min(1+zeros(size(r)), r));
end

% %Fluxes N (MUSCL)
% function [Flux_vec] = Flux(QL,QR,grid)
%
% c = grid.c;
% [NR, UxR, UyR, UzR] = destruct(QR);
% [NL, UxL, UyL, UzL] = destruct(QL);
%
% % compute c for the lax flux
% UR2 = ( UxR.^2 + UyR.^2 + UzR.^2 )/(c^2);
% UL2 = ( UxL.^2 + UyL.^2 + UzL.^2 )/(c^2);
% vRx = UxR./(sqrt(1 + UR2));
% vLx = UxL./(sqrt(1 + UL2));
% c = max( abs(vRx), abs(vLx) ); %%% FIX ABS | lamdba^p|
%
% %Rusanov Flux F = vxQ
% FR = [ NR.*vRx  ; NR.*vRx.*UxR ; NR.*vRx.*UyR  ; NR.*vRx.*UzR ];
% FL = [ NL.*vLx  ; NL.*vLx.*UxL ; NL.*vLx.*UyL  ; NL.*vLx.*UzL ];
% QR_vec = [ NR  ; NR.*UxR ; NR.*UyR  ; NR.*UzR ];
% QL_vec = [ NL  ; NL.*UxL ; NL.*UyL  ; NL.*UzL ];
% Flux_vec = (1/2) * ( FR + FL  - c.*( QR_vec - QL_vec ) );
%
% end