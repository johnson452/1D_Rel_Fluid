function [Ex,Ey,Ez,Bx,By,Bz] = muscl_field_push(Ex,Ey,Ez,Bx,By,Bz,grid)

% Apply the muscl scheme to update the fields
% Build Q
Q = construct(Ex,Ey,Ez,Bx,By,Bz);
Nx = grid.Nx;

%Iterate over the domain (assumes periodic, needs BC call after this
%function)
R = mod( linspace(1,Nx,Nx), Nx-1) + 1;
L = mod( linspace(-1,Nx-2,Nx), Nx-1) + 1;

%Reconstruct Q (soln) within one cell via linear model
dQ = reconstruct(Q,Q(:,L),Q(:,R),grid);

%Update solution, primative variables
%Q_tilde = Q - grid.dt/(2*grid.dx)*AQ(Q,grid).*dQ;
A = AQ(grid.c);
Q_tilde = zeros(6,grid.Nx);
for i = 1:Nx
    Q_tilde(:,i) = Q(:,i) - grid.dt/(2*grid.dx)*A*dQ(:,i);
end

%No need to enforce a positivity correction
[Q_plus_I, Q_minus_I] = edges_linear(Q_tilde,dQ);
[Q_plus_R, ~] = edges_linear(Q_tilde(:,R),dQ(:,R));
[~, Q_minus_L] = edges_linear(Q_tilde(:,L),dQ(:,L));

%Update the fluxes
F_R =  Flux(Q_minus_I,Q_plus_R,grid);
F_L = Flux(Q_minus_L,Q_plus_I,grid);

%Compute the updated Q
Q = Q - grid.dt/(grid.dx)*(F_R - F_L);

% Destruct Q into it's components
[Ex,Ey,Ez,Bx,By,Bz] = destruct(Q);

end

%Locally defined functions for MUSCL-Handcock:
%Construct Q
function [Q] = construct(Ex,Ey,Ez,Bx,By,Bz)
Q = [Ex;Ey;Ez;Bx;By;Bz];

end

%Destruct Q into N, Ux, Uy, Uz
function [Ex,Ey,Ez,Bx,By,Bz] = destruct(Q)
Ex = Q(1,:);
Ey = Q(2,:);
Ez = Q(3,:);
Bx = Q(4,:);
By = Q(5,:);
Bz = Q(6,:);
end

%Fluxes for  fields
function [Fl] = Flux(QL,QR,grid)

clight = grid.c;
[~,EyR,EzR,~,ByR,BzR] = destruct(QR);
[~,EyL,EzL,~,ByL,BzL] = destruct(QL);

% compute c for the lax flux
c = clight; % Always the maximum

Zero = zeros(size(BzL));

%Rusanov Flux
FL = [Zero;BzL*clight*clight;-ByL*clight*clight;Zero;-EzL;EyL];
FR = [Zero;BzR*clight*clight;-ByR*clight*clight;Zero;-EzR;EyR];
Fl = (1/2) * ( FR + FL  - c.*( QR - QL ) );

end

%Reconstruction
function [dW] = reconstruct(Wi,Wm,Wp,grid)

%Option:
option = "eigenvector";

%Average Dw
if option == "eigenvector"

    % Compute eigenvectors time differences
    clight = grid.c;
    L = left_eigenvector(clight);
    R = right_eigenvector(clight);

    %Iterate through the grid
    Nx = grid.Nx;
    dW = zeros(6,Nx);
    for i = 1:Nx

        %Compute the differences with the primative vars
        % (Matrix Multiplications)
        delta_i_minus_1 = L*(Wi(:,i) - Wm(:,i));
        delta_i = L*(Wp(:,i) - Wi(:,i));
        dW(:,i) = R*ave( delta_i_minus_1, delta_i );
    end
else %Standard Option
    dW = ave( Wi - Wm, Wp - Wi );
end

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

% Left eigenvectors (Matrix form)
function [L] = left_eigenvector(c)
L = [[0,	1/2,	0,	0,	0,	-c/2];...
    [0,	0,	1/2,	0,	c/2,	0];...
    [0,	1/2,	0,	0,	0,	c/2];...
    [0,	0,	1/2,	0,	-c/2,	0];...
    [1,	0,	0,	0,	0,	0];...
    [0,	0,	0,	1,	0,	0]];
end

% Right eigenvectors (Matrix form)
function [R] = right_eigenvector(c)
R = [[0,	0,	0,	0,	1,	0];...
    [1,	0,	1,	0,	0,	0];...
    [0,	1,	0,	1,	0,	0];...
    [0,	0,	0,	0,	0,	1];...
    [0,	1/c,	0,	-1/c,	0,	0];...
    [-1/c,	0,	1/c,	0,	0,	0]];
end

function [W_plus, W_minus] = edges_linear(W_tilde,dW)
W_plus = W_tilde - dW/2;
W_minus = W_tilde + dW/2;
end

% Function A(Q)
function [A] = AQ(c)

A = [[0,	0,	0,	0,	0,	0];...
    [0,	0,	0,	0,	0,	c^2];...
    [0,	0,	0,	0,	-c^2,	0];...
    [0,	0,	0,	0,	0,	0];...
    [0,	0,	-1,	0,	0,	0];...
    [0,	1,	0,	0,	0,	0]];
end
