function [Ux,Uy,Uz,Vx,Vy,Vz,N,grid] = fluid_grad_U_FCT(Ux,Uy,Uz,N,grid)

%Grab Q, F vectors
Q = construct(N,Ux,Uy,Uz);
%F = [ N.*u  ; p + N.*u.*u ; u.*(N.*E + p) ];
c = (grid.dt/grid.dx);
Nx = grid.Nx;
grid.R = mod( linspace(1,Nx,Nx), Nx-1) + 1; %mod( linspace(1,Nx,Nx), Nx) + 1; %Good
grid.L = mod( linspace(-1,Nx-2,Nx), Nx-1) + 1; %mod( linspace(-1,Nx-2,Nx), Nx) + 1; %Good

%%%%
% FCT:
% (1) Compute F with with a lower-order monotonic scheme
[FL_i_minus_half, FL_i_plus_half] = FL(Ux,Uy,Uz,N,grid);

% (2) Compute F with a higher-order scheme
[FH_i_minus_half, FH_i_plus_half] = FH(Ux,Uy,Uz,N,grid);

% (3) Compute "anti-diffusive flux":
A_i_plus_half = FH_i_plus_half - FL_i_plus_half;
A_i_minus_half = FH_i_minus_half - FL_i_minus_half;

% (4) Compute updated lower order solution (td= transported and diffused)
Q_td = Q - c*(FL_i_plus_half - FL_i_minus_half);

% (5) Limit A_i_plus_half such that the new step is free of extrema
[C_i_minus_half,C_i_plus_half] = C(A_i_plus_half,A_i_minus_half,Q,Q_td,grid);
A_C_i_plus_half = C_i_plus_half.*A_i_plus_half;
A_C_i_minus_half = C_i_minus_half.*A_i_minus_half;

% (6) Q-update the solution with the corrected fluxes
Q = Q_td - c*(A_C_i_plus_half - A_C_i_minus_half);

%Decompose into original prim-variables
[N, Ux, Uy, Uz] = destruct(Q);

%Save the output
gamma = sqrt(1+(Ux.*Ux + Uy.*Uy + Uz.*Uz)/(grid.c*grid.c));
Vx = Ux./gamma;
Vy = Uy./gamma;
Vz = Uz./gamma;

end


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


% FCT Functions
function [c_val_minus,c_val_plus] = C(A_i_plus_half,A_i_minus_half,Q,Q_td,grid)

%Equation 14:
% for i = 1:grid.Nx
%     for j = 1:min(size(Q))
%         cond = ( (A_i_plus_half(j,i)*(Q_td(j,i) - Q_td(j,i)) < 0) && ...
%             ( (A_i_plus_half(j,i)*(Q_td(j,grid.R(grid.R(i))) - Q_td(j,grid.R(i))) < 0) ...
%             || (A_i_plus_half(j,i)*(Q_td(j,i) - Q_td(j,grid.L(i))) < 0) ) );
%         if cond
%             A_i_plus_half(j,i) = 0;
%         end
%     end
% end
% Q_td2 = Q_td(:,grid.L);
% for i = 1:grid.Nx
%     for j = 1:min(size(Q))
%         cond = ( (A_i_minus_half(j,i)*(Q_td2(j,i) - Q_td2(j,i)) < 0) && ...
%             ( (A_i_minus_half(j,i)*(Q_td2(j,grid.R(grid.R(i))) - Q_td2(j,grid.R(i))) < 0) ...
%             || (A_i_minus_half(j,i)*(Q_td2(j,i) - Q_td2(j,grid.L(i))) < 0) ) );
%         if cond
%             A_i_minus_half(j,i) = 0;
%         end
%     end
% end



%Constants:
c = (grid.dx/grid.dt);

%Compute Pi+, Qi+, and Ri+
Zero = zeros(size(Q));
P_plus = max(Zero,A_i_minus_half) - min(Zero,A_i_plus_half);
Q_plus = (max_w(Q,Q_td,grid) - Q_td)*c;
R_plus = zeros(size(Q));
for i = 1:max(size(Q))
    for j = 1:min(size(Q))
        if P_plus(j,i) > 0
            R_plus(j,i) = min(1,Q_plus(j,i)/P_plus(j,i));
        elseif P_plus(j,i) == 0
            R_plus(j,i) = 0;
        else
            fprintf("Error");
        end
    end
end

%Compute Pi-, Qi-, and Ri-
P_minus = max(Zero,A_i_plus_half) - min(Zero,A_i_minus_half);
Q_minus = (Q_td - min_w(Q,Q_td,grid))*c;
R_minus = zeros(size(Q));
for i = 1:max(size(Q))
    for j = 1:min(size(Q))
        if P_minus(j,i) > 0
            R_minus(j,i) = min(1.0,Q_minus(j,i)/P_minus(j,i));
        elseif P_minus(j,i) == 0
            R_minus(j,i) = 0;

        else
            fprintf("Error");
        end
    end
end

if (min(min(Q_minus)) < 0 ) || (min(min(Q_plus)) < 0 )
    fprintf("goofed\n");

end


%Compute C+ value
c_val_plus = zeros(size(Q));
for i = 1:grid.Nx
    for j = 1:min(size(Q))
        if A_i_plus_half(j,i) >= 0
            c_val_plus(j,i) = min( R_plus(j,grid.R(i)), R_minus(j,i) );
        elseif A_i_plus_half(j,i) < 0
            c_val_plus(j,i) = min( R_plus(j,i),R_minus(j,grid.R(i)) );
        end
    end
end

%Compute C- value
c_val_minus = c_val_plus(:,grid.L);

%Max and min of c: okay
% fprintf("Max C: %f, min C: %f\n",max(max(c_val_plus)),min(min(c_val_plus)));

end

% max w
function [w_max] = max_w(w,w_td,grid)

%Compute improved fluxes
w_td = max(w,w_td);

%Find the max
w_max = max( w_td(:,grid.L) ,max(w_td,w_td(:,grid.R)) );

end

% min w
function [w_min] = min_w(w,w_td,grid)

%Compute improved fluxes
w_td = min(w,w_td);

%Find the max
w_min = min( w_td(:,grid.L) ,min(w_td,w_td(:,grid.R)) );

end