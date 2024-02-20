%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Name: Grant Johnson
%Date: 5/17/2023
%Fluid algroithm, V (divergence U)

%Notes:
%-1D
% F_wave Scheme: Leveque Dust Preprint (2003)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Update the quanitites Ux, Uy, Uz (t -> t + dt)
function [Ux,Uy,Uz,Vx,Vy,Vz,N,grid] = fluid_grad_U_fwaves_eigen(Ux,Uy,Uz,N,grid)

% Build Q
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

%Calulcate fluxes (Between cells i-1 and i)
[amdq,apdq] = Flux_FO(Q,grid,0,0);
apdq_l = apdq;
amdq_r = amdq(:,R);


% TEMP DIAGNOSTIC
Nx = grid.Nx;
Vx = zeros(1,Nx);
Fx = zeros(4,Nx);
conserved = zeros(4,Nx);
for i = 1:Nx
    Vx(i) = Ux(i)/sqrt(1.0 + (Ux(i)*Ux(i) + Uy(i)*Uy(i) + Uz(i)*Uz(i))/(grid.c*grid.c));
    Fx(:,i) = Vx(i)*Q(:,i);
end
% for i = 1:Nx
%     conserved(:,i) = (Fx(:,i) - Fx(:,L(i))) - (amdq(:,i) + apdq(:,i));
%     for j = 1:4
%         if max(abs(Fx(:,i))) > 0
%             if max(abs(conserved(:,i)))/max(abs((Fx(:,i) - Fx(:,L(i))))) > 1e-5 || ...
%                  max(abs(conserved(:,i)))/max(abs(amdq(:,i) + apdq(:,i))) > 1e-5
%                 fprintf("\n\n *********** Index: %d ************ \n",i)
%                 fprintf("f: \n");
%                 disp(Fx(:,i))
%                 fprintf("fl: \n");
%                 disp(Fx(:,L(i)))
%                 fprintf("amdq: \n");
%                 disp(amdq(:,i))
%                 fprintf("apdq: \n");
%                 disp(apdq(:,i))
%                 fprintf("f - fl: \n");
%                 disp((Fx(:,i) - Fx(:,L(i))))
%                 fprintf("amdq + apdq: \n");
%                 disp((amdq(:,i) + apdq(:,i)))
%                 fprintf("(f - fl) - (amdq + apdq): \n");
%                 disp(conserved(:,i))
%             end
%         end
%     end
% end
%     clf()
% for j = 1:4
%     subplot(1,4,j)
%     plot((Fx(j,:) - Fx(j,L(:))),'.')
%     hold on
%     plot( (amdq(j,:) + apdq(j,:)),'.' )
%     legend("F","A")
% end
% TEMP DIAGNOSTIC

%Calculate the corrector Fluxes
[FCorr_L,grid,wv,wv_limited] = Flux_Corr(Q,grid);
FCorr_R = FCorr_L(:,R);

%Compute the updated U
FR = 0.0*FCorr_R;
FL = 0.0*FCorr_L;
Q_new = Q - (grid.dt/grid.dx)*((apdq_l + amdq_r) + (FR - FL));

% clf()
% for j = 1:4
%     subplot(2,2,j)
%     plot(((apdq_l(j,:) + amdq_r(j,L))))
%     hold on
%     plot(apdq_l(j,:))
%     hold on
%     plot(amdq_r(j,L))
%     legend("A + A","A+","A-")
% end
% pause(0.1)


% clf()
% for j = 1:4
%     subplot(2,2,j)
%     plot(Fx(j,:) - Fx(j,L))
%     hold on
%     plot(apdq_l(j,:) + amdq_r(j,:))
%     legend("F","A+")
% end
% pause(1.0)

% Check if the density fails the positivity test
[N_new, ~, ~, ~] = destruct_cons(Q_new);
if min(N_new) < 0

    X = (N_new<0);
    fprintf("Second order failed in iteration: %d\n",grid.iter)
    for i = 1:Nx
        if X(i) == 1
            FCorr_R(:,i) = zeros(4,1); %FG_R(X(i));
            FCorr_L(:,i) = zeros(4,1); %FG_L(X(i));
            FCorr_R(:,i-1) = zeros(4,1); %FG_R(X_L(i));
            FCorr_L(:,i+1) = zeros(4,1); %FG_L(X_R(i));

            % Limit the waves:
            wv(:,:,i) = wv_limited(:,:,i);
            wv(:,:,L(i)) = wv_limited(:,:,L(i));
            wv(:,:,R(i)) = wv_limited(:,:,R(i));
        end
    end
    FR = FCorr_R;
    FL = FCorr_L;

    %RE-Calulcate fluxes (Between cells i-1 and i due to positivity violations)
    [amdq,apdq] = Flux_FO(Q,grid,1,wv);
    apdq_l = apdq;
    amdq_r = amdq(:,R);

    Q_new = Q - (grid.dt/grid.dx)*((apdq_l + amdq_r) + (FR - FL));
    [N_new2, ~, ~, ~] = destruct_cons(Q_new);
    if min(N_new2,[],"all") < 0
        fprintf("Wave limiting did not fix the problem!\n");
    end


    fprintf("Total Fluxes fixed: %d\n",sum(X,"all"));
end


if isreal((apdq_l + amdq_r)) == 0
    fprintf("Complex A+ A-\n");
end

%Push with the fluxes
Q = Q - (grid.dt/grid.dx)*((apdq_l + amdq_r) + (FR - FL));

% Destruct Q into it's components
[N, Ux, Uy, Uz] = destruct_cons(Q);

if min(N,[],"all") < 0
    fprintf("N < 0\n");
end

%Save the output
gamma = sqrt(1+(Ux.*Ux + Uy.*Uy + Uz.*Uz)/(grid.c*grid.c));
Vx = Ux./gamma;
Vy = Uy./gamma;
Vz = Uz./gamma;
end

%Construct Q
function [Q] = construct_cons(N, Ux, Uy, Uz)
Q = [N; (N.*Ux); (N.*Uy); (N.*Uz)];
end

%Destruct Q into N, Ux, Uy, Uz
function [N, Ux, Uy, Uz] = destruct_cons(Q)
N = Q(1,:);
Ux = Q(2,:)./N;
Uy = Q(3,:)./N;
Uz = Q(4,:)./N;

% Positivity Limiter
% Ux(N<=0) = 0;
% Uy(N<=0) = 0;
% Uz(N<=0) = 0;
N(N<=0) = 1;
end

%Destruct Q into N, Ux, Uy, Uz
function [N, Ux, Uy, Uz] = destruct_cons_special(Q)
N = Q(1,:);
N(N<=0) = 1;
Ux = Q(2,:)./N;
Uy = Q(3,:)./N;
Uz = Q(4,:)./N;
end


%Fluxes Godunov ( u_hat(i) | Q(i) )
function [amdq,apdq] = Flux_FO(Q,grid,which_wave,wv_limited)

%Compute wv1 and wv2
[wv,s] = compute_wv(Q,grid);
if which_wave == 1
    wv = wv_limited;
end

amdq = zeros(min(size(Q)),grid.Nx);
apdq = zeros(min(size(Q)),grid.Nx);

% Loop over the grid, making the fluxes
for mw = 1:2
    for i = 1:grid.Nx
        if s(mw,i) < 0
            Z1_1 = wv(mw,1,i);
            Z1_2 = wv(mw,2,i);
            Z1_3 = wv(mw,3,i);
            Z1_4 = wv(mw,4,i);
            amdq(:,i) = amdq(:,i) + [Z1_1;Z1_2;Z1_3;Z1_4];
            Z2_1 = 0;
            Z2_2 = 0;
            Z2_3 = 0;
            Z2_4 = 0;
            apdq(:,i) = apdq(:,i) + [Z2_1;Z2_2;Z2_3;Z2_4];
        elseif s(mw,i) > 0
            Z1_1 = 0;
            Z1_2 = 0;
            Z1_3 = 0;
            Z1_4 = 0;
            amdq(:,i) = amdq(:,i) + [Z1_1;Z1_2;Z1_3;Z1_4];
            Z2_1 = wv(mw,1,i);
            Z2_2 = wv(mw,2,i);
            Z2_3 = wv(mw,3,i);
            Z2_4 = wv(mw,4,i);
            apdq(:,i) = apdq(:,i) + [Z2_1;Z2_2;Z2_3;Z2_4];
        else
            %Split half and half
            Z1_1 = 0.5*wv(mw,1,i);
            Z1_2 = 0.5*wv(mw,2,i);
            Z1_3 = 0.5*wv(mw,3,i);
            Z1_4 = 0.5*wv(mw,4,i);
            amdq(:,i) = amdq(:,i) + [Z1_1;Z1_2;Z1_3;Z1_4];
            Z2_1 = 0.5*wv(mw,1,i);
            Z2_2 = 0.5*wv(mw,2,i);
            Z2_3 = 0.5*wv(mw,3,i);
            Z2_4 = 0.5*wv(mw,4,i);
            apdq(:,i) = apdq(:,i) + [Z2_1;Z2_2;Z2_3;Z2_4];
        end
        if i == grid.Nx-1
%             subplot(1,2,1)
%             plot(squeeze(wv(1,1,:)))
%             hold on
%             plot(squeeze(wv(2,1,:)))
%             subplot(1,2,2)
%             plot(apdq(1,:))
%             hold on
%             plot(amdq(1,:))
%             fprintf("Stop\n");
        end
    end
end

end

%Fluxes Corrective ( Q(i-1) | u_hat(i) | Q(i) )
function [Fl,grid,wv,wv_limited] = Flux_Corr(Q,grid)

%Build Flux Array
Fl = zeros(min(size(Q)),grid.Nx);
L = grid.L;
R = grid.R;


%Compute wv1 and wv2
[wv,s] = compute_wv(Q,grid);

s1 = s(1,:);
s2 = s(2,:);
Z1 = squeeze(wv(1,:,:));
Z2 = squeeze(wv(2,:,:));


Z1_tilde = Z1*0.0;
Z2_tilde = Z2*0.0;
for i = 1:grid.Nx
    %Minmod Limit:
    %     for j = 1:4
    %         phi_Z1 = limiter(Z1(j,L(i)),Z1(j,i),Z1(j,R(i)),s1(i));
    %         phi_Z2 = limiter(Z2(j,L(i)),Z2(j,i),Z2(j,R(i)),s2(i));
    %         Z1_tilde(j,i) = phi_Z1*Z1(j,i);
    %         Z2_tilde(j,i) = phi_Z2*Z2(j,i);
    %     end
    phi_Z1 = limiter(Z1(:,L(i)),Z1(:,i),Z1(:,R(i)),s1(i));
    phi_Z2 = limiter(Z2(:,L(i)),Z2(:,i),Z2(:,R(i)),s2(i));
    Z1_tilde(:,i) = phi_Z1*Z1(:,i);
    Z2_tilde(:,i) = phi_Z2*Z2(:,i);
end

wv_limited = wv;
wv_limited(1,:,:) = Z1_tilde;
wv_limited(2,:,:) = Z2_tilde;

% Compute the fluxes
for i = 1:4
    Fl(i,:) = 0.5*sign(s1).*(1-(grid.dt/grid.dx)*abs(s1)).*Z1_tilde(i,:) +...
        0.5*sign(s2).*(1-(grid.dt/grid.dx)*abs(s2)).*Z2_tilde(i,:);
end

end


% Compute Z1 and Z2:
function [wv,s] = compute_wv(Q,grid)

% mw number of waves
c = grid.c;

%Grab R and L
L = grid.L;
Nx = grid.Nx;

%Grab the elements of Q
[rhol, ulx, uly, ulz] = destruct_cons_special(Q(:,L));
[rhor, urx, ury, urz] = destruct_cons_special(Q);

% compute the constants:
gammal = sqrt(1.0 + (ulx.*ulx + uly.*uly + ulz.*ulz)/(c*c));
gammar = sqrt(1.0 + (urx.*urx + ury.*ury + urz.*urz)/(c*c));
vlx = ulx./gammal;
vrx = urx./gammar;
vly = uly./gammal;
vry = ury./gammar;
vlz = ulz./gammal;
vrz = urz./gammar;

%Primative rho
rhol_prim = rhol./(gammal);
rhor_prim = rhor./(gammar);

if min(rhol_prim,[],"all") < 0
    fprintf("N < 0\n");
end

k = (sqrt(rhol_prim) + sqrt(rhor_prim))/(2.0*c);
w0 = (sqrt(rhol_prim).*gammal + sqrt(rhor_prim).*gammar)/(2.0);
w1 = (sqrt(rhol_prim).*gammal.*vlx/c + sqrt(rhor_prim).*gammar.*vrx/c)/(2.0);
w2 = (sqrt(rhol_prim).*gammal.*vly/c + sqrt(rhor_prim).*gammar.*vry/c)/(2.0);
w3 = (sqrt(rhol_prim).*gammal.*vlz/c + sqrt(rhor_prim).*gammar.*vrz/c)/(2.0);

dQ = Q - Q(:,L);
d0 = dQ(1,:);
d1 = dQ(2,:);
d2 = dQ(3,:);
d3 = dQ(4,:);

wv = zeros(2,4,Nx);
s = zeros(2,Nx);
wv(1,1,:) =  (c .* (d1 .* c .* c .* k .* k .* k - d0 .* c .* c .* k .* k .* w1 + d1 .* k .* w0 .* w0 - d1 .* k .* w1 .* w1 - 2 .* d2 .* k .* w1 .* w2 - 2 .* d3 .* k .* w1 .* w3 + d1 .* k .* w2 .* w2 + d1 .* k .* w3 .* w3 - d0 .* w0 .* w0 .* w1 + d0 .* w1 .* w1 .* w1 + d0 .* w1 .* w2 .* w2 + d0 .* w1 .* w3 .* w3)) ./ (w0 .* (c .* c .* k .* k - w0 .* w0 + w1 .* w1 + w2 .* w2 + w3 .* w3));
wv(1,2,:) =  (2 .* c .* w1 .* (d1 .* c .* c .* k .* k - d0 .* w1 .* c .* c .* k + d1 .* w2 .* w2 - d2 .* w1 .* w2 + d1 .* w3 .* w3 - d3 .* w1 .* w3)) ./ (w0 .* (c .* c .* k .* k - w0 .* w0 + w1 .* w1 + w2 .* w2 + w3 .* w3));
wv(1,3,:) =  (c .* (d2 .* c .* c .* k .* k .* w1 + d1 .* c .* c .* k .* k .* w2 - 2 .* d0 .* c .* c .* k .* w1 .* w2 - d2 .* w0 .* w0 .* w1 + d1 .* w0 .* w0 .* w2 + d2 .* w1 .* w1 .* w1 - d1 .* w1 .* w1 .* w2 - d2 .* w1 .* w2 .* w2 - 2 .* d3 .* w1 .* w2 .* w3 + d2 .* w1 .* w3 .* w3 + d1 .* w2 .* w2 .* w2 + d1 .* w2 .* w3 .* w3)) ./ (w0 .* (c .* c .* k .* k - w0 .* w0 + w1 .* w1 + w2 .* w2 + w3 .* w3));
wv(1,4,:) =  (c .* (d3 .* c .* c .* k .* k .* w1 + d1 .* c .* c .* k .* k .* w3 - 2 .* d0 .* c .* c .* k .* w1 .* w3 - d3 .* w0 .* w0 .* w1 + d1 .* w0 .* w0 .* w3 + d3 .* w1 .* w1 .* w1 - d1 .* w1 .* w1 .* w3 + d3 .* w1 .* w2 .* w2 - 2 .* d2 .* w1 .* w2 .* w3 - d3 .* w1 .* w3 .* w3 + d1 .* w2 .* w2 .* w3 + d1 .* w3 .* w3 .* w3)) ./ (w0 .* (c .* c .* k .* k - w0 .* w0 + w1 .* w1 + w2 .* w2 + w3 .* w3));
s(1,:) = (c.*w1)./w0;
wv(2,1,:) =  -(2 .* c .*  k .* w0 .* (d1 .* c .* c .* k .* k - 2 .* d0 .* c .* c .* k .* w1 + d1 .* w0 .* w0 - d1 .* w1 .* w1 - 2 .* d2 .* w1 .* w2 - 2 .* d3 .* w1 .* w3 + d1 .* w2 .* w2 + d1 .* w3 .* w3)) ./ ((c .* c .* k .* k - w0 .* w0 + w1 .* w1 + w2 .* w2 + w3 .* w3) .* (c .* c .* k .* k + w0 .* w0 + w1 .* w1 + w2 .* w2 + w3 .* w3));
wv(2,2,:) =  -(2 .* c .* w0 .* w1 .* (d1 .* c .* c .* k .* k - 2 .* d0 .* c .* c .* k .* w1 + d1 .* w0 .* w0 - d1 .* w1 .* w1 - 2 .* d2 .* w1 .* w2 - 2 .* d3 .* w1 .* w3 + d1 .* w2 .* w2 + d1 .* w3 .* w3)) ./ ((c .* c .* k .* k - w0 .* w0 + w1 .* w1 + w2 .* w2 + w3 .* w3) .* (c .* c .* k .* k + w0 .* w0 + w1 .* w1 + w2 .* w2 + w3 .* w3));
wv(2,3,:) =  -(2 .* c .* w0 .* w2 .* (d1 .* c .* c .* k .* k - 2 .* d0 .* c .* c .* k .* w1 + d1 .* w0 .* w0 - d1 .* w1 .* w1 - 2 .* d2 .* w1 .* w2 - 2 .* d3 .* w1 .* w3 + d1 .* w2 .* w2 + d1 .* w3 .* w3)) ./ ((c .* c .* k .* k - w0 .* w0 + w1 .* w1 + w2 .* w2 + w3 .* w3) .* (c .* c .* k .* k + w0 .* w0 + w1 .* w1 + w2 .* w2 + w3 .* w3));
wv(2,4,:) =  -(2 .* c .* w0 .* w3 .* (d1 .* c .* c .* k .* k - 2 .* d0 .* c .* c .* k .* w1 + d1 .* w0 .* w0 - d1 .* w1 .* w1 - 2 .* d2 .* w1 .* w2 - 2 .* d3 .* w1 .* w3 + d1 .* w2 .* w2 + d1 .* w3 .* w3)) ./ ((c .* c .* k .* k - w0 .* w0 + w1 .* w1 + w2 .* w2 + w3 .* w3) .* (c .* c .* k .* k + w0 .* w0 + w1 .* w1 + w2 .* w2 + w3 .* w3));
s(2,:) =  2.0.*c.*w0.*w1./( c.*c.*k.*k + w0.*w0 + w1.*w1 + w2.*w2 + w3.*w3 );

% Limit unphysical answers
% Check for NaN and replace with 0
isnan_loc_s = isnan(s);
s(isnan_loc_s) = 0;

% Check for Inf and -Inf and replace with 0
isinf_loc_s = isinf(s);
s(isinf_loc_s) = 0;

% Check for NaN and replace with 0
isnan_loc_wv = isnan(wv);
wv(isnan_loc_wv) = 0;

% Check for Inf and -Inf and replace with 0
isinf_loc_wv = isinf(wv);
wv(isinf_loc_wv) = 0;


% clf()
% if max(isnan_loc_s,[],"all") == 1 || max(isnan_loc_wv,[],"all") == 1 || ...
%     max(isinf_loc_s,[],"all") == 1 || max(isinf_loc_wv,[],"all") == 1
%     for i = 1:4
%         subplot(2,3,i)
%         plot(squeeze(isinf_loc_wv(1,i,:)))
%         hold on
%         plot(squeeze(isnan_loc_wv(1,i,:)))
%         hold on
%         plot(squeeze(isinf_loc_wv(2,i,:)))
%         hold on
%         plot(squeeze(isnan_loc_wv(2,i,:)))
%         legend("wv1_inf","wv1_nan","wv2_inf","wv2_nan")
%     end
%     subplot(2,3,5)
%     plot(isinf_loc_s(1,:))
%     hold on
%     plot(isnan_loc_s(1,:))
%     legend("s1_inf","s1_nan")
%     subplot(2,3,6)
%     plot(isinf_loc_s(2,:))
%     hold on
%     plot(isnan_loc_s(2,:))
%     legend("s2_inf","s2_nan")
%     pause(0.1)
% end


% Check for complex answers (shouldn't be allowed)
if min(isreal(wv),[],"all") == 0
    fprintf("wv are not real!\n");
end
if min(isreal(s),[],"all") == 0
    fprintf("wv are not real!\n");
end


%plot around 6503
% clf()
% for j = 1:4
%     subplot(2,3,j)
%     plot(squeeze(wv(1,j,:)))
%     hold on
%     plot(squeeze(wv(2,j,:)))
%     legend("wv1","wv2")
% end
% subplot(2,3,5)
% plot(squeeze(s(1,:)))
% hold on
% plot(squeeze(s(2,:)))
% subplot(2,3,6)
% plot(squeeze(s(1,:) - s(2,:)))
% pause(1.0)

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