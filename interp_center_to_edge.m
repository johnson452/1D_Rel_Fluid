function [x_interp] = interp_center_to_edge(x)

%interpolate the edge cells to centered
Nx = max(size(x));
x_interp = zeros(1,Nx+1);
for i = 2:Nx-1
    x_interp(i) = (x(i) + x(i+1))/2;
end

%Handle the edge case: Periodic, but 
% is overwritten by BC is non-periodic
x_interp(1) = (x(1)+x(Nx))/2;
x_interp(Nx+1) = (x(1)+x(Nx))/2;


end