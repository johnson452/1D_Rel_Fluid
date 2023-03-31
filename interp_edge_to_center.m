function [x] = interp_edge_to_center(x)

%interpolate the edge cells to centered
Nx = max(size(x));
x_interp = zeros(Nx-1,1);
for i = 1:N-1
    x_interp(i) = (x(i) + x(i+1))/2;
end

end