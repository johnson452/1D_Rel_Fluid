function [N] = density_func(z)

%Specifty the density
N =  1.e10 + 20.e23*((z*5.e4 + -0.5).*(z>10.e-6).*(z<30.e-6)) + 20.e23*((z>30.e-6));
%N = 20.e23+z*0;

end