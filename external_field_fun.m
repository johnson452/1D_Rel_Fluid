function [Ex,Ey,Ez,Bx,By,Bz] = external_field_fun(Ex,Ey,Ez,Bx,By,Bz,grid)

%Specify external fields on FDTD:
if grid.laser_envelope_model == "true" && grid.BC_type == "WFA" && grid.solve_type_field == "FDTD"
    t = grid.time;
    m = grid.m0;
    q = grid.e0;
    c = grid.c;
    z0 = grid.laser1.position;               % This point is on the laser plane
    tau = grid.laser1.profile_duration;        % The duration of the laser (in s)
    a0 = 2.66;
    x = grid.x2;
    a_env = a0*exp( -((x - c*t -z0).^2)/(c*c*tau*tau) );
    rhs = (m.*(x - c.*t -z0).*a_env.^2)./(q.*tau.*tau.*sqrt(1+a_env.^2));
    Ex = Ex + rhs;
end

end