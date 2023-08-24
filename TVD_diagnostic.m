%%%%%%%%%%%%%%%%%%%%%%%%%%
% TVD Test
%%%%%%%%%%%%%%%%%%%%%%%%%%

function [grid] = TVD_diagnostic(t_str,N,NUx,NUy,NUz,grid)

% for the first iteration
if grid.iter == 2
    grid.old_TV_of_N = 0.0;
    grid.old_TV_of_NUx = 0.0;
    grid.old_TV_of_NUy = 0.0;
    grid.old_TV_of_NUz = 0.0;
end

% Sum over the domain
TV_of_N = sum(abs(diff(N)));
TV_of_NUx = sum(abs(diff(NUx)));
TV_of_NUy = sum(abs(diff(NUy)));
TV_of_NUz = sum(abs(diff(NUz)));

% TVD test:
if t_str ~= "start"
    %Print the TVD
    fprintf("\n---- TVD at the %s of the MUSCL ----\n",t_str);
    str = str_for_TVD(TV_of_N,grid.old_TV_of_N);
    fprintf("TVD (N) n+1: %1.5e <= n: %1.5e (%s)\n",TV_of_N,grid.old_TV_of_N,str)
    str = str_for_TVD(TV_of_NUx,grid.old_TV_of_NUx);
    fprintf("TVD (NUx) n+1: %1.5e <= n: %1.5e (%s)\n",TV_of_NUx,grid.old_TV_of_NUx,str)
    str = str_for_TVD(TV_of_NUy,grid.old_TV_of_NUy);
    fprintf("TVD (NUy) n+1: %1.5e <= n: %1.5e (%s)\n",TV_of_NUy,grid.old_TV_of_NUy,str)
    str = str_for_TVD(TV_of_NUz,grid.old_TV_of_NUz);
    fprintf("TVD (NUz) n+1: %1.5e <= n: %1.5e (%s)\n",TV_of_NUz,grid.old_TV_of_NUz,str)
end

%Update the old TV
grid.old_TV_of_N = TV_of_N;
grid.old_TV_of_NUx = TV_of_NUx;
grid.old_TV_of_NUy = TV_of_NUy;
grid.old_TV_of_NUz = TV_of_NUz;

end

function [str] = str_for_TVD(a,b)
    if a > b
        str = "VIOLATED TVD";
    else
        str = "OK";
    end
end


% function [str] = print_str(a,b,str)
% 
% % Create a formatted string using sprintf
% if a <= b
%     formattedString = sprintf('\x1b[32m TVD (%s) n+1: %1.5e <= n: %1.5e \x1b[0m\n',str,a,b);  % "Hi" in green
% else
%     formattedString = sprintf('\x1b[31m TVD (%s) n+1: %1.5e <= n: %1.5e \x1b[0m\n',str,a,b);  % "Hi" in red
% end
% 
% % Print the formatted string
% fprintf('%s\n', formattedString);
% 
% end