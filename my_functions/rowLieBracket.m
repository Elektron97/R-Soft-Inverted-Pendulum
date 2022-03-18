%%%%%%%%%%%%%%% Row Lie Bracket %%%%%%%%%%%
function LfOmega = rowLieBracket(omega, f, x, opt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function compute Directional Drivative of covector   %
% omega along a vector field  f.                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(nargin > 3 && opt == true )
    LfOmega = simplify(f' * (jacobian(omega', x))' + omega*jacobian(f, x));
else
    LfOmega = f' * (jacobian(omega', x))' + omega*jacobian(f, x);
end