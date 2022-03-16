%%%%%%%%%%%%%%% Lie Bracket %%%%%%%%%%%
function LfG = lieBracket(f, g, x, opt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function compute Lie Bracket [f, g] of vectors field %
% f(x), g(x).                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(nargin > 3 && opt == true )
    LfG = simplify(jacobian(g, x)*f - jacobian(f, x)*g);
else
    LfG = jacobian(g, x)*f - jacobian(f, x)*g;
end