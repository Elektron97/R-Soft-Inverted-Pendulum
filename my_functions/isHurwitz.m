%%%%%%%%% isHurwitz %%%%%%%
function is_hurw = isHurwitz(A)

%Output: 
% > 1: As. Stable
% > 0: Marginally Stable
% > -1: Unstable

eig_A = eig(A);

if(any(eig_A > 0))
    is_hurw = -1;
    
elseif(any(eig_A == 0))
    is_hurw = 0;

else
    is_hurw = 1;
end

end