%%%%%%%%% isPositiveDefinitiveness %%%%%%%%%%%%
function is_pos = isPositiveDef(m)
eig_m = eig(m);

if(all(eig_m > 0))
    is_pos = true;
else
    is_pos = false;
end

end