function C = christoffel(B, q, q_dot)
n = length(q);
C = sym(zeros(n, n));

for i = 1:n
    for j = 1:n
        for k = 1:n
            gamma(i, j, k) = 0.5*(diff(B(i, j), q(k)) + diff(B(i, k), q(j)) - diff(B(j, k), q(i)));
            C(i, j) = C(i, j) + gamma(i, j, k)*q_dot(k);
        end
    end
end


end