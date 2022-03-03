%%%%%%%%%% Test Derivata Matrice %%%%%%
% Da aggiustare
function diffA = diffMatrix(A, v)
% Matrix Size
[m, p] = size(A);

% Vector Size
n = length(v);

diffA = sym(zeros(m, p*n));

for i = 1:m
    for j = 1:p
        diffA(i, j:n) = jacobian(A(i, j), v);
    end
end


end