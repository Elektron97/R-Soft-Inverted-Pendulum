%%%%%%%%%% Test Derivata Matrice %%%%%%
% Da aggiustare
function diffA = diffMatrix(A, v)
% Matrix Size
[m, p] = size(A);
% Vector Size
n = length(v);

diffA = [];
row_diffA = [];

for i = 1:m
    for j = 1:p
        row_diffA = [row_diffA, jacobian(A(i, j), v)];
    end
    diffA = [diffA; row_diffA];
    %Reset row_diffA
    row_diffA = [];
end

% Check dimension of Diff Matrix
% size(diffA)

end