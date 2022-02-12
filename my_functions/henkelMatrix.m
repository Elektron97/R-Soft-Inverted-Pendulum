function H = henkelMatrix(dim)
for i = 1:dim
    for j = 1:dim
        H(i, j) = 1/(i+j-1);
    end
end
end