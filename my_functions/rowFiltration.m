%%%%%%%%%% Filtration %%%%%%%%%
function [On, x_sing] = rowFiltration(Delta, Omega, x, opt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m = size(Delta, 2);
p = size(Omega, 1);
n = length(x);

Omegan = Omega;
while rank(Omegan) < n
    
    % Cicle for compute every Lie Bracket of Distributions
    for i = 1:p
        for j = 1:m
            candidateLieBracket = [Omegan; rowLieBracket(Omegan(i, :), Delta(:, j), x)];
            
            % Verify that this LieBracket is linear indipendent
            if(rank(candidateLieBracket) > rank(Omegan))
                Omegan = candidateLieBracket;
            end
    
        end
    end

    % Check if dim(Omegak) == dim(Omegak-1)
    if(p == size(Omegan, 1))
        break;
    else
        % Update number of vectors that spans Omegan
        p = size(Omegan, 1);
    end
end

if(p < n)
    disp("Dimension of distribution: " + num2str(p));

else
    % Check if there is a singular value
    x_sing = struct2cell(solve(det(Omegan) == 0, x));
    isSing = true;
    
    for i = 1:n
        isSing = isSing && (~isempty(x_sing{i}));
        if(isSing)
            disp("Singular Value!")
            break;
        end
    end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
On = Omegan;
end