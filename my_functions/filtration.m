%%%%%%%%%% Filtration %%%%%%%%%
function [Dn, x_sing] = filtration(Delta, Delta0, x, opt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m = size(Delta, 2);
p = size(Delta0, 2);
n = length(x);

Deltan = Delta0;
while rank(Deltan) < n
    
    % Cicle for compute every Lie Bracket of Distributions
    for i = 1:p
        for j = 1:m
            candidateLieBracket = [Deltan lieBracket(Deltan(:, i), Delta(:, j), x)];
            
            % Verify that this LieBracket is linear indipendent
            if(rank(candidateLieBracket) > rank(Deltan))
                Deltan = candidateLieBracket;
            end
    
        end
    end

    % Check if dim(Deltak) == dim(Deltak-1)
    if(p == size(Deltan, 2))
        break;
    else
        % Update number of vectors that spans Deltan
        p = size(Deltan, 2);
    end
end

if(p < n)
    disp("Dimension of distribution: " + num2str(p));

else
    % Check if there is a singular value
    x_sing = struct2cell(solve(det(Deltan) == 0, x));
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
Dn = Deltan;
end