%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison of graphs via normalizec count statistics/moments of the adjacency
% matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = moments(A, k, method)
    % A = adjacency matrix of the desired graph
    % k = no. of moments requied
    % method = 'exact' or 'approx': 'exact' returns exact counts and
    % 'approx' returns normalized traces of powers of adjacency matrix
    y = zeros(k,1);
    
    if (strcmp(method,'approx'))
        n = size(A, 1);
        temp = A;
        for t = 1:k
            y(t) = trace(temp)/exp((1+t/2)*log(n)); % this normalization is for dense graphs
            temp = temp*temp;
        end
        
    elseif (strcmp(method, 'exact'))
        n = size(A, 1);
        % need to call those count functions Purna di downloaded
    end 
end