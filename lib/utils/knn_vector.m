function nearest = knn_vector(inp_grid,k)
    n = size(inp_grid,1);
    
    % absolute distance between all test and training data
    dist = abs(repmat(inp_grid,1,n) - repmat(inp_grid(:,1)',n,1));

    % indicies of nearest neighbors
    [~,nearest] = sort(dist,2);
    
    % k nearest
    nearest = nearest(:,1:k);
end