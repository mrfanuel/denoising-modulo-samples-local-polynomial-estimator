

%% Beware: not vectorized and only for 1d

function g_est = localPoly_denoise(z,x,h)
    %localPoly_denoise - Description
    %
    % Syntax: g_est = localPoly_denoise(x,l)
    % input: z complex n x 1 vector (signal value on the grid nodes)
    %        x real n x 1 vector (grid nodes positions) 
    %        h positive real number (window width)
    % output: g_est complex n x 1 vector (local polynomial estimate on the grid)
    
    n = size(z,1);
        g_est = zeros(n,1);    
    
    
        for j=1:n
            %% calculate B_nx_j matrix
            B_nxj = zeros(l,l);
    
            for i=1:n
                xi = x(i);
                xj = x(j);
                sij = (xj-xi)/h;
                U_ixj = polynomialFeatureMap(sij,l);
                B_nxj = B_nxj + U_ixj * U_ixj' * kernel(si);
            end
            B_nxj = B_nxj/(n*h);
    
            %% calculate a_nx vector
            a_nxj = zeros(l,1);
    
            for i=1:n
                xi = x(i);
                xj = x(j);
                sij = (xj-xi)/h;
                zi = z(i);
                U_ixj = polynomialFeatureMap(sij,l);
                a_nxj = a_nxj + zi * U_ixj * kernel(si);
            end
            a_nxj = a_nxj/(n*h);
    
            %% calculate estimator
    
            theta_hat_nx = lsqminnorm(B_nx,a_nx);
    
            g_est(j) = theta_hat_nx(1);
            % Beware: rounding on the product manifold is done elwehere
        end
    
    end
    
    
    
    function output = polynomialFeatureMap(x,l)
    %polynomialFeatureMap - Description
    %
    % Syntax: output = polynomialFeatureMap(x,l)
    %
    % input:
    % x is a real
    % l is an integer
    % output: l x 1 vector, each entry being a monomial
    
        output = zeros(l,1);
        for i=1:l
            output(i) = (x.^i)/factorial(i);
        end
        
    end
    
    
    function output = kernel(x)
    %kernel - Description
    %
    % Syntax: output = kernel(x)
    % 
    % input: x is a real 
    %
    % output 1 if |x|<1 and zero otherwise
    % i.e., rectangular kernel
    
        output = 0.5*(heaviside(x-1)-heaviside(x+1));
        
    end