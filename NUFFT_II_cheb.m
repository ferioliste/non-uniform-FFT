function f = NUFFT_II_cheb(c, x, tol)
% f = NUFFT_II_cheb(c, x, tol)
% Computes the Non-uniform Fast Fourier Transform (NUFFT) of type II using a bivariate Chebyshev expansion.
%
% Inputs:
%   c   - Complex vector of data points to transform.
%   x   - Vector of sample points. Has to be a perturbated equispaced grid
%   tol - Positive scalar that controls the accuracy of the approximation.
%
% Output:
%   f   - Complex vector of the DFT results.
    
    % check all inputs
    if ~isnumeric(c) || ~isnumeric(x)
        error('All inputs must be numeric.');
    end
    if ~isvector(c) || ~isvector(x)
        error('All inputs must be vectors.');
    end
    if length(c) ~= length(x)
        error('All inputs must be of the same length.');
    end
    if ~isnumeric(tol) || ~isscalar(tol) || tol <= 0
        error('The tollerance must be a positive scalar.');
    end
    
    % make vectors vertical
    c = c(:);
    x = x(:);
    
    % compute discrete Fourier transform
    N = length(c);
    w = (0:N-1)';
    gamma = max(0.05, max(abs(x - (0:N-1)'/N))*N);
    
    K = max(3, ceil(5*gamma*exp(lambertw(log(140/tol)/(5*gamma)))));

    e = (0:N-1)'/N;
    u = (repmat(exp(-1i*pi*N*(x-e)), 1, K) .* chebT_eval_1_n((N/gamma)*(x-e), K-1)) * a_coefs(K, gamma);
    v = chebT_eval_1_n((2*w)/N - 1, K-1);
    v(:,1) = v(:,1)/2;

    f = zeros(N,1);
    for r = 1:K
        f = f + (v(:,r) .* fft(u(:,r) .* c));
    end
end

function a = a_coefs(K, gamma)
% Generates the matrix of coefficients a_{pr}, the coefficients a_{0r} are halved

    a = zeros(K,K);
    
    for p = 0:K-1
        for r = 0:K-1
            if mod(abs(p-r),2) == 0
                a(p+1,r+1) = 4*(1i^r)*besselj((r+p)/2, -gamma*pi/2)*besselj((r-p)/2, -gamma*pi/2) * (1 - 0.5*(p==0));
            end
        end
    end
end

function T = chebT_eval_1_n(x, n)
% Evaluate Chebyshev polynomials of degree 0,...,n at points in x,
% using the three-term recurrence relation. Function partially taken from
% https://viewer.mathworks.com/?viewer=plain_code&url=https%3A%2F%2Fjp.mathworks.com%2Fmatlabcentral%2Fmlc-downloads%2Fdownloads%2Fe5abb292-4a80-11e4-9553-005056977bd0%2Fd3323a9d-1c0f-f373-1ffc-69bab360c2c3%2Ffiles%2F%40chebfun%2Fnufft2.m&embed=web

    N = length(x);
    T = zeros(N, n+1);
    
    T(:,1) = ones(N,1);
    if n == 0
        return
    end
    T(:,2) = x;
    
    twoX = 2*x;
    for k = 2:n
        T(:,k+1) = twoX.*T(:,k) - T(:,k-1);
    end
end