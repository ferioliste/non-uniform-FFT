function f = NUFFT_II_tayl(c, x, tol)
% f = NUFFT_II_tayl(c, x, tol)
% Computes the Non-uniform Fast Fourier Transform (NUFFT) of type II using a Taylor expansion.
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
    gamma = max(abs(x - (0:N-1)'/N))*N;

    K = 20; % hard coded value, an adaptive formula could exist
    
    u = ones(N,K);
    v = ones(N,K);
    e = (0:N-1)'/N;
    for r = 1:K-1
        u(:,r+1) = (N*(x-e)).^r;
        v(:,r+1) = ((-1i)^r)*((2*pi*w/N).^r)/factorial(r);
    end

    f = zeros(N,1);
    for r = 1:K
        f = f + (v(:,r) .* fft(u(:,r) .* c));
    end
end