function f = naiveDFT(c, x, w)
% f = naiveDFT(c, x, w)
% Computes the naive Discrete Fourier Transform (DFT) of a sequence in O(n^2) operations.
%
% Inputs:
%   c - Complex vector of data points to transform.
%   x - Vector of sample points.
%   w - Vector of sample frequencies.
%
% Output:
%   f - Complex vector of the DFT results.
    
    % check all inputs
    if ~isnumeric(c) || ~isnumeric(x) || ~isnumeric(w)
        error('All inputs must be numeric.');
    end
    if ~isvector(c) || ~isvector(x) || ~isvector(w)
        error('All inputs must be vectors.');
    end
    if length(c) ~= length(x) || length(c) ~= length(w)
        error('All inputs must be of the same length.');
    end
    
    % make vectors vertical
    c = c(:);
    x = x(:);
    w = w(:);
    
    % compute discrete Fourier transform
    N = length(c);

    f = zeros(N,1);
    for j = 1:N
        f(j) = sum(c .* exp(-2*pi*1i*x*w(j)));
    end
    % Optional: compact version for reference
    % f = exp(-2*pi*1i*w*x') * c;
end