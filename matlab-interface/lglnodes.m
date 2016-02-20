function [x, w, P] = lglnodes( N )
% [x, w, P] = lglnodes( N )
%
% Computes the Legendre-Gauss-Lobatto nodes, weights and the LGL Vandermonde
% matrix. The LGL nodes are the zeros of (1-x^2)*P'_N(x). Useful for numerical
% integration and spectral methods.
%
% Reference on LGL nodes and weights:
%   C. Canuto, M. Y. Hussaini, A. Quarteroni, T. A. Tang, "Spectral Methods
%   in Fluid Dynamics," Section 2.3. Springer-Verlag 1987
%
% Written by Greg von Winckel - 04/17/2004
% Contact: gregvw@chtm.unm.edu
%
% Takes 1 argument:
%
%   N - Order of the Legendre-Gauss-Lobatto polynomial.
%
% Returns 3 values:
%
%   x - Vector, of length N + 1, of the computed Legendre-Gauss-Lobatto nodes.
%   w - Vector, of length N + 1, of the computed Legendre-Gauss-Lobatto weights.
%   P - Matrix, of size N + 1, containing the Legendre Vandermonde
%       coefficients.

% Truncation plus one.
N1 = N + 1;

% Use the Chebyshev-Gauss-Lobatto nodes as the first guess.
x = cos( pi * (0:N) / N )';

% The Legendre Vandermonde Matrix.
P = zeros( N1, N1 );

% Compute P_(N) using the recursion relation.  Compute its first and second
% derivatives and update x using the Newton-Raphson method.
xold = 2;

while max( abs( x - xold ) ) > eps

    xold = x;

    P(:, 1) = 1;
    P(:, 2) = x;

    for k = 2:N
        P(:, k+1) = ((2*k-1) * x .* P(:, k) - (k-1) * P(:, k-1)) / k;
    end

    x = xold - (x .* P(:, N1) - P(:, N)) ./ (N1 * P(:, N1));
end

w = 2 ./ (N * N1 * P(:, N1).^2);

return
