% [w,U,x] = waveguide(k, n, M, L, h, a1, b1)
%
% Usage: compute the n smallest-frequency (w) waveguide modes
% of b div(a grad u) = -w^2 u, where u = exp(iky) u_k(x)
% and b,a = b1,a1 in a region of width h around x=0,
% and = 1 otherwise.
%
% Use a computational cell of M points and total width L in the
% x direction, with Dirichlet boundary conditions.  Note that
% the boundary conditions are of negligible importance for the
% guided modes if L is big enough, since the guided modes are
% exponentially localized.  The guided modes are the ones for
% which w < k.
%
% Return a row vector of n eigenvalues w, and a matrix U whose
% n columns are the eigenvectors, for the column vector x of
% coordinates.
%
% k can be a vector, in which case w is a length(k) by n matrix
% of eigenvalues at each k, and U are the eigenvectors for k(end).

function [w,U,x] = waveguide(k, n, M, L, h, a1, b1)
  dx = L / (M+1);
  x = [1:M]'*dx - L/2;
  a = 1 + (a1 - 1) * (abs(x) < h/2);
  b = 1 + (b1 - 1) * (abs(x) < h/2);
  x2 = [0:M]'*dx + dx/2 - L/2;
  a2 = 1 + (a1 - 1) * (abs(x2) < h/2);
  D = sdiff1(M);
  A = spdiags(-b,0,M,M) * D' * spdiags(a2,0,M+1,M+1) * D;
  U = zeros(length(x), n);
  w = zeros(length(k), n);
  for i = 1:length(k)
    Ak = A - k(i)^2 * spdiags(a.*b, 0, M, M);
    [U,S] = eigs(Ak, n, 'sm');
    [lam,ilam] = sort(-diag(S)); % make sure frequencies are in order
    w(i,:) = sqrt(lam);
    U = U(:,ilam);
    % eigs normalizes vectors to length 1, but picks a random sign
    % ... it is convenient to pick a deterministic sign based
    % on average U in first half (1:M/2).
    for j = 1:n
      U(:,j) = U(:,j) * sign(mean(U(1:round(M/2),j)));
    end
  end
