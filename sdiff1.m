% Usage: D = sdiff1(n)
%
% Return an (n+1) x n sparse matrix for center-difference differentiation
% on a grid of n points with spacing 1 and Dirichlet (0) boundary conditions,
% where the n+1 resulting points (the column space of D) correspond to
% the approximate derivatives midway between the grid points, including
% the boundaries (i.e. the points at 0.5, 1.5, ...., n+0.5).
%
% The n x n discrete 1d Laplacian matrix is -diff1(n)' * diff1(n)
function D = sdiff1(n)
  o = ones(n,1);
  D = spdiags([-o,o], [-1,0], n+1,n);
