% Usage: D = diff1(n)
%
% Return an (n+1) x n matrix representing center-difference differentiation
% on a grid of n points with spacing 1 and Dirichlet (0) boundary conditions,
% where the n+1 resulting points (the column space of D) correspond to
% the approximate derivatives midway between the grid points, including
% the boundaries (i.e. the points at 0.5, 1.5, ...., n+0.5).
%
% The n x n discrete 1d Laplacian matrix is -diff1(n)' * diff1(n)
function D = diff1(n)
  D = -diag(ones(n,1)) + diag(ones(n-1,1),1);
  D = [ [1,zeros(1,n-1)]; D];
