% Usage: g = green1d(x, L, xp [, dx])
% 
% Return the Green's function g(x,xp) of
% -d^2/dx^2 on [0,L] with Dirichlet boundary
% conditions.  If the optional argument
% dx is supplied, then a finite-width
% approximation for the delta function is
% used (= 1/dx from xp to xp + dx, 0 elsewhere).
%
% x and xp can be vectors, in which case g
% is a length(x) by length(xp) matrix.
function g = green1d(x, L, xp, dx)
  g = zeros(length(x), length(xp));

  if nargin <= 3
    for i = 1:length(xp)
      y = xp(i);
      g(:,i) = (1-y/L) * x .* (x <= y) - y/L * (x - L) .* (x > y);
    end
  else
    for i = 1:length(xp)
      y = xp(i);
      m = y / dx;
      kappa = -y^2 / (2*dx);
      beta = - 1/L * (y + dx/2);
      gamma = (m+1) + beta;
      alpha = 1 + beta;
      g(:,i) = (alpha*x) .* (x <= y) + (-x.^2/(2*dx) + gamma*x + kappa) .* (x > y) .* (x < (y + dx)) + (beta * (x - L)) .* (x >= (y + dx));
    end
  end
  