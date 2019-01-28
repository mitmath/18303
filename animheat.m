% Usage: [u,x] = animheat(u0, L, t, dtdx2, scheme)
%
% Animate the solution to the heat equation u'' = du/dt with
% Dirichlet boundaries u(0) = u(L) = 0, from u0 at time 0 to
% time t, returning the final solution.
%
% Uses a center-difference discretization of u'' in space, with
% a number of points given by length(u0).
%
% The time-step size is dt = dx^2 * dtdx2, where dx = L / (length(u0)+1).
%
% The time-stepping scheme is determined by the "scheme" argument:
%    scheme = 1: forward-difference (explicit, requires dtdx2 < 0.5)
%    scheme = 2: backward-difference (implicit, unconditionally stable)
%    scheme = 3: Crank-Nicolson (implicit, unconditionally stable)
%
function [u,x] = animheat(u0, L, t, dtdx2, scheme)
  if (nargin < 5)
    scheme = 1
  end
  n = length(u0);
  dx = L / (n+1);
  D = diff1(n);
  A = sparse(-D' * D / dx^2);
  x = [1:n]' * dx;
  u = u0;
  dt = dx^2 * dtdx2;
  nt = round(t / dt);
  hold off;
  for i = 0:nt-1
    plot(x, u);
    legend(sprintf('t = %g', i * dt));
    drawnow;
    if scheme == 1        % forward difference
      u = u + A*(dt*u);
    elseif scheme == 2    % backward difference
      u = (speye(n,n) - A*dt) \ u;
    elseif scheme == 3    % Crank-Nicolson
      u = (speye(n,n) - A*(dt/2)) \ (u + A*(dt/2)*u);
    end
  end
