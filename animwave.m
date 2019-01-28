% Usage: [u,v,x] = animwave(f, M, L, ct, cdtdx)
%
% Animate a solution to the 1d wave equation with speed c,
% discretized using a staggered-grid/leap-frog scheme
% on M points from [0,L) with periodic boundary conditions,
% running for a total time c*t = ct.
%
% Uses timestep c*dt = cdtdx * dx.  For stability, cdtdx
% should be < 1.
%
% f is a function f(x).  The initial conditions approximate
% those of the solution u(x,t) = f(x - ct), v(x,t) = -f(x - ct).

function [u,v,x] = animwave(f, M, L, ct, cdtdx)
  dx = L / M;
  x = [0:M-1]'*dx;
  u = feval(f, x);
  umin = min(u) - (max(u)-min(u)) * 0.2;
  umax = max(u) + (max(u)-min(u)) * 0.2;
  v = -feval(f, x + 0.5*dx * (1 + cdtdx));
  nt = round(ct / (cdtdx * dx));
  hold off;
  for i = 1:nt
    plot(x, u);
    % plot(x, circshift(u,-round((i-1)*cdtdx)));
    axis([0,L,umin,umax]);
    legend(sprintf('ct = %g', (i-1)*cdtdx*dx));
    drawnow;
    v = v + cdtdx * [diff(u); u(1)-u(n)];
    u = u + cdtdx * [v(1)-v(n); diff(v)];
  end

