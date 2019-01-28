% Usage: [V,S] = mmdemo(w0,p)
%
% Min-max theorem demo: compute the smallest p eigenvalues of -1/w \nabla^2
% in a 2d "L" shape with Dirichlet boundary conditions, where w(x) is a
% weight function which is 1 everywhere except that it is w0 in a circle
% within the upper-left branch of the 'L'.
%
% Returns the eigenvectors (columns of V) and eigenvalues (diagonals of S).
% Also plots them in a side-by-side surface plot.
%
% w0 can be an array of weights, in which case an animation is performed.

function [V,S,A,G,w2] = mmdemo(w0,p)

  if (nargin < 2)
    p = 2;
  end

  N=20;
  G = numgrid('L', N+2);
  dx = 1 / (N+2);

  x2 = ones(N+2,1) * linspace(0,1, N+2);
  y2 = x2';

  for iw = 1:length(w0)
    w2 = ((x2-0.75).^2 + (y2-0.75).^2 < 0.2.^2) * (w0(iw) - 1) + 1;
    w = w2(G>0);
    
    n = sum(G(:)>0);
    A = spdiags(1./w, 0, n,n) * delsq(G) / dx^2;
    
    [V,S] = eigs(A,p,'SM');
    [lam,i] = sort(diag(S));
    
    clf
    for j = 1:p
      subplot(1,p,j)
      U=G; U(G>0) = V(G(G>0),i(j));
      U = U * sign(U(round(3.1*N/4),round(3.2*N/4))); % pick a consistent sign
      pcolor(U);
      colormap(cool);
      shading interp;
      axis square;
      drawnow;
    end
  end

