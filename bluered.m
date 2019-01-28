function bluered(n)
  if nargin < 1
    n = 64;
  end
  r = linspace(0,1,n/2)';
  map = [r, r, ones(size(r))];
  map1 = flipud(fliplr(map(1:end-1,:)));
  map = [map; map1];
  colormap(map);
  c = max(abs(caxis));
  caxis([-c c]);
