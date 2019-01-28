% Usage y = cosines(x, L, a)
% 
% Computes y = sum of cosine series(x),
% where the term with amplitude a(n+1) is cos(n pi x / L),
% except that the a(1) term multiples 1/2.
function y = cosines(x, L, a)
  y = a(1)/2;
  for n = 2:length(a)
    y = y + a(n) * cos((n-1)*pi*x/L);
  end
