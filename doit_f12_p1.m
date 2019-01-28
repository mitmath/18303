dx = logspace(-4,1,50);
x = 1;
u = @(x) sin(x);
upp = (-u(x+2*dx)+16*u(x+dx)-30*u(x)+16*u(x-dx)-u(x-2*dx)) ./ (12*dx.^2);
loglog(dx, abs(upp + sin(x)), 'r.-', dx, dx.^4, 'k--')
xlabel('{\Delta}x')
legend('|error|', '{\Delta}x^4', 'Location', 'NorthWest')
