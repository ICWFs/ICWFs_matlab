
% Branched Error-function:

function[value] = za_error_function_over_r(x,rc)

value = erf(x/rc)./x;

%indices = find(x < rc);
indices = find(x == 0);
value(indices) = 2/(sqrt(pi) * rc) - 2*x(indices).^2/(3*(sqrt(pi)*rc^3)) + 2*x(indices).^4/(10*sqrt(pi)*rc^5) - 2*x(indices).^6/(42*sqrt(pi)*rc^7);



