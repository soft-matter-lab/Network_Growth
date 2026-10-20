function [c, E] = chainSpring(r, n, b, kT, model)
% Chain tension as a force coefficient c (force vector = c * bond vector)
% and energy E per bond.
%
%   'gaussian': f = 3 kT r / (n b^2)
%   'cohen'   : f = (kT/b) * x(3-x^2)/(1-x^2), x = r/(n b)
%               (Cohen Pade approximation of the inverse Langevin function)
%               E = kT n (x^2/2 - ln(1-x^2))
%
% Beyond x = xc the Cohen force is continued linearly (tangent at xc) so that
% over-stretched chains stay finite.

xc = 0.95;
x = r ./ (n * b);

switch lower(model)
    case 'gaussian'
        c = 3 * kT ./ (n * b^2);
        E = 3 * kT * r.^2 ./ (2 * n * b^2);

    otherwise
        xs = min(x, xc);
        g  = xs .* (3 - xs.^2) ./ (1 - xs.^2);
        gp = 1 + 2 * (1 + xs.^2) ./ (1 - xs.^2).^2;
        dx = max(x - xc, 0);

        gx = g + gp .* dx;                          % dimensionless force kT*gx/b
        c = kT * (gx ./ max(x, 1e-12)) ./ (n * b^2);
        c(x < 1e-12) = 3 * kT ./ (n(x < 1e-12) * b^2);

        E = kT * n .* (xs.^2 / 2 - log(1 - xs.^2) + g .* dx + 0.5 * gp .* dx.^2);
end
end
