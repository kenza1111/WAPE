function [H,Hp,Hpf] = terrain_profile(terrain,x)
% Terrain height and slope on the computational range grid.
%
% terrain.x      : abscissae of the profile points (m), increasing
% terrain.h      : ground height at terrain.x (m), piecewise-linear between
%                  points, constant beyond the first/last point
% terrain.smooth : standard deviation of a Gaussian smoothing window (m),
%                  0 for no smoothing (sharp slope kinks are allowed)
% x              : uniform range grid (m), x(1) = 0
%
% H, Hp          : height and slope dH/dx at x
% Hpf            : slope on the half-step grid (0:2*(nb_x-1))*dx/2,
%                  input of main.m

x = x(:).';
nb_x = numel(x);
dx = x(2)-x(1);

% half-step grid: node ix is at index 2*ix-1
xf = x(1) + (0:2*(nb_x-1))*dx/2;
xc = min(max(xf,terrain.x(1)),terrain.x(end));      % constant extrapolation
Hf = interp1(terrain.x(:).',terrain.h(:).',xc,'linear');

if isfield(terrain,'smooth') && terrain.smooth > 0
    s = terrain.smooth;
    hw = ceil(4*s/(dx/2));                            % half width (points)
    g = exp(-0.5*((-hw:hw)*(dx/2)/s).^2);
    g = g/sum(g);
    Hpad = [Hf(1)*ones(1,hw) Hf Hf(end)*ones(1,hw)]; % replicate padding
    Hf = conv(Hpad,g,'valid');
end

% slope on the half-step grid (centred differences, one-sided at the ends)
Hpf = zeros(size(Hf));
Hpf(2:end-1) = (Hf(3:end)-Hf(1:end-2))/dx;
Hpf(1) = (Hf(2)-Hf(1))/(dx/2);
Hpf(end) = (Hf(end)-Hf(end-1))/(dx/2);

H  = Hf(1:2:end);
Hp = Hpf(1:2:end);

end
