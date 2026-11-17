function pp2D = pressure(Phi,Mx,k0,delta_x,delta_z,x,Hp)
% 2D acoustic pressure amplitude from the velocity potential, eq 84 of
% [Ostashev et al 2020] (eq 6 of [Kayser 2023]):
%       p = exp(i k0 x) [ (1-Mx) phi + (i Mx/k0) d(phi)/dx ]
%
% In terrain-following coordinates the solver envelope w relates to the
% physical potential (without exp(i k0 x)) by
%   phi = w * exp(i k0 sin(thetaG) zeta) * exp(i int k0 (sqrt(1+H'^2)-1) dx)
% with tan(thetaG) = H'(x) (see main.m), and the x-derivative at constant
% height z is   d/dx|z = d/dx|zeta - H'(x) d/dzeta.
% With Hp = 0 this reduces exactly to the flat-ground formula.
%
% Phi : nb_z x nb_x field from main.m (column ix at x(ix))
% Mx  : nb_z x 1 Mach number profile
% Hp  : 1 x nb_x terrain slope (optional, flat if omitted)

[nb_z,nb_x] = size(Phi);
if nargin < 7 || isempty(Hp), Hp = zeros(1,nb_x); end
Hp = Hp(:).';
zeta = (0:nb_z-1)'*delta_z;

sinG = Hp./sqrt(1+Hp.^2);
acc = zeros(1,nb_x);                                   % extra phase along the slope
acc(2:end) = cumsum(k0*delta_x*((sqrt(1+Hp(1:end-1).^2)+sqrt(1+Hp(2:end).^2))/2 - 1));
phi = Phi.*exp(1i*k0*zeta*sinG).*exp(1i*acc); % physical potential (without exp(i k0 x))

% d/dx at constant zeta (centred, one-sided at the ends)
dphidx = zeros(nb_z,nb_x);
dphidx(:,2:end-1) = (phi(:,3:end)-phi(:,1:end-2))/(2*delta_x);
dphidx(:,1)   = (phi(:,2)-phi(:,1))/delta_x;
dphidx(:,end) = (phi(:,end)-phi(:,end-1))/delta_x;

% d/dzeta (centred, one-sided at the ends), only needed on slopes
if any(Hp ~= 0)
    dphidz = zeros(nb_z,nb_x);
    dphidz(2:end-1,:) = (phi(3:end,:)-phi(1:end-2,:))/(2*delta_z);
    dphidz(1,:)   = (phi(2,:)-phi(1,:))/delta_z;
    dphidz(end,:) = (phi(end,:)-phi(end-1,:))/delta_z;
    dphidx = dphidx - dphidz.*Hp;
end

pp2D = abs(exp(1i*k0*x(:).').*((1-Mx).*phi + (1i.*Mx./k0).*dphidx)); % (Pa)

end
