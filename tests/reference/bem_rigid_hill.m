function p = bem_rigid_hill(k0,src,xh,hh,rcv,npw)
% Independent 2D reference for tests: line source above a rigid plane z = 0
% carrying a rigid hill z = h(x) >= 0 (h -> 0 at both ends of xh).
% Boundary-element method with the half-plane (Neumann) Green's function
% G_N = (i/4)[H0(k|r-r'|) + H0(k|r-r'_image|)], so only the hill surface is
% discretised. Constant elements, collocation at element centres,
% 4-point Gauss quadrature, adaptive composite rule for near pairs.
%   (1/2) p(r) = p_inc(r) - PV int_hill p dG_N/dn' ds'   (r on the hill)
%         p(r) = p_inc(r) -    int_hill p dG_N/dn' ds'   (r in the air)
% with n' the normal pointing into the ground.
%
% k0  : wavenumber;  src : [xs zs];  xh, hh : hill profile (fine sampling)
% rcv : M x 2 receiver positions.  npw : elements per wavelength (default 10)
% p   : M x 1 complex pressure.

lam = 2*pi/k0;
if nargin < 6, npw = 10; end
% resample the hill uniformly in arc length, ~lambda/10 elements
s = [0 cumsum(hypot(diff(xh),diff(hh)))];
N = ceil(s(end)/(lam/npw));
sn = linspace(0,s(end),N+1);
xn = interp1(s,xh,sn); zn = interp1(s,hh,sn);
xa = xn(1:end-1); xb = xn(2:end); za = zn(1:end-1); zb = zn(2:end);
len = hypot(xb-xa,zb-za);
tx = (xb-xa)./len; tz = (zb-za)./len;
nx = tz; nz = -tx;                 % normal into the ground (fluid above)
xc = (xa+xb)/2; zc = (za+zb)/2;

gp = [-0.8611363116 -0.3399810436 0.3399810436 0.8611363116];
gw = [0.3478548451 0.6521451549 0.6521451549 0.3478548451];

Ginc = @(x,z) 1i/4*(besselh(0,1,k0*hypot(x-src(1),z-src(2))) + ...
                    besselh(0,1,k0*hypot(x-src(1),z+src(2))));

    function K = kernel(x0,z0,self)
        % K(i,j) = int_element_j dG_N/dn'(r0_i, r') ds'
        % 4-point Gauss rule, refined (composite rule) for near pairs, in
        % particular image terms at the hill feet where h -> 0.
        x0 = x0(:); z0 = z0(:);
        M = numel(x0); K = zeros(M,N);
        for b0 = 1:500:M                     % row blocks to limit memory
            r = b0:min(M,b0+499);
            for q = 1:4
                xq = xc + gp(q)*len/2.*tx; zq = zc + gp(q)*len/2.*tz;
                for img = [1 -1]
                    dx = xq - x0(r); dz = img*zq - z0(r);      % r'(image) - r0
                    R = hypot(dx,dz);
                    dRdn = (dx.*nx + dz.*(img*nz))./R;
                    K(r,:) = K(r,:) + gw(q)*len/2 .* (1i/4*(-k0)*besselh(1,1,k0*R)).*dRdn;
                end
            end
        end
        for b0 = 1:500:M
            r = b0:min(M,b0+499);
            Rmin = min(hypot(xc - x0(r), zc - z0(r)), hypot(xc - x0(r), -zc - z0(r)));
            [ii,jj] = find(Rmin < 3*len | (self & (r(:) == 1:N)));
            for m = 1:numel(ii)
                i = r(ii(m)); j = jj(m);
                K(i,j) = elemint(x0(i),z0(i),j, self && i == j, Rmin(ii(m),j));
            end
        end
    end

    function v = elemint(x0,z0,j,skip_direct,Rmin)
        % composite Gauss rule over element j; the direct part of the
        % principal value over the element itself is 0 (straight element)
        nsub = min(256, max(8, ceil(8*len(j)/max(Rmin,1e-6))));
        t = ((0:nsub-1)' + (gp+1)/2)/nsub - 0.5;       % positions in [-1/2,1/2]
        t = t(:); w = repmat(gw/2/nsub,nsub,1); w = w(:);
        xq = xc(j) + t*len(j)*tx(j); zq = zc(j) + t*len(j)*tz(j);
        v = 0;
        for img = [1 -1]
            if img == 1 && skip_direct, continue; end
            dx = xq - x0; dz = img*zq - z0; R = hypot(dx,dz);
            dRdn = (dx*nx(j) + dz*(img*nz(j)))./R;
            v = v + len(j)*sum(w.*(1i/4*(-k0)*besselh(1,1,k0*R)).*dRdn);
        end
    end

Kc = kernel(xc,zc,true);
ps = (0.5*eye(N) + Kc) \ Ginc(xc,zc).';
p = Ginc(rcv(:,1),rcv(:,2)) - kernel(rcv(:,1),rcv(:,2),false)*ps;
end
