% /*
%  * Copyright 2020-2021 UMRAE
%  *
%  * Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the
%  * European Commission - subsequent versions of the EUPL (the "Licence");
%  *
%  * You may not use this work except in compliance with the Licence.
%  * You may obtain a copy of the Licence at:
%  *
%  * https://joinup.ec.europa.eu/software/page/eupl5
%  *
%  * Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis,
%  * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%  * See the Licence for the specific language governing permissions and limitations under the Licence.
%  */

function [Phi] = main(U0,gamma_x,tau_bar1,k0,delta_x,delta_z,nb_x,nb_z,epsilon0,epsilon1,Beta,haut_a,coeff_a,Hpf,slope_tol)
% Range-marching WAPE solver (Crank-Nicholson, Pade(1,1), moving medium)
% over flat or uneven terrain.
%
% Terrain is handled in terrain-following (sheared) coordinates
% zeta = z - H(x), zeta being the vertical height above the local ground.
% For a local slope s = H'(x) = tan(thetaG), the exact transformation of the
% Helmholtz equation shows that the envelope w defined by
%       phi = w * exp(1i*k0*sin(thetaG)*zeta) * exp(1i*int(k0*sqrt(1+s^2) dx))
% obeys the FLAT-ground WAPE written in the frame aligned with the slope,
% i.e. with a marching step delta_x*sqrt(1+s^2) (distance along the slope)
% and a vertical step delta_z/sqrt(1+s^2) (distance perpendicular to the
% slope). Hence, at each step:
%   - the Crank-Nicholson matrices are those of the flat ground, built with
%     these effective steps (cached, rebuilt only when the slope changes by
%     more than slope_tol);
%   - slope changes (terrain curvature) are applied as phase screens
%     exp(-1i*k0*d(sin(thetaG))*zeta) before and after each CN step
%     (Strang splitting). For small slopes this is the GT-PE curvature term
%     eps_eff = eps - 2*zeta*H'' of [Sack & West 1995].
% For a uniform slope this is the rotated PE ([Collins 1990]); for flat
% ground it reduces exactly to the original solver.
%
% Inputs (new / changed with respect to the flat version)
%   U0        : starting field at x = 0 (stored in column 1 of Phi)
%   Beta      : ground admittance, scalar or 1 x nb_x vector (one value per
%               range column, for mixed ground)
%   Hpf       : terrain slope H'(x) on the half-step grid
%               x = (0:2*(nb_x-1))*delta_x/2 (see terrain_profile.m).
%               Optional: flat ground if omitted or empty.
%   slope_tol : slope change that triggers a rebuild of the CN matrices
%               (default 1e-3). Phase screens always use the exact slope.
%
% Output
%   Phi       : nb_z x nb_x, column ix is the envelope at x = (ix-1)*delta_x
%               (column 1 is the starter). See pressure.m for the
%               conversion to acoustic pressure.

if nargin < 14 || isempty(Hpf), Hpf = zeros(1,2*nb_x-1); end
if nargin < 15 || isempty(slope_tol), slope_tol = 1e-3; end
if isscalar(Beta), Beta = Beta*ones(1,nb_x); end
Hpf = Hpf(:).';

Phi = zeros(nb_z,nb_x); % initialization of variable that stores the velocity potential

zeta = (0:nb_z-1)'*delta_z;       % vertical height above local ground (m)
sinG = Hpf./sqrt(1+Hpf.^2);       % sine of the ground slope angle, half-step grid

%% CN matrices for the first step
s_cur = Hpf(2);
[Amat,Bmat] = build_matrices(s_cur,Beta(2),Beta(1));
BetaA = Beta(2); BetaB = Beta(1);

%% Main calculation
Phi(:,1) = U0; % starter at x = 0

for ix = 2:nb_x

    % slope at mid-step and admittances: rebuild matrices only if needed
    s_mid = Hpf(2*ix-2);
    if abs(s_mid - s_cur) > slope_tol || Beta(ix) ~= BetaA || Beta(ix-1) ~= BetaB
        s_cur = s_mid; BetaA = Beta(ix); BetaB = Beta(ix-1);
        [Amat,Bmat] = build_matrices(s_cur,BetaA,BetaB);
    end

    % terrain curvature, first half step (phase screen)
    dsin = sinG(2*ix-2) - sinG(2*ix-3);
    if dsin ~= 0
        U0 = U0.*exp(-1i*k0*dsin*zeta);
    end

    % first matrix resolution
    U0=Bmat*U0;
    
    % second matrix resolution
    U0 = Amat\U0;

    % terrain curvature, second half step (phase screen)
    dsin = sinG(2*ix-1) - sinG(2*ix-2);
    if dsin ~= 0
        U0 = U0.*exp(-1i*k0*dsin*zeta);
    end
    
    % absorbing layer
    U0 = amort(U0,haut_a,coeff_a);

    % save velocity field after 1 iteration
    Phi(:,ix) = U0; % It is "Phi_hat" of the EWAPE given in eq 79 from [Ostashev et al 2020]
       
end

%% nested function: CN matrices of the flat-ground WAPE in the slope frame
    function [Amat,Bmat] = build_matrices(s,BetaA,BetaB)

        dx = delta_x*sqrt(1+s^2);   % step along the slope
        dz = delta_z/sqrt(1+s^2);   % grid spacing perpendicular to the slope

        %% Amat
        a_vect = (2 - 1i*k0*dx .*(2* gamma_x.^2 - tau_bar1)   )./(8*k0.^2* gamma_x.^2*dz.^2 );

        b_vect = 1 +   epsilon1/4 ...
            -  1i*k0.*dx/2.*(gamma_x.^2 .* epsilon1./2  - (1+ epsilon1/4).* tau_bar1 )...
            -  (2 - 1i*k0*dx .*(2* gamma_x.^2 - tau_bar1)   )./(4*k0.^2* gamma_x.^2*dz.^2 ) ;

        c_vect = a_vect;

        % boundary conditions
        b_vect(1) = b_vect(1) + a_vect(1)* (2*1i*dz*k0*BetaA);
        c_vect(1) = 2*c_vect(1);

        Amat = sparse([1:nb_z 2:nb_z 1:nb_z-1],[1:nb_z 1:nb_z-1 2:nb_z],[b_vect; a_vect(2:nb_z); c_vect(1:nb_z-1)],nb_z,nb_z); % same as b_vect.*A_b + a_vect.*A_a + c_vect.*A_c (Octave compatible)

        %% Bmat
        d_vect = (2 + 1i*k0*dx .*(2* gamma_x.^2 - tau_bar1)   )./(8*k0.^2* gamma_x.^2*dz.^2 );

        e_vect = 1 +   epsilon0/4 ...
            +  1i*k0.*dx/2.*(gamma_x.^2 .* epsilon0./2  - (1+ epsilon0/4).* tau_bar1 )...
            - (2 + 1i*k0*dx .*(2* gamma_x.^2 - tau_bar1)   )./(4*k0.^2* gamma_x.^2*dz.^2 ) ;

        f_vect = d_vect;

        % boundary conditions
        e_vect(1) = e_vect(1) + d_vect(1)* (2*1i*dz*k0*BetaB);
        f_vect(1) = 2*f_vect(1);

        Bmat = sparse([1:nb_z 2:nb_z 1:nb_z-1],[1:nb_z 1:nb_z-1 2:nb_z],[e_vect; d_vect(2:nb_z); f_vect(1:nb_z-1)],nb_z,nb_z); % same as e_vect.*B_e + d_vect.*B_d + f_vect.*B_f (Octave compatible)
    end

end
