% Regression test: on flat ground the terrain-aware solver must reproduce
% the original solver. The original stored the field after step ix in
% column ix; the new solver stores the starter in column 1, so
%   Phi_new(:,ix+1) == Phi_original(:,ix).
% Run from the repository root:  addpath tests tests/reference; test_flat_regression

freq = 100; hS = 10; T = 10; Tlog = 0.2; shear_exp = 0.15; v_ref = 5; z_ref = 80; theta = 0;
Rg = 286.6896552; gam = 1.41; temp0 = 273.15+T;
c0 = sqrt(gam*Rg*temp0); lambda = c0/freq; k0 = 2*pi*freq/c0;
delta_x = lambda/10; delta_z = lambda/10;
nb_x = floor(400/delta_x)+1; nb_z = floor(100/delta_z)+1;
z = ((0:nb_z-1)*delta_z)';
z0 = 0.00001; d = 0;
Beta = Miki(freq,k0,300,0);
Rp = (1-Beta)/(1+Beta);
U0 = sourceSalomons_order2(z,hS,nb_z,k0,Rp);
[Mx,eps0,eps1] = Mach(v_ref,z_ref,z,z0,d,shear_exp,c0,Tlog,temp0,Rg,gam,theta);
gamma_x = 1./sqrt(1-Mx.^2);
tau_bar1 = (Mx.*gamma_x.^2).*(sqrt(1+eps1) - Mx);

Phi_old = main_original(U0,gamma_x,tau_bar1,k0,delta_x,delta_z,nb_x,nb_z,eps0,eps1,Beta,0.8,10);
Phi_new = main(U0,gamma_x,tau_bar1,k0,delta_x,delta_z,nb_x,nb_z,eps0,eps1,Beta,0.8,10);

err = max(max(abs(Phi_new(:,2:end)-Phi_old(:,1:end-1)))) / max(abs(Phi_old(:)));
printf('flat regression: max relative difference = %.3e\n', err);
assert(err < 1e-12, 'flat-ground regression failed');

% terrain with zero slope everywhere but non-zero offset must not change psi
terrain.x = [0 400]; terrain.h = [25 25]; terrain.smooth = 0;
x = (0:nb_x-1)*delta_x;
[H,Hp,Hpf] = terrain_profile(terrain,x);
Phi_off = main(U0,gamma_x,tau_bar1,k0,delta_x,delta_z,nb_x,nb_z,eps0,eps1,Beta,0.8,10,Hpf);
assert(isequal(Phi_off,Phi_new), 'constant-height terrain changed the field');

% pressure.m with Hp = 0 must equal the original pressure formula
n = 2:nb_x-1;
p_ref = abs(exp(1i*k0*x(n)).*((1-Mx).*Phi_new(:,n) + (1i.*Mx./(2.*k0.*delta_x)).*(Phi_new(:,n+1)-Phi_new(:,n-1))));
p_new = pressure(Phi_new,Mx,k0,delta_x,delta_z,x,Hp);
errp = max(max(abs(p_new(:,n)-p_ref)))/max(p_ref(:));
printf('pressure regression: max relative difference = %.3e\n', errp);
assert(errp < 1e-12, 'pressure regression failed');
disp('test_flat_regression: PASSED');
