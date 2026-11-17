% Validation of the terrain-following WAPE over hills, against an
% independent 2D boundary-element solution of the Helmholtz equation
% (tests/reference/bem_rigid_hill.m; converged reference values precomputed
% by tests/reference/make_bem_reference.m and stored in bem_reference.mat).
% Homogeneous motionless air, rigid ground, receivers z_r above the local
% ground. Compared quantity: insertion loss of the hill
%   IL = 20 log10 |p_hill(rcv)| - 20 log10 |p_flat(rcv)|
% computed with each method, which removes the source normalisation.
% Run from the repository root:
%   addpath tests tests/reference; test_terrain_vs_bem

ref = load('tests/reference/bem_reference.mat');
IL_bem_all = ref.IL(:,:,:,end);                         % finest BEM
IL_bem_err = ref.IL(:,:,:,end) - ref.IL(:,:,:,end-1);   % BEM self-convergence

freq = ref.freq; hS = ref.hS; dim_x = 1000; dim_z = 250;
T = 10; Rg = 286.6896552; gam = 1.41; temp0 = 273.15+T;
c0 = sqrt(gam*Rg*temp0); lambda = c0/freq; k0 = 2*pi*freq/c0;
delta_x = lambda/10; delta_z = lambda/10;
nb_x = floor(dim_x/delta_x)+1; nb_z = floor(dim_z/delta_z)+1;
x = (0:nb_x-1)*delta_x; zeta = ((0:nb_z-1)*delta_z)';
[Mx,eps0,eps1] = Mach(0,80,zeta,1e-5,0,0.15,c0,0,temp0,Rg,gam,0);
gamma_x = 1./sqrt(1-Mx.^2); tau_bar1 = (Mx.*gamma_x.^2).*(sqrt(1+eps1) - Mx);
U0 = sourceSalomons_order2(zeta,hS,nb_z,k0,1);

% flat-ground PE (for the insertion loss)
Phi_flat = main(U0,gamma_x,tau_bar1,k0,delta_x,delta_z,nb_x,nb_z,eps0,eps1,0,0.8,10);
p_flat = pressure(Phi_flat,Mx,k0,delta_x,delta_z,x);

xh = 0:1:dim_x; xr = ref.xr;
res = zeros(size(ref.hills,1),numel(ref.z_rs));
out = struct();
for k = 1:size(ref.hills,1)
    terrain.x = xh; terrain.h = ref.hills(k,1)*exp(-0.5*((xh-500)/ref.hills(k,2)).^2);
    terrain.smooth = 0;
    [H,Hp,Hpf] = terrain_profile(terrain,x);
    Phi = main(U0,gamma_x,tau_bar1,k0,delta_x,delta_z,nb_x,nb_z,eps0,eps1,0,0.8,10,Hpf);
    pp = pressure(Phi,Mx,k0,delta_x,delta_z,x,Hp);
    Hr = interp1(x,H,xr);
    for j = 1:numel(ref.z_rs)
        z_r = ref.z_rs(j);
        L_hill = 20*log10(interp2(x,zeta,pp,xr,z_r*ones(size(xr))));
        L_fl   = 20*log10(interp2(x,zeta,p_flat,xr,Hr+z_r));
        IL_pe = L_hill - L_fl;
        IL_bem = squeeze(IL_bem_all(k,j,:)).';
        d = IL_pe - IL_bem;
        ok = IL_bem > -30 & xr <= x(end);  % ignore deep interference nulls, stay inside the PE grid
        res(k,j) = sqrt(mean(d(ok).^2));
        printf('%s, z_r = %2d m: rms IL difference %.2f dB (BEM self-convergence %.2f dB)\n', ...
            ref.names{k}, z_r, res(k,j), sqrt(mean(squeeze(IL_bem_err(k,j,ok)).^2)));
        out(k,j).IL_pe = IL_pe; out(k,j).IL_bem = IL_bem; out(k,j).H = Hr;
        out(k,j).name = ref.names{k}; out(k,j).z_r = z_r;
    end
end
save('-v7','tests/bem_results.mat','xr','out','res');
assert(all(res(:) < 1.0), 'terrain-following WAPE and BEM disagree');
disp('test_terrain_vs_bem: PASSED (rms insertion-loss difference < 1 dB)');
