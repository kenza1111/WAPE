% Validation of the terrain-following WAPE on a uniform slope.
% Homogeneous, motionless atmosphere above a rigid plane z = s*x (s = slope).
% Exact reference: point source + image source across the tilted plane.
% Run from the repository root:  addpath tests; test_slope_analytic

freq = 100; hS = 30; z_r = 10; dim_x = 1000; dim_z = 200;
T = 10; Rg = 286.6896552; gam = 1.41; temp0 = 273.15+T;
c0 = sqrt(gam*Rg*temp0); lambda = c0/freq; k0 = 2*pi*freq/c0;
delta_x = lambda/10; delta_z = lambda/10;
nb_x = floor(dim_x/delta_x)+1; nb_z = floor(dim_z/delta_z)+1;
z = ((0:nb_z-1)*delta_z)'; x = (0:nb_x-1)*delta_x;
Beta = 0; Rp = 1;
[Mx,eps0,eps1] = Mach(0,80,z,1e-5,0,0.15,c0,0,temp0,Rg,gam,0);
gamma_x = 1./sqrt(1-Mx.^2);
tau_bar1 = (Mx.*gamma_x.^2).*(sqrt(1+eps1) - Mx);

slopes = [0 0.05 0.1 -0.1 0.2 0.3];
iz = round(z_r/delta_z)+1; zr = z(iz);
win = x >= 200;                       % PE not valid in the near field
results = zeros(numel(slopes),3);
for k = 1:numel(slopes)
    s = slopes(k);
    terrain.x = [0 dim_x]; terrain.h = [0 s*dim_x]; terrain.smooth = 0;
    [H,Hp,Hpf] = terrain_profile(terrain,x);
    cG = 1/sqrt(1+Hp(1)^2);   % starter defined perpendicular to the slope at the source
    U0 = sourceSalomons_order2(z*cG,hS*cG,nb_z,k0,Rp);
    Phi = main(U0,gamma_x,tau_bar1,k0,delta_x,delta_z,nb_x,nb_z,eps0,eps1,Beta,0.8,10,Hpf);
    pp = pressure(Phi,Mx,k0,delta_x,delta_z,x,Hp);
    % receiver at height zr above the local ground
    xr = x; zrabs = H + zr;
    R1 = sqrt(xr.^2 + (zrabs-hS).^2);
    dL_pe = 20*log10(pp(iz,:)) + 20*log10(R1) - 10*log10(xr);
    % exact solution: image of S=(0,hS) across the line z = s*x
    nrm = [-s 1]/sqrt(1+s^2); dS = hS/sqrt(1+s^2);
    Si = [0 hS] - 2*dS*nrm;
    R2 = sqrt((xr-Si(1)).^2 + (zrabs-Si(2)).^2);
    dL_ex = 20*log10(abs(1 + (R1./R2).*exp(1i*k0*(R2-R1))));
    d = dL_pe(win) - dL_ex(win);
    lit = dL_ex(win) > -6;           % outside the interference nulls
    results(k,:) = [mean(d) sqrt(mean(d(lit).^2)) sqrt(mean(d.^2))];
    printf('slope %+5.2f (%5.1f deg): mean %+5.2f dB, rms outside nulls %4.2f dB, rms all %4.2f dB\n', ...
        s, atand(s), results(k,1), results(k,2), results(k,3));
    curves{k} = [dL_pe; dL_ex];
end
save('-v7','tests/slope_results.mat','x','slopes','curves','results');

% Outside the interference nulls the error must stay small up to 16.7 deg.
% (In the nulls, a fringe shift of a fraction of a metre gives large dB
% differences; this is reported as "rms all".)
assert(all(results(:,2) < 0.5), 'uniform-slope error too large outside the nulls');
disp('test_slope_analytic: PASSED (rms < 0.5 dB outside interference nulls, slopes up to 16.7 deg)');
