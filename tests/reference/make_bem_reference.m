% Generates tests/reference/bem_reference.mat, the converged boundary-element
% reference used by tests/test_terrain_vs_bem.m (insertion loss of three
% rigid Gaussian hills, line source, 100 Hz). Slow (tens of minutes, ~8 GB RAM).
% Run from the repository root:  addpath tests/reference; make_bem_reference

freq = 100; hS = 10; dim_x = 1000;
T = 10; Rg = 286.6896552; gam = 1.41; temp0 = 273.15+T;
c0 = sqrt(gam*Rg*temp0); k0 = 2*pi*freq/c0;
xh = 0:1:dim_x;
names = {'gentle hill  (20 m, sigma 150 m, max slope  4.6 deg)'
         'Gaussian hill (40 m, sigma 100 m, max slope 13.6 deg)'
         'steep hill   (30 m, sigma  50 m, max slope 19.9 deg)'};
hills = [20 150; 40 100; 30 50];        % height, sigma (m), centred at 500 m
z_rs = [2 10];
xr = 150:5:dim_x;
npw = [30 40];                          % elements per wavelength (error estimate, reference)
G = @(x,z) 1i/4*(besselh(0,1,k0*hypot(x,z-hS)) + besselh(0,1,k0*hypot(x,z+hS)));
IL = zeros(numel(names),numel(z_rs),numel(xr),numel(npw));
for k = 1:numel(names)
    hh = hills(k,1)*exp(-0.5*((xh-500)/hills(k,2)).^2);
    keep = hh > 1e-4*max(hh);
    rcv = []; 
    for j = 1:numel(z_rs)
        rcv = [rcv; xr(:) interp1(xh,hh,xr(:))+z_rs(j)];
    end
    for m = 1:numel(npw)
        tic; p = bem_rigid_hill(k0,[0 hS],xh(keep),hh(keep),rcv,npw(m));
        il = 20*log10(abs(p./G(rcv(:,1),rcv(:,2))));
        IL(k,:,:,m) = reshape(il,numel(xr),numel(z_rs)).';
        printf('%s, lambda/%d done in %.0f s\n', names{k}, npw(m), toc); fflush(stdout);
    end
end
save('-v7','tests/reference/bem_reference.mat','freq','hS','hills','names','z_rs','xr','npw','IL');
