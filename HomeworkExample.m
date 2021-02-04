rho0 = 1024; delrho=4; h=25; mu=5;
rho=@(z) (rho0-1/2*delrho*(1+tanh((z+h)/mu)));
z = linspace(-100,0,100)';
im = InternalModes(rho,[-100 0],z,0,'method','spectral');
im.normalization = Normalization.wMax;
im.ShowLowestModesAtWavenumber(10^-6);