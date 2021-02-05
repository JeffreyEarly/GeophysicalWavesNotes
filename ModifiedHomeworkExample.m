rho0 = 1024; delrho=4; h=25; mu=5;
rho=@(z) (rho0-1/2*delrho*(1+tanh((z+h)/mu)));
z = linspace(-100,0,100)';
im = InternalModes(rho,[-100 0],z,0,'method','spectral');
im.normalization = Normalization.wMax;

Nmax = max(sqrt(im.N2));
omega = 0.85*Nmax;

[F,G,h] = im.ModesAtFrequency(omega);

figure
subplot(1,3,1)
plot(F(:,1:4),z, 'LineWidth', 2)
ylabel('depth (meters)');
xlabel('(u,v)-modes');

b = subplot(1,3,2);
plot(G(:,1:4),z, 'LineWidth', 2)
title(b, sprintf('Internal Modes h = (%.2g, %.2g, %.2g, %.2g)', h(1) , h(2), h(3), h(4) ));
xlabel('w-modes');
yticks([]);

subplot(1,3,3)
plot(sqrt(im.N2),im.z, 'LineWidth', 2), hold on
plot([omega omega],[min(z) max(z)],'k','LineWidth', 2)
xlim([0.0 1.1*max(sqrt(im.N2))])
xlabel('buoyancy frequency');

yticks([]);

print('ModesHighOmega.eps','-depsc2')