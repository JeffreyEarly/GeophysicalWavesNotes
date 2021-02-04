% This is a quick little demo to show how you can add up functions in a
% complete basis to represent any other function of your choosing.
% Jeffrey J. Early, February 2021

N = 100; % number of points
L = 200; % length of domain

f = @(z) exp(-(z-L/2).^2/(10^2)); % function we are going to represent

% Cosine basis
z_cosine = L*(0:(N-1))'/(N-1);
m = (pi/L)*(0:(N-1));
T_cosine = cos(m.*z_cosine); % we could normalize this by (2/L)
a_cosine = T_cosine\f(z_cosine);
% figure, plot(T_cosine(:,1:4),z_cosine)

% Chebyshev basis
z_cheb = (L/2)*( cos(((0:N-1)')*pi/(N-1)) + 1);
T_cheb = cos( (0:(N-1)) .* ((0:N-1)')*pi/(N-1) );
a_cheb = T_cheb\f(z_cheb);
% figure, plot(T_cheb(:,1:4),z_cheb)

% Internal mode basis---let's use a goofy stratification profile
k = 0;
[rhoFunc, ~, zIn] = InternalModes.StratificationProfileWithName('pycnocline-constant');
z = linspace(min(zIn),max(zIn),1024)';
im = InternalModesSpectral(rhoFunc,zIn,z,33,'nEVP',512);
% Now find the guass-quadrature points, and use those instead.
z_im = im.GaussQuadraturePointsForModesAtWavenumber(N,k);
im_gauss = InternalModesSpectral(rhoFunc,zIn,z_im,33,'nEVP',512,'nModes',N);
[F,G,h] = im_gauss.ModesAtWavenumber(0);
G = G(:,1:(N-2)); % Last modes are not linearly independent
a_im = G\f(z_im+L);
% figure, plot(G(:,1:4),z)

figure('Position',[1 1 600 200])
N_modes = 4;
subplot(1,N_modes+1,1)
plot(f(z_cosine),z_cosine, 'k', 'LineWidth', 2), set(gca, 'XTickLabel', []), set(gca, 'YTickLabel', [])
for iMode=1:N_modes
    subplot(1,N_modes+1,1+iMode)
    plot(T_cosine(:,iMode),z_cosine, 'k', 'LineWidth', 2), set(gca, 'XTickLabel', []), set(gca, 'YTickLabel', [])
    title(sprintf('a=%.2g',a_cosine(iMode)));
end
print('cosine-representation.eps','-depsc2')

figure('Position',[1 1 600 200])
subplot(1,N_modes+1,1)
plot(f(z_cheb),z_cheb, 'k', 'LineWidth', 2), set(gca, 'XTickLabel', []), set(gca, 'YTickLabel', [])
for iMode=1:N_modes
    subplot(1,N_modes+1,1+iMode)
    plot(T_cheb(:,iMode),z_cheb, 'k', 'LineWidth', 2), set(gca, 'XTickLabel', []), set(gca, 'YTickLabel', [])
    title(sprintf('a=%.2g',a_cheb(iMode)));
end
print('chebyshev-representation.eps','-depsc2')

figure('Position',[1 1 600 200])
subplot(1,N_modes+1,1)
plot(f(z_im+L),z_im, 'k', 'LineWidth', 2), set(gca, 'XTickLabel', []), set(gca, 'YTickLabel', [])
for iMode=1:N_modes
    subplot(1,N_modes+1,1+iMode)
    plot(G(:,iMode),z_im, 'k', 'LineWidth', 2), set(gca, 'XTickLabel', []), set(gca, 'YTickLabel', [])
    title(sprintf('a=%.2g',a_im(iMode)));
end
print('internal-mode-representation.eps','-depsc2')

figure
plot(a_cosine.^2, 'LineWidth', 1.5), hold on
plot(a_cheb.^2, 'LineWidth', 1.5)
plot(a_im.^2, 'LineWidth', 1.5)
ylog
legend('cosine', 'chebyshev', 'internal mode')
