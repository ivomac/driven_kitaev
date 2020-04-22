N = 50;

g = pi;
d = 0.1;
% g = Jx+Jy;
% d = Jy-Jx;

Jy = (g+d)/2;
Jx = (g-d)/2;

T = 9.4*pi;

B0 = 0;
B1 = linspace(0,pi/T,500);

En = zeros(2*N,numel(B1));

for n = 1:numel(B1)

    H0 = maj_H(Jx,Jy,B0,N);

    ex_H0 = expm(H0 * T);

    H1 = maj_H(0,0,B1(n)*T,N);

    U = expm(H1)*ex_H0;

    [~,ex] = eig(U,'vector');

    En(:,n) = angle(ex)/T;

end

En = sort(En);

% En = table_ordering4(En);

fig = figure; ax = axes; hold on;

plot(B1,En,'.','MarkerSize',2);

axis([0 pi/T 0 pi/T])

xlabel('$\mu_1$')
ylabel('$\varepsilon$','Interpreter','latex')

% std_fig(fig,ax)

%title(...
%['$N=' num2str(N) ', J_x=' num2str(Jx) ', J_y=' num2str(Jy) ...
%', \mu_0=' num2str(B0,'%1.1f') ', T=' num2str(T) '$'] ...
%,'Interpreter','latex','FontSize',9)
%', \gamma=' num2str(g) ', \Delta=' num2str(d) ...

title(...
['$(\Delta,\mu_0) = (' num2str(d) ', ' num2str(B0,'%1.1f') ')$'] ...
,'Interpreter','latex','FontSize',9)

ax.XTick = [-1 0 1]*pi/T;
ax.XTickLabel = {'$\frac{-\pi}{T}$','$0$','$\frac{\pi}{T}$'};

ax.YTick = [-1 0 1]*pi/T;
ax.YTickLabel = {'$-\frac{\pi}{T}$','$0$','$\frac{\pi}{T}$'};

%export_plot()
