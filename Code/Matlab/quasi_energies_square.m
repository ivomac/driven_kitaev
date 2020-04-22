N = 10;

Jy = 0.8;
Jx = 1-Jy;

g = Jx+Jy;
d = Jy-Jx;

T = 1;
B0 = 0;
B1 = linspace(0.001,8*pi,2000);

a = 0.5;

En = zeros(2*N,numel(B1));

for n = 1:numel(B1)

    H0 = maj_H(Jx,Jy,B0,N);
    
    H1 = maj_H(Jx,Jy,B1(n),N);
    
    U = expm(H1* T * (1-a))*expm(H0 * T *a);
    
    [~,ex] = eig(U,'vector');
    
    En(:,n) = angle(ex)/T;
    
end

En = sort(En);

En = table_ordering4(En);

fig = figure; ax = axes; hold on;

plot(B1,En,'.','MarkerSize',2);

axis([0 8*pi/T -pi/T pi/T])

xlabel('$\mu_1$')
ylabel('$\varepsilon$','Interpreter','latex')

std_fig(fig,ax)

%title(...
%['$N=' num2str(N) ', J_x=' num2str(Jx) ', J_y=' num2str(Jy) ...
%', \mu_0=' num2str(B0,'%1.1f') ', T=' num2str(T) '$'] ...
%,'Interpreter','latex','FontSize',9) 
%', \gamma=' num2str(g) ', \Delta=' num2str(d) ...

title(...
['$(\Delta,\mu_0) = (' num2str(d) ', ' num2str(B0,'%1.1f') ')$'] ...
,'Interpreter','latex','FontSize',9) 

ax.XTick = [0 1 2 3 4]*pi/T;
ax.XTickLabel = {'$0$','$\frac{\pi}{T}$','$\frac{2\pi}{T}$','$\frac{3\pi}{T}$','$\frac{4\pi}{T}$'};

ax.YTick = [-pi 0 pi]/T;
ax.YTickLabel = {'$-\frac{\pi}{T}$','$0$','$\frac{\pi}{T}$'};

%export_plot()
