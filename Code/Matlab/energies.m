N = 10;
B = 0:0.001:2;

Jx = 0.8;
Jy = -0.1;

En = zeros(N,numel(B));

for n = 1:numel(B)
    
    T = B(n)*(Jy+Jx);
    U = Jx*Jy;
    
    A = T*diag(ones(1,N-1),1) + U*diag(ones(1,N-2),2);
    
    A = A + A';
    
    A(1,1) = -Jx^2;
    
    A(N,N) = -Jy^2;
    
    [V,E] = eig(A,'vector');
    
    E = E + B(n)^2 + Jx^2 + Jy^2;
    
    En(:,n) = sort(sqrt(abs(E)));
    
end

En = [-En;En];

En = table_ordering4(En,1);

fig = figure; ax = axes; hold on; 
fig.Color = 'w';
fig.Renderer='painters';
fig.Position = [100 100 300 240];
box on; grid on; grid minor;

plot(B,En,'-','LineWidth',1.2);

xlabel('$\mu_1$','Interpreter','latex','FontSize',12)
ylabel('$\varepsilon$','Interpreter','latex','FontSize',12)

ax.TickDir = 'in';
ax.YLabel.Rotation = 90;
%axis tight
axis([0 pi -pi/2 pi/2])
box on
grid on
grid minor

%%
%export_fig /home/iaguiar/Dropbox/PhD-1-Edge_spin_coherence/Paper_edge_spin_coherence/Figures/Fig2.png -png -r280

export_fig C:\Dropbox\PhD-1-Edge_spin_coherence\Paper_edge_spin_coherence\v4\Figures\Fig2.png -png -r280
