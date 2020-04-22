N = 10;

Jx = 0.6;
Jy = 0.4;

T = 1;
B0 = linspace(0.0001,1,300);
B1 = linspace(0.0001,1,30);

fig = figure; ax = axes; hold on;

for n = 1:numel(B1)
    
    En = zeros(2*N,numel(B0));
    
    for m = 1:numel(B0)
        
        H0 = maj_H(Jx,Jy,B1(n),N);
        H1 = maj_H(0 ,0 ,B0(m),N);
        
        U = expm(2*H1)*expm(2*H0*T);
        
        [~,ex] = eig(U,'vector');
        
        En(:,m) = angle(ex);
        
    end
    
    En = sort(En);
    
    En = En([N N+1],:);
    
    En = table_ordering4(En);
    
    ind=find(En(:,end)>0);
    
    zer = interp1(En(ind,:),B0,0);
    
    plot3(zer,B1(n)*ones(size(zer)),zeros(size(zer)),'*r','MarkerSize',4); 
    
    plot3(B0,B1(n)*ones(size(B0)),En(ind,:),'.k','MarkerSize',3);  
    
end

xlabel('$\Gamma_0$','Interpreter','latex','FontSize',16)
ylabel('$\Gamma_1$','Interpreter','latex','FontSize',16)
zlabel('$\varepsilon$','Interpreter','latex','FontSize',16)

fig.Color = 'w';
fig.Renderer='painters';
fig.Position = [300 300 360 200];
ax.TickDir = 'in';
ax.YLabel.Rotation = 0;
box on
grid on
grid minor

%%
%export_fig /home/iaguiar/Dropbox/PhD-1-Edge_spin_coherence/Paper_edge_spin_coherence/Figures/Fig2.png -png -r280
%export_fig C:\Dropbox\PhD-1-Edge_spin_coherence\Paper_edge_spin_coherence\v4\Figures\Fig2.png -png -r280
