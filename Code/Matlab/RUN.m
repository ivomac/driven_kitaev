close all;clear all;clc

N = 10;

d = 0.8;
t = 1;

B0 = 1.4;
T = 1;

ord = 2;

k = 0.5;
B = linspace(-k*pi/T,k*pi/T,1000);

quasienergies_exact(N,t,d,B0,B,T)

quasienergies_expansion(N,t,d,B0,B,T,ord);

Path = 'C:\Dropbox\PhD\Floquet_Kitaev\Notes\Figures\New';
export_plots(Path)

function std_fig_quasienergies(fig,ax,d,B0,B,En,T)

En = sort(En);
%En = table_ordering4(En);

plot(ax,B,En,'.');

axis([B(1) B(end) B(1) B(end)])

xlabel('$\mu$')
ylabel('$\varepsilon$','Interpreter','latex')

title(['$(\Delta,\mu_0) = (' num2str(d) ', ' num2str(B0,'%1.1f') ')$'])

ax.XTick = [-1 -0.5 0 0.5 1]*pi/T;
ax.XTickLabel = {'$-\frac{\pi}{T}$','$-\frac{\pi}{2T}$','$0$' ...
    ,'$\frac{\pi}{2T}$','$\frac{\pi}{T}$'};

ax.YTick = ax.XTick;
ax.YTickLabel = ax.XTickLabel;

std_fig(fig,ax)

end

function std_fig_amplitudes(fig,ax,d,B0,B,A,T)

plot(B,A,'-');

ax.YAxis.Exponent = round( log10( max(max(A))-min(min(A)) ) );
axis([B(1) B(end) min(min(A)) max(max(A))])

xlabel('$\mu$')
ylabel('$\Vert H_n \Vert / \Vert H_0 \Vert $','Interpreter','latex')

title(['$(\Delta,\mu_0) = (' num2str(d) ', ' num2str(B0,'%1.1f') ')$']) 
ax.XTick = [-1 -0.5 0 0.5 1]*pi/T;
ax.XTickLabel = {'$-\frac{\pi}{T}$','$-\frac{\pi}{2T}$','$0$' ...
    ,'$\frac{\pi}{2T}$','$\frac{\pi}{T}$'};

std_fig(fig,ax)

txt = text(B0,2e-2*max(max(A)),'$\leftarrow \mu_0$');
txt.Rotation = 90;
txt.Interpreter = 'latex';

legendflex(ax,{num2str((1:4)')}...
    , 'ref', ax, 'xscale', 0.6...
    , 'anchor',{'ne','ne'}, 'buffer', [-10 -10]...
    ,'fontsize', 9,'Interpreter', 'latex'...
    ,'box', 'on', 'title', 'order'); 

end

function quasienergies_exact(N,t,d,B0,B,T)

Jy = (t+d)/2;
Jx = (t-d)/2;

En = zeros(2*N,numel(B));

M0 = M_mat(B0,Jx,Jy,N,0);
ex_M0 = expm(M0 * T);

for n = 1:numel(B)

    M1 = M_mat(B(n)-B0,0,0,N,0);
    
    U = expm(M1*T)*ex_M0;
    
    [~,ex] = eig(U,'vector');
    En(:,n) = angle(ex)/T;
    
end

fig = figure('Name','Quasienergies_Exact'...
    ,'NumberTitle','off'); ax = axes; hold on;
std_fig_quasienergies(fig,ax,d,B0,B,En,T)

end

function quasienergies_expansion(N,t,d,B0,B,T,ord)

Jy = (t+d)/2;
Jx = (t-d)/2;

En = zeros(2*N,numel(B));

M0 = M_mat(B0,Jx,Jy,N,0);

H = cell(1,4);
A = zeros(numel(B),4);

for n = 1:numel(B)

    M1 = M_mat(B(n)-B0,0,0,N,0);
    
    H0 = (M0+M1);
    
    H{1} = C(M0,M1)*T/2;
    H{2} = (C(M0,H{1})-C(M1,H{1}))*T^2/12;
    H{3} = -C(M1,M0,H{1})*T^3/24;
    H{4} = T^4*(-1/720*( C(M0,M0,M0,M0,M1)  +  C(M1,M1,M1,M1,M0) )...
        + 1/360*( C(M1,M0,M0,M0,M1)  +  C(M0,M1,M1,M1,M0) )...
        +1/120*( C(M0,M1,M0,M1,M0)  + C(M1,M0,M1,M0,M1 ) ) ) ;
    
    H_F = H0;
    A0 = sqrt(trace(H0'*H0));
    for m = 1:ord
        H_F = H_F + H{m};
    end
    
    for m = 1:4
        A(n,m) = sqrt(trace(H{m}'*H{m}))/A0;   
    end
    
    H_F = 1i*H_F;

    [~,E] = eig(H_F,'vector');
    E = real(E);
    E = mod(E+pi/T,2*pi/T)-pi/T;
    En(:,n) = E;
    
end

fig1 = figure('Name','Quasienergies_Expansion'...
    ,'NumberTitle','off'); ax = axes; hold on;
std_fig_quasienergies(fig1,ax,d,B0,B,En,T)

fig2 = figure('Name','Expansion_Amplitude'...
    ,'NumberTitle','off'); ax = axes; hold on;
std_fig_amplitudes(fig2,ax,d,B0,B,A,T)

% for n = 1:4
%     figure(n); colormap('gray')
%     imagesc(abs(H{n})) %H{n}*H{n}'
%     title(num2str(n))
% end

end

function quasienergies_D_effective(N,t,d,B0,B,T)

B1 = B - B0;

t_ef = t*(1); %-d^2*B_d.^2*T^4/180
d_ef = d*(1-B1.*B*T^2/3 ); %+B_d.^2*T^4/720.*(6*B_d.^2+t^2-9*d^2)

Jx_ef= (t_ef-d_ef)/2;
Jy_ef= (t_ef+d_ef)/2;

En = zeros(2*N,numel(B));

for n = 1:numel(B)

    H = 1i*M_mat(B(n),Jx_ef(n),Jy_ef(n),N,0);
    
    [~,E] = eig(H,'vector');
    
    E = mod(E+pi/T,2*pi/T)-pi/T;
    
    En(:,n) = E;
    
end

fig = figure('Name','Quasienergies_D_effective'...
    ,'NumberTitle','off'); ax = axes; hold on;
std_fig_quasienergies(fig,ax,d,B0,B,En,T)

%figure; plot(B,d_ef)

end


