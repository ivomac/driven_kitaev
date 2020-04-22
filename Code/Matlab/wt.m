
N = 4;

d = 0.8;
t = 1;

B0 = 1.4;
T = 1;

B = 0.1;

Jy = (t+d)/2;
Jx = (t-d)/2;

En = zeros(2*N,1);

M0 = M_mat(B0,Jx,Jy,N,0);
M1 = M_mat(B-B0,0,0,N,0);

U2 = expm(M1*T)*expm(M0 * T);

%[~,ex] = eig(U,'vector');
%En(:,n) = angle(ex)/T;

