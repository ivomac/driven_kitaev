function H = maj_H(Jx,Jy,B,N)

v1 = repmat([B -Jx],[1 N]);
v1(end) = [];

v2 = repmat([Jy 0],[1 N]);
v2((end-2):(end)) = [];

H = diag(v1,1)-diag(v1,-1)+diag(v2,3)-diag(v2,-3);
