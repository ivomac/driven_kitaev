function M = M_mat(B,Jx,Jy,N,R)

if R==0
    
    l1 = repmat([B,Jx]',[N,1]);
    l2 = repmat([-Jy,0]',[N,1]);
    M = diag( l2(1:end-3),3) + diag( l1(1:end-1) ,1);
    
    M = M-M';
    
elseif R==1

    l1 = [B*(2*rand(1,N)-1) ; Jx*(2*rand(1,N)-1)];
    l1 = l1(:);

    l2 = [-Jy*(2*rand(1,N)-1) ; zeros(1,N)];
    l2 = l2(:);
    
    M = diag(l2(1:end-3),3)+diag(l1(1:end-1),1);
    
    M = M-M';
    
    
end