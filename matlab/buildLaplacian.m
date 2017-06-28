function L = buildLaplacian(x,y);

n = x*y;
L = zeros(n,n);

ii = 1;
for j=1:y;
    for i=1:x;
        l = zeros(x,y);
        l(i,j) = -1;
        
        if i>1; i1=i-1; else; i1=i+1; end;
        if j>1; j1=j-1; else; j1=j+1; end;
        if i<x; i2=i+1; else; i2=i-1; end;
        if j<y; j2=j+1; else; j2=j-1; end;
        
        l(i1,j) = l(i1,j)+0.25;
        l(i2,j) = l(i2,j)+0.25;
        if j1>0; l(i,j1) = l(i,j1)+0.25; end;
        if j2>0; l(i,j2) = l(i,j2)+0.25; end;
        
        L(ii,:) = l(:);
        ii = ii+1;
    end;
end;
        