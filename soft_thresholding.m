function sv = soft_thresholding(v,T)

[m,n]=size(v);
mn = m*n;
sv = zeros(m,n);
for i = 1:mn
    if abs(v(i))<=T
        sv(i) = 0;
    else
        sv(i) = sign(v(i))*(abs(v(i))-T);
    end
end
