function fv = firm_thresholding(v,T,t)

[m,n]=size(v);
mn = m*n;
fv = zeros(m,n);
for i = 1:mn
    if abs(v(i))<=T
        fv(i) = 0;
    elseif abs(v(i)) > T && abs(v(i)) <=t
        fv(i) = sign(v(i))*((t*(abs(v(i))-T))/(t-T));
    elseif abs(v(i)) > t
        fv(i) = v(i);
    end
end
