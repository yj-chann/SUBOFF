function K = Get_K(x)
% Get XIELV OF SUBOFF HULL
[sx,sy]=size(x);
K=zeros(sx,sy);
A=1.126395101;B=0.442874707;
RATIO=0.2;
R_MAX = 5/6;
alpha_Minus_One=1/2.1-1;
x=x/RATIO/0.3048;
for j=1:sy
    for i=1:sx
        xij=x(i,j);
        if xij<=10/3
            K(i,j)= 1./2.1.*R_MAX .* (A.*xij.*(0.3.*xij-1).^4 + B.*xij.^2.*(0.3.*xij-1).^3 + 1 - (0.3.*xij-1).^4.*(1.2.*xij+1)).^alpha_Minus_One...
                .*(4.*A.*xij.*(0.3.*xij-1).^3.*0.3+A.*(0.3.*xij-1).^4+ 3.*B.*xij.^2.*(0.3.*xij-1).^2.*0.3 +2.*B.*xij.*(0.3.*xij-1).^3 - 4.*(0.3.*xij-1).^3.*0.3.*(1.2.*xij+1)-1.2.*(0.3.*xij-1).^4);
        else
            K(i,j)=0;
        end
    end
end
end