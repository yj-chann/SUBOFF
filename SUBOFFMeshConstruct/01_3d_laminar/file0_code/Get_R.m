function R = Get_R(x)
% SUBOFF HULL in (m)
% RATIO=0.2
[sx,sy]=size(x);
R=zeros(sx,sy);
RATIO=0.2;
R_MAX = 5/6;
x=x/RATIO/0.3048;
for j=1:sy
    for i=1:sx
        xij=x(i,j);
        if xij <= 10/3
            A=1.126395101;B=0.442874707;
            R(i,j) = R_MAX .* (A.*xij.*(0.3.*xij-1).^4 + B.*xij.^2.*(0.3.*xij-1).^3 + 1 - (0.3.*xij-1).^4.*(1.2.*xij+1)).^(1/2.1);
            R(i,j) = R(i,j)*0.3048*RATIO;
        else
            R(i,j)=R_MAX*0.3048*RATIO;
        end
    end
end
end