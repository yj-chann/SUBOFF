function [XP,YP,ZP]=Get_P(XH,YH,ZH,K,DY)
R=sqrt(YH*YH+ZH*ZH);
XP=XH - DY*K/sqrt(1.0+K^2);
YP=( R + DY/sqrt(1.0+K^2) )*YH/R;
ZP=( R + DY/sqrt(1.0+K^2) )*ZH/R;
end