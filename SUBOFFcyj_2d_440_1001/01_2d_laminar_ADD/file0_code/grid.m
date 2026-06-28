clc;clear;
% L=0.8712 R=0.0508

NJ1=150; NJ2=450; NK=640; 
NL=NK/4; NJ=NJ1+NJ2;NJ_ADD=300;
NI=400;NI_ADD=40;
% ----------------------%


NN=NI+NI_ADD;

YL=0.55;
% Define the sub-region parameters
YlenRatios  = [0.005 , 0.18 , 1 , 5.75, 10];
YcellRatios = [10 , 70 , 70 , 50 , 20];
YexpRatios  = [1 , 10.7, 3.9, 10.9 , 1];
Y= splitEdge(YL, NN, YlenRatios, YcellRatios, YexpRatios);

% Total Length and Total Cells
SL = 0.65;
NS = NL/2 + NJ + NJ_ADD;

% Define the sub-region parameters
SlenRatios  = [0.35 , 0.95 , 9 , 12];
ScellRatios = [NL/2 , NJ1 , NJ2, NJ_ADD];
SexpRatios  = [1, 1.1 , 4 , 1];



% 子午面
% Generate the mesh nodes
S = splitEdge(SL, NS, SlenRatios, ScellRatios, SexpRatios);
X = Get_x_from_S(S);

figure
plot(S)
title('S')

fid= fopen('Tecplot_InputFiles\suboff_mesh_2d.plt','w');
fprintf(fid,'ZONE T = "suboff_mesh_2d", I = %d, J = %d, F=POINT\n',NL/2+NJ+1,NI+1);
for I=NI+1:-1:1
for J=1:NL/2+NJ+1  
    if J==1
        fprintf(fid, '%22.15E %22.15E\n', -Y(I) , 0.0);
    else
        XIELV = Get_K(X(J));
        [XP, YP, ZP] = Get_P(X(J), Get_R(X(J)), 0, XIELV, Y(I)); 
        fprintf(fid, '%22.15E %22.15E\n', XP, YP);
    end
end
end
fclose(fid);

fid= fopen('Tecplot_InputFiles\suboff_mesh_2d_ADD.plt','w');
fprintf(fid,'ZONE T = "suboff_mesh_2d_ADD", I = %d, J = %d, F=POINT\n',NS+1,NN+1);
for I=NN+1:-1:1
for J=1:NS+1  
        if J==1
            fprintf(fid, '%22.15E %22.15E\n', -Y(I) , 0.0);
        else
            XIELV = Get_K(X(J));
            [XP, YP, ZP] = Get_P(X(J), Get_R(X(J)), 0, XIELV, Y(I));
            fprintf(fid, '%22.15E %22.15E\n', XP, YP);
        end
end
end
fclose(fid);



