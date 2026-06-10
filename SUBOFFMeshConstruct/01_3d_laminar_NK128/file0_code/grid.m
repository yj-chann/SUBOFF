clc;clear;
% L=0.8712 R=0.0508

NJ1=150; NJ2=450; NK=640; 
NL=NK/4; NJ=NJ1+NJ2;NJ_ADD=200;
NI=200;NI_ADD=20;
% ----------------------%
NJ1=NJ1/5; NJ2=NJ2/5; NK=NK/5; 
NL=NL/5; NJ=NJ/5;NJ_ADD=NJ_ADD/5;
NI=NI/5;NI_ADD=NI_ADD/5;


NN=NI+NI_ADD;

YL=0.55;
% Define the sub-region parameters
YlenRatios  = [0.005 , 0.18 , 1 , 5.75, 10];
YcellRatios = [10 , 70 , 70 , 50 , 20];
YexpRatios  = [1 , 10.7, 3.9, 10.9 , 1];
Y= splitEdge(YL, NN, YlenRatios, YcellRatios, YexpRatios);

% Total Length and Total Cells
SL = 0.5;
NS = NL/2 + NJ + NJ_ADD;

% Define the sub-region parameters
SlenRatios  = [0.35 , 0.95 , 9 , 8];
ScellRatios = [NL/2 , NJ1 , NJ2, NJ_ADD];
SexpRatios  = [1, 1.1 , 4 , 1];

% 子午面
% Generate the mesh nodes
S = splitEdge(SL, NS, SlenRatios, ScellRatios, SexpRatios);
X = Get_x_from_S(S);

% RECT_Y(J,K) RECT_Z(J,K) RECT_R(J,K)
RECT_L = 2 * Get_R(X(NL/2+1));
RECT_Y = transpose(linspace(-0.5*RECT_L,0.5*RECT_L,NL+1))*ones(1,NL+1);
RECT_Z = ones(NL+1,1) * linspace(-0.5*RECT_L,0.5*RECT_L,NL+1);
RECT_Y(1+NL/2,:)=0.0;RECT_Z(:,1+NL/2)=0.0;
RECT_R=sqrt(RECT_Y.^2+RECT_Z.^2);
RECT_X=Get_x_from_R(RECT_R);

fid = fopen('Tecplot_InputFiles/RECT_surface.plt', 'w');
fprintf(fid, 'TITLE = "RECT_Wall"\n');
fprintf(fid, 'VARIABLES = "X", "Y", "Z"\n');
fprintf(fid, 'ZONE T="RECT_Wall", I=%d, J=%d, F=POINT\n', NL+1, NL+1);
for j = 1:NL+1
    for k = 1:NL+1
        fprintf(fid, '%22.15E %22.15E %22.15E\n', RECT_X(j,k) , RECT_Y(j,k) , RECT_Z(j,k) );
    end
end
fclose(fid);


% CY1
S0_L=zeros(1,NK/8);
for i=1:NK/4
S0_L(i)=Get_S(RECT_X(NL+1,i));
end

SREF_L1=S(NL/2+1);
SREF_L2=S(NL/2+NJ1+1);

S1_S=zeros(NK/8,NJ1+1);
for i=1:NK/8
S1_S(i,:)=(SREF_L2-S0_L(i)) * (S(NL/2+1 : NL/2+NJ1+1)- SREF_L1)/(SREF_L2-SREF_L1) +S0_L(i);
end

CY_X1 = zeros(NK,NJ1+1);
CY_X1(1:NK/8,:) = Get_x_from_S(S1_S);
CY_X1(NL/2+1,:) = X(NL/2+1:NL/2+1+NJ1);
CY_X1(NL:-1:NL/2+2,:)=CY_X1(2:NK/8,:);
CY_X1(NL+1:2*NL,:)=CY_X1(1:NL,:);
CY_X1(2*NL+1:3*NL,:)=CY_X1(1:NL,:);
CY_X1(3*NL+1:4*NL,:)=CY_X1(1:NL,:);


THETA_REF=linspace(pi/4,0,NL/2+1);
THETA_REF=THETA_REF(1:NL/2);
THETA_RECT=acos(RECT_Y(NL+1,1:NL/2)./sqrt(RECT_Y(NL+1,1:NL/2).^2+RECT_Z(NL+1,1:NL/2).^2));
THETA_DIST=THETA_REF-THETA_RECT;
CY_Y1 =zeros(NK,NJ1+1);
CY_Z1 =zeros(NK,NJ1+1);
for i=1:NK/8
    THETA_DIST_X=THETA_RECT(i)+(CY_X1(i,:)-CY_X1(i,1))/(CY_X1(i,NJ1+1)-CY_X1(i,1))*THETA_DIST(i);
    CY_Y1(i,:)=Get_R(CY_X1(i,:)).*cos(THETA_DIST_X);
    CY_Z1(i,:)=-Get_R(CY_X1(i,:)).*sin(THETA_DIST_X);
end
CY_Y1(NL/2+1,:)=Get_R(X(NL/2+1:NL/2+NJ1+1));
CY_Z1(NL/2+1,:)=0;
CY_Y1(NL:-1:NL/2+2,:)=CY_Y1(2:NK/8,:);
CY_Z1(NL:-1:NL/2+2,:)=-CY_Z1(2:NK/8,:);
CY_Y1(NL+1:2*NL,:)=-CY_Z1(1:NL,:);
CY_Z1(NL+1:2*NL,:)=CY_Y1(1:NL,:);
CY_Y1(2*NL+1:3*NL,:)=-CY_Y1(1:NL,:);
CY_Z1(2*NL+1:3*NL,:)=-CY_Z1(1:NL,:);
CY_Y1(3*NL+1:4*NL,:)=CY_Z1(1:NL,:);
CY_Z1(3*NL+1:4*NL,:)=-CY_Y1(1:NL,:);

fid = fopen('Tecplot_InputFiles/CY1_surface.plt', 'w');
fprintf(fid, 'TITLE = "CY1_Wall"\n');
fprintf(fid, 'VARIABLES = "X", "Y", "Z"\n');
fprintf(fid, 'ZONE T="CY1_Wall", I=%d, J=%d, F=POINT\n', NK+1, NJ1+1);
for j = 1:NJ1+1
    for k = 1:NK
        fprintf(fid, '%22.15E %22.15E %22.15E\n', CY_X1(k,j) , CY_Y1(k,j) , CY_Z1(k,j) );
    end
    fprintf(fid, '%22.15E %22.15E %22.15E\n', CY_X1(1,j) , CY_Y1(1,j) , CY_Z1(1,j) );
end
fclose(fid);


% CY2
THETA_REF2=linspace(-pi/4,7*pi/4,NK+1);
THETA_REF2=THETA_REF2(1:NK);
CY_X2=ones(NK,1)*X(NL/2+1+NJ1:NL/2+1+NJ);
CY_Y2=(ones(NK,1)*Get_R(X(NL/2+1+NJ1:NL/2+1+NJ))).*(transpose(cos(THETA_REF2))*ones(1,NJ2+1));
CY_Z2=(ones(NK,1)*Get_R(X(NL/2+1+NJ1:NL/2+1+NJ))).*(transpose(sin(THETA_REF2))*ones(1,NJ2+1));


fid = fopen('Tecplot_InputFiles/CY2_surface.plt', 'w');
fprintf(fid, 'TITLE = "CY2_Wall"\n');
fprintf(fid, 'VARIABLES = "X", "Y", "Z"\n');
fprintf(fid, 'ZONE T="CY2_Wall", I=%d, J=%d, F=POINT\n', NK+1, NJ2+1);
for j = 1:NJ2+1
    for k = 1:NK
        fprintf(fid, '%22.15E %22.15E %22.15E\n', CY_X2(k,j) , CY_Y2(k,j) , CY_Z2(k,j) );
    end
    fprintf(fid, '%22.15E %22.15E %22.15E\n', CY_X2(1,j) , CY_Y2(1,j) , CY_Z2(1,j) );
end
fclose(fid);

CY_X=[CY_X1(:,1:NJ1) CY_X2];
CY_Y=[CY_Y1(:,1:NJ1) CY_Y2];
CY_Z=[CY_Z1(:,1:NJ1) CY_Z2];
fid = fopen('Tecplot_InputFiles/CY_surface.plt', 'w');
fprintf(fid, 'TITLE = "CY_Wall"\n');
fprintf(fid, 'VARIABLES = "X", "Y", "Z"\n');
fprintf(fid, 'ZONE T="CY_Wall", I=%d, J=%d, F=POINT\n', NK+1, NJ+1);
for j = 1:NJ+1
    for k = 1:NK
        fprintf(fid, '%22.15E %22.15E %22.15E\n', CY_X(k,j) , CY_Y(k,j) , CY_Z(k,j) );
    end
    fprintf(fid, '%22.15E %22.15E %22.15E\n', CY_X(1,j) , CY_Y(1,j) , CY_Z(1,j) );
end
fclose(fid);



% Calculate total points
P_NUM = (NI+1)*(NJ+1)*NK + (NL-1)^2 * (NI+1);

% Open file for writing
fid = fopen('points', 'w');

% Write OpenFOAM Header
fprintf(fid, '/*--------------------------------*- C++ -*----------------------------------*\\\n');
fprintf(fid, '  =========                 |\n');
fprintf(fid, '  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox\n');
fprintf(fid, '   \\\\    /   O peration     | Website:  https://openfoam.org\n');
fprintf(fid, '    \\\\  /    A nd           | Version:  v2312\n');
fprintf(fid, '     \\\\/     M anipulation  |\n');
fprintf(fid, '\\*---------------------------------------------------------------------------*/\n');
fprintf(fid, 'FoamFile\n');
fprintf(fid, '{\n');
fprintf(fid, '    version     2.0;\n');
fprintf(fid, '    format      ascii;\n');
fprintf(fid, '    class       vectorField;\n');
fprintf(fid, '    location    "constant/polyMesh";\n');
fprintf(fid, '    object      points;\n');
fprintf(fid, '}\n');
fprintf(fid, '// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n');
fprintf(fid, '\n\n');

% Write total points and opening parenthesis
fprintf(fid, '%d\n', P_NUM);
fprintf(fid, '(\n');

% First main loop block
for K = 1:(NL*4)
    for J = 1:(NJ+1)
        for I = NI+1:-1:1
            % Assuming GET_K outputs XIELV and GET_P outputs XP, YP, ZP.
            XIELV = Get_K(CY_X(K,J));
            [XP, YP, ZP] = Get_P(CY_X(K,J), CY_Y(K,J), CY_Z(K,J), XIELV, Y(I));
            
            fprintf(fid, '( %22.15E %22.15E %22.15E )\n', XP, YP, ZP);
        end
    end
end

% Second main loop block
for J = NL:-1:2
    for K = 2:NL
        if J == (NL/2 + 1) && K == (NL/2 + 1)
            for I = NI+1:-1:1
                XP = -Y(I);
                YP = 0.0;
                ZP = 0.0;

                fprintf(fid, '( %22.15E %22.15E %22.15E )\n', XP, YP, ZP);
            end
        else
            for I = NI+1:-1:1
                XIELV = Get_K(RECT_X(J,K));
                [XP, YP, ZP] = Get_P(RECT_X(J,K), RECT_Y(J,K), RECT_Z(J,K), XIELV, Y(I));
                
                fprintf(fid, '( %22.15E %22.15E %22.15E )\n', XP, YP, ZP);
            end
        end
    end
end

% Write closing and footer
fprintf(fid, ')\n\n\n');
fprintf(fid, '// ************************************************************************* //\n');
% Close the file
fclose(fid);



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



