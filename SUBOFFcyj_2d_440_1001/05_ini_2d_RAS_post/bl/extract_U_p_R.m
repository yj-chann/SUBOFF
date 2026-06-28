clc;clear;
load('ReadData/Mesh.mat')
NS=981;NN=440;
filename='../../05_ini_2d_RAS/50000/U';
fid = fopen(filename, 'r');
dataU=textscan(fid, '(%f %f %*f)', 'HeaderLines',23);
fclose(fid);
Ux=permute(flipud(reshape(dataU{1},NN,NS)),[2,1]);
Uy=permute(flipud(reshape(dataU{2},NN,NS)),[2,1]);

filename='../../05_ini_2d_RAS/50000/p';
fid = fopen(filename, 'r');
datap=textscan(fid, '%f', 'HeaderLines',23);
fclose(fid);
p=permute(flipud(reshape(datap{1},NN,NS)),[2,1]);

Us=Ux.*(COSPHI'*ones(1,NN))+Uy.*(SINPHI'*ones(1,NN));
Un=Ux.*(-SINPHI'*ones(1,NN))+Uy.*(COSPHI'*ones(1,NN));
flow_field.Us=Us;
flow_field.Un=Un;
flow_field.p=p;


filename='../../05_ini_2d_RAS/50000/R';
fid = fopen(filename, 'r');
% OpenFOAM symmTensor syntax: (Rxx Rxy Rxz Ryy Ryz Rzz)
dataR=textscan(fid, '(%f %f %*f %f %*f %f)', 'HeaderLines',23); 
fclose(fid);

% Extract and reshape the active 2D components
Rxx=permute(flipud(reshape(dataR{1},NN,NS)),[2,1]);
Rxy=permute(flipud(reshape(dataR{2},NN,NS)),[2,1]);
Ryy=permute(flipud(reshape(dataR{3},NN,NS)),[2,1]);
Rzz=permute(flipud(reshape(dataR{4},NN,NS)),[2,1]);


% Transform Reynolds stresses to stream-normal coordinates (Optional but typical)
% Rss = Rxx*cos^2 + Ryy*sin^2 + 2*Rxy*cos*sin
% Rnn = Rxx*sin^2 + Ryy*cos^2 - 2*Rxy*cos*sin
% Rsn = (Ryy - Rxx)*cos*sin + Rxy*(cos^2 - sin^2)
cos_mat = COSPHI'*ones(1,NN);
sin_mat = SINPHI'*ones(1,NN);

Rss = Rxx.*(cos_mat.^2) + Ryy.*(sin_mat.^2) + 2.*Rxy.*(cos_mat.*sin_mat);
Rnn = Rxx.*(sin_mat.^2) + Ryy.*(cos_mat.^2) - 2.*Rxy.*(cos_mat.*sin_mat);
Rsn = (Ryy - Rxx).*(cos_mat.*sin_mat) + Rxy.*(cos_mat.^2 - sin_mat.^2);


flow_field.Rss=Rss;
flow_field.Rnn=Rnn;
flow_field.Rsn=Rsn;
flow_field.Rtt=Rzz;
save('flow_field.mat',"flow_field")

%% 
contourf(X/0.8712,Y/0.8712,Us,200,'LineStyle','none')
axis equal

%%
plot(N,Us(662,:))
hold on
plot(N,Us(66,:))
