clc;clear;
load('Mesh.mat')
NS=980;NN=440;
filename='../../01_2d_laminar_ADD/10000/U';
fid = fopen(filename, 'r');
dataU=textscan(fid, '(%f %f %*f)', 'HeaderLines',23);
fclose(fid);
Ux=permute(flipud(reshape(dataU{1},NN,NS)),[2,1]);
Uy=permute(flipud(reshape(dataU{2},NN,NS)),[2,1]);

filename='../../01_2d_laminar_ADD/10000/p';
fid = fopen(filename, 'r');
datap=textscan(fid, '%f', 'HeaderLines',23);
fclose(fid);
p=permute(flipud(reshape(datap{1},NN,NS)),[2,1]);

Us=Ux.*(COSPHI'*ones(1,NN))+Uy.*(SINPHI'*ones(1,NN));
Un=Ux.*(-SINPHI'*ones(1,NN))+Uy.*(COSPHI'*ones(1,NN));
flow_field.Us=Us;
flow_field.Un=Un;
flow_field.p=p;
save('flow_field.mat',"flow_field")