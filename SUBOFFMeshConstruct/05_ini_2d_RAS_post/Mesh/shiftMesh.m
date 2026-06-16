clc;clear;
NN = 440;
NS = 981;
NN_PLUS_ONE = NN+1;
NS_PLUS_ONE = NS+1 ;
NK = 640;
L=4.356*0.2;

data=readmatrix('suboff_mesh_2d_ADD.plt','FileType','text','NumHeaderLines',1);
X0 = reshape(data(:,1), NS_PLUS_ONE, NN_PLUS_ONE)';
Y0 = reshape(data(:,2), NS_PLUS_ONE, NN_PLUS_ONE)';
X0_WALL=X0(NN_PLUS_ONE,:);
X_WALL_CENTER=(X0_WALL(1:NS)+X0_WALL(2:NS_PLUS_ONE))/2;
[kappa0,kappa0_s]=calculate_suboff_curvature(X_WALL_CENTER);
S0_WALL_CENTER=Get_S(X_WALL_CENTER);
COSPHI=(1 + Get_K(X_WALL_CENTER).^2).^(-1/2);
SINPHI=sqrt(abs(1 - COSPHI.^2));
r0=Get_R(X_WALL_CENTER);

% Calculate cell centers by averaging the 4 corners using vectorization
X = (X0(1:end-1, 1:end-1) + ... % Top-left corners
     X0(2:end,   1:end-1) + ... % Bottom-left corners
     X0(1:end-1, 2:end)   + ... % Top-right corners
     X0(2:end,   2:end)) / 4;   % Bottom-right corners

Y = (Y0(1:end-1, 1:end-1) + ...
     Y0(2:end,   1:end-1) + ...
     Y0(1:end-1, 2:end)   + ...
     Y0(2:end,   2:end)) / 4;
X=flipud(X)';
Y=flipud(Y)';
S=S0_WALL_CENTER;
N=flipud(-(X0(1:end-1,1)+X0(2:end,1))/2)';


save('Mesh.mat','X','Y','NN','NS','kappa0','kappa0_s','SINPHI','COSPHI','S','N','r0')