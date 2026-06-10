% =========================================================================
NN = 220;
NS = 880;
NN_PLUS_ONE = NN+1;
NS_PLUS_ONE = NS+1 ;
% =========================================================================

% Calculate total points
P_NUM = NN_PLUS_ONE * NS_PLUS_ONE + NN_PLUS_ONE * NS;

% --- Read 2D mesh data ---
data=readmatrix('suboff_mesh_2d_ADD.plt','FileType','text','NumHeaderLines',1);


X0 = reshape(data(:,1), NS_PLUS_ONE, NN_PLUS_ONE)';
Y0 = reshape(data(:,2), NS_PLUS_ONE, NN_PLUS_ONE)';

% --- Write OpenFOAM points file ---
% Ensure the output directory exists
if ~exist('../polyMesh', 'dir')
    mkdir('../polyMesh');
end

fid2 = fopen('../polyMesh/points', 'w');
if fid2 == -1
    error('Cannot open ../polyMesh/points for writing.');
end

% Write OpenFOAM Header (escaping backslashes properly for MATLAB)
fprintf(fid2, '/*--------------------------------*- C++ -*----------------------------------*\\\n');
fprintf(fid2, '  =========                 |\n');
fprintf(fid2, '  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox\n');
fprintf(fid2, '   \\\\    /   O peration     | Website:  https://openfoam.org\n');
fprintf(fid2, '    \\\\  /    A nd           | Version:  v2312\n');
fprintf(fid2, '     \\\\/     M anipulation  |\n');
fprintf(fid2, '\\*---------------------------------------------------------------------------*/\n');
fprintf(fid2, 'FoamFile\n');
fprintf(fid2, '{\n');
fprintf(fid2, '    version     2.0;\n');
fprintf(fid2, '    format      ascii;\n');
fprintf(fid2, '    class       vectorField;\n');
fprintf(fid2, '    location    "constant/polyMesh";\n');
fprintf(fid2, '    object      points;\n');
fprintf(fid2, '}\n');
fprintf(fid2, '// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n');
fprintf(fid2, ' \n');
fprintf(fid2, ' \n');
fprintf(fid2, '%d\n', P_NUM);
fprintf(fid2, '(\n');

% Calculate angles
theta1 = 0.5 * 2.0 * pi / NK;
theta2 = (NK - 0.5) * 2.0 * pi / NK;  % 0~360

formatSpec = '( %22.15E %22.15E %22.15E )\n';

% Loop 1: J = 1
j = 1;
for i = 1:NN_PLUS_ONE
    x = X0(i, j);
    y = 0.0;
    z = 0.0;
    fprintf(fid2, formatSpec, x, y, z);
end

% Loop 2: J = 2 to NS_PLUS_ONE (theta1 rotation)
for j = 2:NS_PLUS_ONE
    for i = 1:NN_PLUS_ONE
        x = X0(i, j);
        y = Y0(i, j) * cos(theta1);
        z = -Y0(i, j) * sin(theta1);
        fprintf(fid2, formatSpec, x, y, z);
    end
end

% Loop 3: J = 2 to NS_PLUS_ONE (theta2 rotation)
for j = 2:NS_PLUS_ONE
    for i = 1:NN_PLUS_ONE
        x = X0(i, j);
        y = Y0(i, j) * cos(theta2);
        z = -Y0(i, j) * sin(theta2);
        fprintf(fid2, formatSpec, x, y, z);
    end
end

% Write Footer
fprintf(fid2, ')\n');
fprintf(fid2, ' \n');
fprintf(fid2, ' \n');
fprintf(fid2, '// ************************************************************************* //\n');

fclose(fid2);

disp('Successfully exported OpenFOAM points file.');