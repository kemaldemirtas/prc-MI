% clc
% clear all
% close all




Trajectory
atom_num = 5029;
xyz = readdcd('6wx4 100ns.dcd', 1:atom_num); %read trajectory file

x=xyz(:,1:3:end); %
y=xyz(:,2:3:end); 
z=xyz(:,3:3:end);

 %% Crystal Structure %required to extract coordinates of alpha carbons in trajectory data
PDB = pdbread('6xw4.pdb');
PDB_CA = strcmp({PDB.Model.Atom.AtomName},'CA');

X_coords = {PDB.Model.Atom.X};
Y_coords = {PDB.Model.Atom.Y};
Z_coords = {PDB.Model.Atom.Z};

X = {X_coords{PDB_CA}};
Y = {Y_coords{PDB_CA}};
Z = {Z_coords{PDB_CA}};

%% Extracting C-alphas from dcd file
i = 1;
j = 1;
for i = 1:length(PDB_CA)
    if PDB_CA(i) == 1
       x_CA(:,j) = x(:,i-6);
       y_CA(:,j) = y(:,i-6);
       z_CA(:,j) = z(:,i-6);
       j = j + 1;
    end
    i = i + 1;
end

%%  Mean Structure
for i = 1:width(x_CA)
   X_m(:,i) = sum(x_CA(:,i)) / length(x_CA); 
   Y_m(:,i) = sum(y_CA(:,i)) / length(x_CA);
   Z_m(:,i) = sum(z_CA(:,i)) / length(x_CA);
end

%% Calculation of Deviation
for i = 1:size(x_CA,2)
    x_dev(:,i) = x_CA(:,i) - X_m(i);
    y_dev(:,i) = y_CA(:,i) - Y_m(i);
    z_dev(:,i) = z_CA(:,i) - Z_m(i);
end



