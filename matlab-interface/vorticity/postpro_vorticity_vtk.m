%% This programm calculates omega and saves it as a vtk paraview file
% it can also be easily modified to just save the velocity etc.
% Start measure of time
tic

%% Clear Command Window and Workspace
%  -------------------------------------------------------------------
clc
clear variables

%% Add subforlder ./base/ into MATLAB search paths
%  -------------------------------------------------------------------
addpath('./../base/');
addpath('./../plot/');
%% Provide path of the data

%Foldername
% filepath_short='Re28000/nx-256-ny-500-nz-512-alfa0-4.0/';
% filepath_short='Re3000/nx-32-ny-128-nz-128-alfa0-alpha0-8.0-omega-0.5/';
filepath_short='Re3000/nx-32-ny-128-nz-128-alfa0-alpha0-2.0/';
% filepath_short='Re3000/nx-32-ny-128-nz-128-alfa0-alpha0-8.0-omega-2/';
% filepath_short='Re10000/nx-256-ny-500-nz-512-alfa0-4.0/'

% Full path
filepath=['../../../DATA/',filepath_short];
% Get all names of files in the destination
filenamelist=dir([filepath,'/vfield*.dat']);
% sort them (time)
filenamelist= sort_nat({filenamelist.name}.');

%% Set temporal limits of the upcoming evaluation
iF_start=1555;
iF_end=1555%numel(filenamelist);
% spacing between 2 evaluated timesteps
iF_space=2;
size_t=numel(iF_start:iF_space:iF_end);

%% Read field header, set up coordinates and read initial field
%  -------------------------------------------------------------------
iF=iF_start;
filename = strcat(filepath,filenamelist{iF});
[dns,field] = read_header(filename);
tic
setup_derivatives
disp(['setup_derivatives took ',num2str(toc)]);

for iF=iF_start:iF_space:iF_end
    disp(['read field  ',filenamelist{iF}]);
    filename = strcat([filepath,filenamelist{iF}]);
    tic
    [field] = read_field_alternate(filename,dns,field,dc);
    disp(['read field took ',num2str(toc)]);
    tic
    [field] = velocity_gradient(dns, field, derivatives, dc)
    [field] = vorticity(dns, field, derivatives, dc)
    disp(['Vorticity calc took ',num2str(toc)]);
    [field]=rmfield(field,'dU');
    % [field]=rmfield(field,'V');
    [field] = write_vtk_binary(dns, field,{'omega','V'},'Omega',filepath_short,filenamelist{iF},[0.3,1])
end
