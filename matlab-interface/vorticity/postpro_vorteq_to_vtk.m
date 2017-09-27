%% This function calculates the different terms of the vorticity transport equation and saves them as vtk
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

%% Set temporal limits of the upcoming evaluation
iF_start=1500;
iF_end=1500%numel(filenamelist);
% spacing between 2 evaluated timesteps
iF_space=2;
size_t=numel(iF_start:iF_space:iF_end);
%% Read field header, set up coordinates and read initial field
%  -------------------------------------------------------------------

%Foldername
% filepath_short='Re28000/nx-256-ny-500-nz-512-alfa0-4.0/';
% filepath_short='Re3000/nx-32-ny-128-nz-128-alfa0-alpha0-8.0-omega-0.5/';
filepath_short='Re3000/nx-32-ny-128-nz-128-alfa0-alpha0-2.0/';
% filepath_short='Re3000/nx-32-ny-128-nz-128-alfa0-alpha0-8.0-omega-2/';
% filepath_short='Re10000/nx-256-ny-500-nz-512-alfa0-4.0/'
filepath=['../../../DATA/',filepath_short]
filenamelist=dir([filepath,'/vfield*.dat']);
filenamelist= sort_nat({filenamelist.name}.');
filename = strcat(filepath,filenamelist{iF_start});
[dns,field] = read_header(filename);

tic
setup_derivatives
disp(['setup_derivatives took ',num2str(toc)]);
% [field] = poisson(dns, field, derivatives, dc);
% return
tic

%% Load the fields, calculate the different terms of the equation, save them as vtk
for iF=iF_start:iF_space:iF_end
 disp(['read field  ',filenamelist{iF}]);
 filename = strcat([filepath,filenamelist{iF}]);
[field] = read_field_alternate(filename,dns,field,dc);
disp(['read field took ',num2str(toc)]);
% [field] = continuumequation(dns, field, derivatives, dc)
tic
[field] = vorticitytransport(dns, field, derivatives, dc)
[field] = write_vtk_binary(dns, field,{'omega','SR','CV','DIF'},'Vorteq',filepath_short,filenamelist{iF},[0.3,1])

end
