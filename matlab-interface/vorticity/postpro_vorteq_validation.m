%% This programm validates if the vorticity transport equation is fullfilled
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
%% Check which timestep
%use this exemplary timestep for the validation. can be varied
t_valid=1500;

%% Read field header, set up coordinates and read initial field
%  -------------------------------------------------------------------
filepath_short='Re3000/nx-32-ny-128-nz-128-alfa0-alpha0-2.0/';
filepath=['../../../DATA/',filepath_short]
filenamelist=dir([filepath,'/vfield*.dat']);
filenamelist= sort_nat({filenamelist.name}.');
filename = strcat(filepath,filenamelist{t_valid});
[dns,field] = read_header(filename);
setup_derivatives
disp(['setup_derivatives took ',num2str(toc)]);
% [field] = poisson(dns, field, derivatives, dc);
% return

%% Calculate the temporal derivative of the equation with central differences
% omega(t_valid-1)
iF=t_valid-1;
% get the time
disp(['read field  ',filenamelist{iF}]);
filename = strcat([filepath,filenamelist{iF}]);
f = fopen(filename);
header=fread(f,1024,'*char')';
i=find(header=='=');
fclose(f);
t1 = str2double(header(i(9)+1 : end));
% load the field and calculate vorticity
[field] = read_field_alternate(filename,dns,field,dc);
disp(['read field took ',num2str(toc)]);
[field] = velocity_gradient(dns, field, derivatives, dc)
[field] = vorticity(dns, field, derivatives, dc)
vort1=field.omega;
% omega(t_valid-1)
iF=t_valid+1;
disp(['read field  ',filenamelist{iF}]);
filename = strcat([filepath,filenamelist{iF}]);
f = fopen(filename);
header=fread(f,1024,'*char')';
i=find(header=='=');
fclose(f);
t3 = str2double(header(i(9)+1 : end));
% load the field and calculate vorticity
[field] = read_field_alternate(filename,dns,field,dc);
disp(['read field took ',num2str(toc)]);
[field] = velocity_gradient(dns, field, derivatives, dc)
[field] = vorticity(dns, field, derivatives, dc)
vort3=field.omega;

%temporal derivative of omega
dOmegadt=cell(dns.ny+1,1);
dt=t3-t1;
for iy=1+(0:dns.ny)
    dOmegadt{iy}=(vort3{iy}-vort1{iy})/dt;
end
%% Calculate the right side of the equation
for iF=t_valid
    disp(['read field  ',filenamelist{iF}]);
    filename = strcat([filepath,filenamelist{iF}]);
    [field] = read_field_alternate(filename,dns,field,dc);
    disp(['read field took ',num2str(toc)]);
    tic
    [field] = vorticitytransport(dns, field, derivatives, dc)
end
%% Calculate the residuum
restVorteq=cell(dns.ny+1,1);
for iy=1+(0:dns.ny)
    restVorteq{iy}=field.y(iy)^2*dOmegadt{iy}-squeeze(field.CV{iy})-squeeze(field.DIF{iy})-squeeze(field.SR{iy});
end

% Plot the results
addpath('./../plot/');
[plotit] = total_ift_nat(dns, field,restVorteq, dc);
plotwhat.name{1}='rphi';
plotwhat.plane{1}=19;
plotwhat.lvl{1}=NaN;
plotwhat.dim{1}=[1,1];
plotfield(dns, field, plotit,plotwhat)


