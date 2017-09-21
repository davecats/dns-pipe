% This program reads the files produced by pipe.cpl
% Throughout this program I use following convenctions for indexing
% variables. Indexing variables written in CAPITAL LETTERS, for instance
%
% IZ (index correspoding to a wavenumber in azimutal direction)
%
% are indices carying a logical meaning, for  instance in this case IZ 
% varies between -nz and nz. However, MATLAB Matrices always start from
% index 1. Therefore I define indices written in small letters, which are
% the corresponding MATLAB indices. In the example then
%
% iz 
%
% would go from 1 to 2*nz+1. The transoformation from one to the other
% would then be: iz = IZ + nz + 1
%%
% This program calculates the statistics:
% mean velocity, mean vorticity, reynolds stresses, tke, mke, ke, y+, u_tau etc. etc.

%% Clear Command Window and Workspace
%  -------------------------------------------------------------------
clc
clear variables
restoredefaultpath

% Start measure of time
tic

%% Options
% Calculates the mean vorticity, enstrophy and the vorticity variance:
% Needs some RAM
vorticityfull=0;

% Give the maximum angular speed
omegamax=1;
stats.omegamax=omegamax;

%% Add subforlder ./base/ into MATLAB search paths
%  -------------------------------------------------------------------
addpath('./../base/');
addpath('./utilities/');
%% Provide path of the data

%Foldername
% filepath_short='Re28000/nx-256-ny-500-nz-512-alfa0-4.0/';
% filepath_short='Re3000/nx-32-ny-128-nz-128-alfa0-alpha0-8.0-omega-0.5/';
filepath_short='Re3000/nx-32-ny-128-nz-128-alfa0-alpha0-2.0/';
% filepath_short='Re3000/nx-32-ny-128-nz-128-alfa0-alpha0-8.0-omega-2/';
% filepath_short='Re10000/nx-256-ny-500-nz-512-alfa0-4.0/'

% Full path
filepath=['../../../DATA/',filepath_short]
%Get all names of files in the destination
filenamelist=dir([filepath,'/vfield*.dat']);
% sort them (time)
filenamelist= sort_nat({filenamelist.name}.');

%% Read field header, set up coordinates and read initial field
%  -------------------------------------------------------------------
filename = strcat(filepath,filenamelist{150});
[dns,field] = read_header(filename);
% read out the radius of the pipe
fieldymax=field.y(end);
stats.fieldymax=fieldymax;
% normalize the radius
field.y=field.y/fieldymax;

% setup the derivatives
setup_derivatives
disp(['setup_derivatives took ',num2str(toc)]);

%% Set temporal limits of the upcoming evaluation
iF_start=1;
iF_end=numel(filenamelist);
% spacing between 2 evaluated timesteps
iF_space=2;
size_t=numel(iF_start:iF_space:iF_end);

%% Prelocate quantities
stats.t=NaN(size_t,1); %time
stats.Vm=NaN(3,dns.ny+1,size_t); %mean velocity
stats.var=NaN(3,3,dns.ny+1,size_t); %variance/Reynolds stress tensor (3 velocity components, 3 spatial derivatives)

stats.tke=NaN(dns.ny+1,size_t); %TKE(r,t)
stats.tke_int=NaN(size_t,1); %TKE(t)
stats.mke=NaN(dns.ny+1,size_t); %mean kinetic energy (r,t)
stats.mke_int=NaN(size_t,1); %mean kinetic energy (t)
stats.ke=NaN(dns.ny+1,size_t); %total kinetic energy (r,t)
stats.ke_int=NaN(size_t,1); %total kinetic energy (t)

stats.ypluspos=NaN(3,3,size_t); %position of the peak of  the maximum of the Reynolds stress tensor in wall units
stats.u_tau=NaN(size_t,1); %friction/shear velocity 
stats.tau_w=NaN(size_t,1); %wall shear stress
stats.yplus=NaN(dns.ny+1,size_t); %normalization of radial direction into wall units
stats.d99=NaN(size_t,1); %boundary layer thickness d99
% stats.d1=NaN(size_t); 
% stats.d2=NaN(size_t);
% stats.dfitpot=NaN(size_t);
% stats.dfitstraight=NaN(size_t);

stats.Re_TKE=NaN(size_t,1); %Reynolds number based on max turbulent kinetic energy
stats.Re_TKEmean=NaN(size_t,1); %Reynolds number based on mean turbulent kinetic energy
stats.Re_TC=NaN(size_t,1); %Reynolds number as calculated for TC (with d99 as radius of the inner cylinder)
stats.Re_tau=NaN(size_t,1); %Reynolds number based on on the wall shear stress

stats.omegam=NaN(3,dns.ny+1,size_t); %mean vorticity
stats.omega_int=NaN(3,size_t); %vorticity


%% Calculates Stats 
%-----------------------------------------------------------------------
counter=0;
colorcount=0;
for iF=iF_start:iF_space:iF_end
    counter=counter+1;
    %load data into 'field'
    filename = strcat(filepath,filenamelist{iF});
    disp(['read field  ',filenamelist{iF}]);
    [field] = read_field_alternate(filename,dns,field,dc); %reads V of the particular timestep
    %Normalize velocity
    for iy=1:dns.ny+1;
        field.V{iy}=field.V{iy}/(omegamax*fieldymax);
    end
    %get time
    f = fopen(filename);
    header=fread(f,1024,'*char')';
    i=find(header=='=');
    fclose(f);
    t(iF) = str2double(header(i(9)+1 : end));
    stats.t(counter)=t(iF);
    
    %Calculate statistics
    field=velocity_mean(dns, field, derivatives, dc); %mean_velocities
    field=velocity_variance(dns, field, derivatives, dc); %variance of velocities
    field=bl_thinkness(dns, field, derivatives, dc); %bl thickness
    if vorticityfull
        field=velocity_gradient(dns, field, derivatives, dc);
        field=vorticity(dns, field, derivatives, dc);
    end
    stats.Vm(:,:,counter)=field.Vm;
    stats.var(:,:,:,counter)=field.var;
    % Turbulent kinetic energy
    stats.tke(:,counter)=0.5*(field.var(1,1,1,:)+field.var(2,2,1,:)+field.var(3,3,1,:));
    stats.tke_int(counter)=sum(0.25*(field.y(1:end-1)+field.y(2:end)).*(stats.tke(1:end-1,counter)+stats.tke(2:end,counter)) ...
        .*(field.y(2:end)-field.y(1:end-1)));
    % Mean kinetic energy
    stats.mke(:,counter)=0.5*(field.Vm(1,:).^2+field.Vm(2,:).^2+field.Vm(3,:).^2);
    stats.mke_int(counter)=sum(0.25*(field.y(1:end-1)+field.y(2:end)).*(stats.mke(1:end-1,counter)+stats.mke(2:end,counter)) ...
        .*(field.y(2:end)-field.y(1:end-1)));
    % Total kinetic engergy
    for IY=0:dns.ny
        iy=IY+1;
        ke_y=complex(zeros(2*field.nzN(iy)+1,dns.nx+1),2*dns.nxd);
        tmpVr  =complex(zeros(1,3,field.nzd(iy),2*dns.nxd),0);
        tmpVr  =plane_ift(reshape(field.V{iy},[1,3,2*field.nzN(iy)+1,dns.nx+1]),tmpVr);
        ke_y(:,:)=plane_fft( ...
            0.5*(tmpVr(1,1,:,:).^2+tmpVr(1,2,:,:).^2+tmpVr(1,3,:,:).^2), ...
            reshape(ke_y(:,:),[1,1,2*field.nzN(iy)+1,dns.nx+1]));
        ke(iy)=ke_y(field.nzN(iy)+1,1);
    end
    stats.ke(:,counter)=ke;
    stats.ke_int(counter)=sum(0.25*(field.y(1:end-1)+field.y(2:end)).*(stats.ke(1:end-1,counter)+stats.ke(2:end,counter)) ...
        .*(field.y(2:end)-field.y(1:end-1)));
    
    % Law of the wall
    field.dVm=fieldymax*derivatives.d0{field.nzN(end)+1}\(derivatives.d1{field.nzN(end)+1}*(omegamax*fieldymax*field.Vm(3,:)'));
    nu=omegamax*fieldymax^2/dns.Re;
    tau_w=nu*field.dVm(end);
    u_tau=sqrt(abs(tau_w));
    yplus=(fieldymax-field.y*fieldymax)*u_tau/nu;
    stats.tau_w(counter)=tau_w;
    stats.u_tau(counter)=u_tau;
    stats.yplus(:,counter)=yplus;
    %Position of the maxima of the Reynolds stresses
    for iV1=1:3; for iV2=iV1:3
            maxpos=find(field.var(iV1,iV2,1,:)==max(field.var(iV1,iV2,1,:)));
            F = griddedInterpolant(-yplus,squeeze(field.var(iV1,iV2,1,:)),'spline','none');
            xmaxs = -arrayfun(@(xx)fminsearch(@(x)-F(x),xx),-yplus(maxpos));
            stats.ypluspos(iV1,iV2,counter)=max(xmaxs);
     end;end
    
    stats.d99(counter)=field.d99;
    stats.d99_pos(counter)=field.d99_pos;
    %Re_tau=u_tau*d99/nu (POPE S.270)
    stats.Re_tau(counter)=u_tau*1*dns.Re*(1-stats.d99(counter));
    %Maximum of TKE *R/nu
    stats.Re_TKE(counter)=max(squeeze(sqrt(field.var(1,1,1,:)+field.var(2,2,1,:)+field.var(3,3,1,:))))*1*dns.Re;
    % d99*(1-d99)*Omega/nu
    stats.Re_TC(counter)=1*field.d99*(1-field.d99)*dns.Re; %omega*Ri*d/nu
    % MeanTKE *R/nu
    stats.Re_TKEmean(counter)=sqrt(stats.tke_int(counter))*dns.Re;

    % Mean vorticity in axial direction 
    stats.omegam(1,:,counter)=squeeze(field.Vm(3,:))'./field.y+squeeze(field.dVm(:));
    for iV=1:3
        stats.omega_int(iV,counter)=sum(0.25*(field.y(1:end-1)+field.y(2:end)).*(stats.omegam(iV,1:end-1,counter)+stats.omegam(iV,2:end,counter))' ...
            .*(field.y(2:end)-field.y(1:end-1)));
    end
    
    
    
end
stats.y=field.y;
mkdir(['../../../Results/',filepath_short]);
save(['../../../Results/',filepath_short,'Statistics.mat'],'stats');


