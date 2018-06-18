% This program reads the files produced by pipe.cpl

% MATLAB has lot of built-in function which you can use and which I also
% used in this program. By typing
%
% help <command> 
%
% in the Command window below you can read what a comamnd does and
% its proper syntax.

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


% Start measure of time
tic

%% Clear Command Window and Workspace
%  -------------------------------------------------------------------
clc
clear variables

%% Add subforlder ./base/ into MATLAB search paths
%  -------------------------------------------------------------------
addpath('./base/');

%% Read field header, set up coordinates and read initial field
%  -------------------------------------------------------------------
[dns,field] = read_header('vfield.dat');
setup_derivatives
[field] = read_field_alternate('vfield.dat',dns,field,dc);

%% Compute velocity gradient and pressure
[field] = poisson(dns, field, derivatives, dc);

%% Fourier transorm the pipe data to space
%  -------------------------------------------------------------------
field.Vr=cell(dns.ny+1,1);  % declare cell array for velocity in physical space
for IY=0:dns.ny
    iy=IY+1; field.Vr{iy}=complex(zeros(3,field.nzd(iy),2*dns.nxd),0);  % for every iy in the cell array declare matrix of velocity (3 velocity components, azimultal mode, axial mode)
    field.Vr{iy} = plane_ift(field.V{iy},field.Vr{iy});
end

% % Plot one (r,x) slice of vertical velocity
% for iy=1:dns.ny+1; for ix=1:2*dns.nxd
%         Vz(iy,ix)=field.Vr{iy}(1,10,ix); 
%         Vr(iy,ix)=field.Vr{iy}(2,10,ix);
%         Vtheta(iy,ix)=field.Vr{iy}(3,10,ix);
% end; end
% contourf(linspace(0,2*pi/dns.alfa0,2*dns.nxd), field.y, real(Vtheta))
% hold on; quiver(linspace(0,2*pi/dns.alfa0,2*dns.nxd),field.y,real(Vz),real(Vr),'k')
% xlabel('z')
% ylabel('r')
% axis equal


toc