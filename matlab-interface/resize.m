%
% This script resizes a velocity field
%

%% Clear Command Window and Workspace
%  -------------------------------------------------------------------
clc
clear variables

%% INPUTS
%  -----------------------------------
old.filename=('vfield341.dat');
new.filename=('resize.dat');
new.dns.nx=400;   new.dns.ny=500; new.dns.nz=400;
new.dns.alfa0=2;  new.dns.htcoef=1.2;
new.dns.meanpx=0; new.dns.Re=28000;
new.dns.u_conv=0; new.dns.time=27.9;
new.dns.ymin=0;   new.dns.ymax=1;
%  -----------------------------------

% Start measure of time
tic

%% Add subforlder ./base/ into MATLAB search paths
%  -------------------------------------------------------------------
addpath('./base/');

%% Read field header, set up coordinates and read field to be resized
%  -------------------------------------------------------------------
[old.dns,old.field] = read_header(old.filename);
[old.der,old.dc]=setup_derivatives(old.dns,old.field,3);
[old.field] = read_field_alternate(old.filename,old.dns,old.field,old.dc);

%% Set up coordinates for the new field
%  -------------------------------------------------------------------
[new.dns,new.field] = setup_field(new.dns);

%1d spektrum x
old.psdx =zeros(old.dns.nx+1,old.dns.ny+1);
old.psdz =zeros(2*old.dns.nz+1,old.dns.ny+1);
for iy=1:old.dns.ny+1
    old.psdx(:,iy)=2*sum(abs(old.field.V{iy}(1,:,:)),2); old.psdx(:,iy)=old.psdx(:,iy)/max(old.psdx(:,iy));
    n=(numel(old.field.V{iy}(1,:,1))-1)/2; N=2*old.dns.nz+1;
    old.psdz((N-1)/2+1+(-n:n),iy)=sum(squeeze(abs(old.field.V{iy}(1,:,[1:end 2:end]))),2); old.psdz(:,iy)=old.psdz(:,iy)/max(old.psdz(:,iy));
end
figure(); contour(old.field.y,(-old.dns.nz:old.dns.nz),log(old.psdz),[-10 -7 -5 -3],'ShowText','on'); hold on; plot(old.field.y,old.field.nzN+1,'k--','Linewidth',2); plot(old.field.y,new.field.nzN+1,'r--','Linewidth',2); xlabel('r/R'); ylabel('nz(r)'); ylim([0 old.dns.nz]); 
figure(); contour(old.field.y,(0:old.dns.nx),log(old.psdx),[-10 -7 -5 -3],'ShowText','on'); hold on;  plot(old.field.y,new.dns.nx*ones(size(old.field.y)),'k','Linewidth',2); ylim([0 old.dns.nx]); xlabel('r/R'); ylabel('nx(r)');


%% Resize
%  -------------------------------------------------------------------
IXmin=min([old.dns.nx new.dns.nx]); IZmin=min([old.dns.nz new.dns.nz]);
for IX=0:IXmin; ix=IX+1;
    disp([num2str(IX) ' of ' num2str(IXmin)])
    for IZ=-IZmin:IZmin; oldiz=IZ+old.dns.nz+1; newiz=IZ+new.dns.nz+1;
        old.Vpencil=complex(zeros(3,old.dns.ny-old.field.iy0(oldiz)+1),0);
        old.ypencil=zeros(1,old.dns.ny-old.field.iy0(oldiz)+1);
        new.Vpencil=complex(zeros(3,new.dns.ny-new.field.iy0(newiz)+1),0);
        new.ypencil=zeros(1,new.dns.ny-new.field.iy0(newiz)+1);
        % Load iy-pencil from old data and do cell->vector
        for IY=old.field.iy0(oldiz):old.dns.ny; iy=IY+1; jz=old.field.nzN(iy)+1+IZ; jy=iy-old.field.iy0(oldiz);
            old.Vpencil(:,jy)=old.field.V{iy}(:,jz,ix);
        end
        old.ypencil=old.field.y(old.field.iy0(oldiz)+1:end);
        new.ypencil=new.field.y(new.field.iy0(newiz)+1:end);
        % Interpolate iy-pencil to new data
        for iV=1:3
            new.Vpencil(iV,:)=interp1(squeeze(old.ypencil),squeeze(old.Vpencil(iV,:)),squeeze(new.ypencil),'pchip','extrap');
        end
        % Save iy-pencil into new data and do vector->cell
        for IY=new.field.iy0(newiz):new.dns.ny; iy=IY+1; jz=new.field.nzN(iy)+1+IZ; jy=iy-new.field.iy0(newiz);
            new.field.V{iy}(:,jz,ix)=new.Vpencil(:,jy);
        end
    end
end

%% Save data to file
%  -------------------------------------------------------------------
% Open file
f=fopen(new.filename,'w+');
header = sprintf('nx=%d\tny=%d\tnz=%d\nalpha0=%g\thtcoef=%g\tRe=%g\nmeanpx=%g\nu_conv=%g\ntime=%g\n', new.dns.nx,new.dns.ny,new.dns.nz,new.dns.alfa0,new.dns.htcoef,new.dns.Re,new.dns.meanpx,new.dns.u_conv,new.dns.time);
header = [header, repmat(' ',[1, 1024-numel(header)])];
fprintf(f,header);
% Go past header in the file
fseek(f,1024,'bof');
% Write y-coordinate & iy0(0..nz)
fwrite(f,new.field.y,'double'); fwrite(f,new.field.iy0(new.dns.nz+1:end),'int32');
% Write velocity field
for IY=1:new.dns.ny; iy=IY+1;
    fseek(f,new.field.startpos(IY),'bof');
    buf=zeros(2,3,new.field.nzN(iy)*2+1,new.dns.nx+1);
    buf(1,:,:,:)=real(new.field.V{iy}(1:3,:,:)); buf(2,:,:,:)=imag(new.field.V{iy}(1:3,:,:));
    fwrite(f,buf,'double');
end
fclose(f);

% Stop measure of time
toc


