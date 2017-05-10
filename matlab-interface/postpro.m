%
% This program reads, hopefully in an intelligent and smart way, the fields
% produced by pipe.cpl! 
%

% super shlau pipe :)

tic

clc
clear variables
addpath('./base/');

[dns,field] = read_header('vfield.dat');
setup_derivatives
[field] = read_field_alternate('vfield.dat',dns,field,dc);

% Fourier transorm the pipe data to space
field.Vr=cell(dns.ny+1,1); 
for IY=0:dns.ny
    iy=IY+1; field.Vr{iy}=complex(zeros(3,field.nzd(iy),2*dns.nxd),0);
    field.Vr{iy}(:,end+(-field.nzN(iy):-1)+1,(0:dns.nx)+1) = field.V{iy}(:,(-field.nzN(iy):-1)+field.nzN(iy)+1,:);
    field.Vr{iy}(:,1+(0:field.nzN(iy)),(0:dns.nx)+1) = field.V{iy}(:,(0:field.nzN(iy))+field.nzN(iy)+1,:);
    for iV=1:3
        for ix=1:2*dns.nxd;     field.Vr{iy}(iV,:,ix) = ifft(squeeze(field.Vr{iy}(iV,:,ix))); end
        field.Vr{iy}(iV,:,end+(-dns.nx:-1)+1) = conj(field.Vr{iy}(iV,:,(dns.nx:-1:1)+1));
        for iz=1:field.nzd(iy); field.Vr{iy}(iV,iz,:) = ifft(squeeze(field.Vr{iy}(iV,iz,:)))*(2*field.nzd(iy)*dns.nxd); end
    end
end

% Plot one (r,x) slice of vertical velocity
for iy=1:dns.ny+1; for ix=1:2*dns.nxd;  V(iy,ix)=mean(field.Vr{iy}(3,10,ix)); end; end
contourf(linspace(0,2*pi/dns.alfa0,2*dns.nxd), field.y, real(V))
xlabel('x')
ylabel('r')

toc