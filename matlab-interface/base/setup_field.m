function  [dns,field] = setup_field(dns)

% Compute extended number of modes
dns.nxd=3*dns.nx/2; while ~fftfit(dns.nxd); dns.nxd=dns.nxd+1; end   
dns.nzd=3*dns.nz;   while ~fftfit(dns.nzd); dns.nzd=dns.nzd+1; end

% y-coordinate
field.y=(dns.ymin +(dns.ymax-dns.ymin)*tanh(dns.htcoef*(0:dns.ny)/dns.ny)/tanh(dns.htcoef))';

% index (radius) iy0 at which the mode iz starts to appear
field.iy0=zeros(2*dns.nz+1,1,'int32');
for IZ=0:dns.nz; iz=dns.nz+1+IZ;
    field.iy0(iz)=0; 
    while field.y(field.iy0(iz)+3+1)*(dns.nz*(field.y(field.iy0(iz)+3+1)>=0.2)+100*(field.y(field.iy0(iz)+3+1)<0.2)-0.5) < (IZ-1)*dns.ymax
        field.iy0(iz)=field.iy0(iz)+1; 
    end
end
for m=1:dns.nz; field.iy0(dns.nz+1-m)=field.iy0(dns.nz+1+m); end

% Define positions in the file from which to read correct velocities
field.startpos=zeros(dns.ny,1,'uint64');
field.startpos(1)=1024+numel(field.y)*8+(dns.nz+1)*4;
iz=1;
for iy=1:dns.ny-1
    while iz<dns.nz && iy>=field.iy0(dns.nz+1+iz+1); iz=iz+1; end
    field.startpos(iy+1)=field.startpos(iy)+3*8*2*(dns.nx+1)*(2*iz+1);
end

% Define positions in the file from which to read correct pressure
field.startpos_p=zeros(dns.ny,1,'uint64');
field.startpos_p(1)=1024+numel(field.y)*8+(dns.nz+1)*4;
iz=1;
for iy=1:dns.ny-1
    while iz<dns.nz && iy>=field.iy0(dns.nz+1+iz+1); iz=iz+1; end
    field.startpos_p(iy+1)=field.startpos_p(iy)+8*2*(dns.nx+1)*(2*iz+1);
end

% Count the variable number of Fourier modes in nz that we have
% and declare variables
field.V=cell(dns.ny+1,1); field.nzN=zeros(dns.ny+1,1);
for IY=0:dns.ny
    iy=IY+1; IZ=0;
    while IZ<=dns.nz && IY>=field.iy0(IZ+dns.nz+1); IZ=IZ+1; end
    field.nzN(iy)=IZ-1;
    field.V{iy}=complex(zeros(3,2*field.nzN(iy)+1,dns.nx+1),0);
end

% Define nzd(iy)
field.nzd = zeros(dns.ny+1,1);
for IY=0:dns.ny
    iy=IY+1;
    field.nzd(iy)=dns.nz; 
    while field.nzd(iy)~=0 && IY<field.iy0(dns.nz+1+field.nzd(iy)); field.nzd(iy)=field.nzd(iy)-1; end
    field.nzd(iy)=3*field.nzd(iy);
    while ~fftfit(field.nzd(iy)); field.nzd(iy)=field.nzd(iy)+1; end
end

% Constant number of points in radial direction
%for IY=0:dns.ny
%    iy=IY+1; field.nzd(iy)=field.nzd(end);
%end
end