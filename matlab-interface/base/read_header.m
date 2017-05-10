function  [dns,field] = read_field(filename)

f = fopen(filename);
header=fread(f,1024,'*char')';
i=find(header=='=');
dns.nx=str2double(header(i(1)+1 : i(2)-3));
dns.ny=str2double(header(i(2)+1 : i(3)-3));
dns.nz=str2double(header(i(3)+1 : i(4)-7));
dns.alfa0=str2double(header(i(4)+1 : i(5)-7));
dns.htcoef=str2double(header(i(5)+1 : i(6)-3));
if numel(strfind(header,'meanpx'))>0
    dns.Re=str2double(header(i(6)+1 : i(7)-7));
    dns.meanpx=str2double(header(i(7)+1 : i(8)-7));
else
    dns.Re=str2double(header(i(6)+1 : i(7)-10));
    dns.meanflowx=str2double(header(i(7)+1 : i(8)-7));
end
dns.u_conv=str2double(header(i(8)+1 : i(9)-5));
dns.time=str2double(header(i(9)+1 : end));

dns.nxd=3*dns.nx/2; while ~fftfit(dns.nxd); dns.nxd=dns.nxd+1; end
dns.nzd=3*dns.nz;   while ~fftfit(dns.nzd); dns.nzd=dns.nzd+1; end

fseek(f,1024,'bof');

field.y=fread(f,dns.ny+1,'double');
field.iy0=zeros(2*dns.nz+1,1,'int32'); field.iy0(dns.nz+1:end)=fread(f,dns.nz+1,'int32'); 
for m=1:dns.nz; field.iy0(dns.nz+1-m)=field.iy0(dns.nz+1+m); end

field.startpos=zeros(dns.ny,1,'uint64');
field.startpos(1)=1024+numel(field.y)*8+(dns.nz+1)*4;
iz=1;
for iy=1:dns.ny-1
    while iz<dns.nz && iy>=field.iy0(dns.nz+1+iz+1); iz=iz+1; end
    field.startpos(iy+1)=field.startpos(iy)+3*8*2*(dns.nx+1)*(2*iz+1);
end

% Count the number of Fourier modes in nz that we have
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

fclose(f);