function dataOUT=plane_ift(dataIN,dataOUT)

dataIN=squeeze(dataIN);    dataOUT=squeeze(dataOUT);
[nz,nx] = size(dataIN);    nz=(nz-1)/2; nx=nx-1;
[nzd,nxd] = size(dataOUT); nxd=nxd/2;

% Swap the two half of the vector from (-n..n) --> (0..n, -n..-1)
dataOUT(end+(-nz:-1)+1,(0:nx)+1) = dataIN((-nz:-1)+nz+1,:);
dataOUT(1+(0:nz),(0:nx)+1)       = dataIN((0:nz)+nz+1,:);

% Backward fourtier transform
for ix=1:nx+1; dataOUT(:,ix) = ifft(dataOUT(:,ix)); end  
% reconstruct missing coefficients in x direction which were not
% saved due to simmetry
dataOUT(:,end+(-nx:-1)+1) = conj(dataOUT(:,(nx:-1:1)+1));
% 1D transform in the axial direction
for iz=1:nzd; dataOUT(iz,:) = ifft(dataOUT(iz,:))*(2*nzd*nxd); end


end