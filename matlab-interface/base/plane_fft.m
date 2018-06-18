function dataOUT=plane_fft(dataIN,dataOUT)

sout = size(dataOUT);  
ndims = length(sout); idx=repmat({':'}, 1, ndims); dims=[ndims ndims-1];
nz=sout(dims(2));   nx=sout(dims(1));   nz=(nz-1)/2; nx=nx-1;
sin= size(dataIN); nzd=sin(dims(2));   nxd=sin(dims(1)); nxd=nxd/2;

% 1D transform in the axial direction
dataIN = fft(dataIN,[],dims(1))/(2*nzd*nxd);

% Backward fourtier transform
idxIN = idx;                  idxOUT=idx; 
idxIN{dims(1)}=(0:nx)+1;      idxOUT{dims(1)}=(0:nx)+1;
dataIN(idxIN{:}) = fft(dataIN(idxOUT{:}),[],dims(2));

% Swap the two half of the vector from (-n..n) --> (0..n, -n..-1)
idxIN = idx;                               idxOUT=idx; 
idxIN{dims(2)}=[nzd+(-nz:-1)+1 1+(0:nz)];  idxOUT{dims(2)}=[(-nz:-1)+nz+1 (0:nz)+nz+1];
idxIN{dims(1)}=(0:nx)+1;                   idxOUT{dims(1)}=(0:nx)+1;
dataOUT(idxOUT{:}) = dataIN(idxIN{:});

end