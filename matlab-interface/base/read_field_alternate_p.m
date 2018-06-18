function field = read_field_alternate_p(filename,dns,field,dc)

%
% This function reads the velocity field out of the field data produced by
% pipe.cpl
%

% Declare array
field.p=cell(dns.ny+1,1); 
for IY=0:dns.ny; iy=IY+1;
    field.p{iy}=complex(zeros(2*field.nzN(iy)+1,dns.nx+1),0);
end

% Open file
f = fopen(filename);

% Read data from file
for IY=1:dns.ny
    iy=IY+1;
    fseek(f,field.startpos_p(IY),'bof');
    dum = reshape(fread(f,2*(dns.nx+1)*(field.nzN(iy)*2+1),'double'),[2,field.nzN(iy)*2+1,dns.nx+1]); 
    field.p{iy}(:,:) = -0.5*1j*complex(dum(1,:,:),dum(2,:,:));
end

% Use continuity condition at the centerline of the pipe to compute the
% pressure at index 1
for IX=0:dns.nx
    ix=IX+1;
    for IZ=-field.nzN(1):field.nzN(1); iz=IZ+field.nzN(1)+1;
            field.p{1}(iz,ix)=complex(0,0);
            for i=2:3
              jz=IZ+field.nzN(i)+1;
              field.p{1}(iz,ix) = field.p{1}(iz,ix) - dc(abs(IZ)+1,2,i).*field.p{i}(jz,ix);
            end
    end
end

fclose(f);
