function field = read_field(filename,dns,field,dc)

f = fopen(filename);

% Declare variables
field.V=cell(2*dns.nz+1,dns.nx+1);
for IX=0:dns.nx
    ix=IX+1;
    for IZ=-dns.nz:dns.nz
        iz=IZ+dns.nz+1;
        field.V{iz,ix} = complex(zeros(3,dns.ny-max([field.iy0(iz) -1])+1));
    end
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

% Read data from file
for IY=1:dns.ny-1
    iy=IY+1;
    fseek(f,field.startpos(IY),'bof');
    for IX=0:dns.nx
        ix=IX+1;
        for IZ=-dns.nz:dns.nz
            iz=IZ+dns.nz+1; 
            if IY>=field.iy0(iz)
              %disp(field.startpos(IY)); %disp(IX); disp(IZ); disp(IY); return
              dum = reshape(fread(f,6,'double'),[2,3]);
              field.V{iz,ix}(1:3,iy) = complex(dum(1,:),dum(2,:));
            end
        end
    end
end

% 
for IX=0:dns.nx
    ix=IX+1;
    for IZ=-dns.nz:dns.nz
        iz=dns.nz+IZ+1;
        if field.iy0(iz)>0
            field.V{iz,ix}(1,1) = - sum( squeeze(dc(abs(IZ)+1,1,2:3)).*squeeze(field.V{iz,ix}(1,2:3))' );
            field.V{iz,ix}(2,1) = - sum( squeeze(dc(abs(IZ)+1,2,2:3)).*squeeze(field.V{iz,ix}(2,2:3))' );
            field.V{iz,ix}(3,1) = - sum( squeeze(dc(abs(IZ)+1,3,2:3)).*squeeze(field.V{iz,ix}(3,2:3))' );
        end
    end
end

fclose(f);