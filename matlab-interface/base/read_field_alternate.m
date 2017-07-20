function field = read_field_alternate(filename,dns,field,dc)

%
% This function reads the velocity field out of the field data produced by
% pipe.cpl
%

% Open file
f = fopen(filename);

% Read data from file
for IY=1:dns.ny-1
    iy=IY+1;
    fseek(f,field.startpos(IY),'bof');
    
    dum = reshape(fread(f,6*(dns.nx+1)*(field.nzN(iy)*2+1),'double'),[2,3,field.nzN(iy)*2+1,dns.nx+1]);
    field.V{iy}(1:3,:,:) = complex(dum(1,:,:,:),dum(2,:,:,:));
%     
%     for IX=0:dns.nx
%         ix=IX+1;
%         for IZ=-dns.nz:dns.nz
%             iz=IZ+field.nzN(iy)+1; 
%             if IY>=field.iy0(IZ+dns.nz+1)
%               dum = reshape(fread(f,6,'double'),[2,3]);
%               field.V{iy}(1:3,iz,ix) = complex(dum(1,:),dum(2,:));
%             end
%         end
%     end
%     
    
    
end

% Use continuity condition at the centerline of the pipe to compute the
% velocity at index 1
for IX=0:dns.nx
    ix=IX+1;
    for IZ=-field.nzN(1):field.nzN(1)
            field.V{1}(1:3,IZ+field.nzN(1)+1,ix)=complex(0,0);
            for i=2:3
              field.V{1}(1,IZ+field.nzN(1)+1,ix) = field.V{1}(1,IZ+field.nzN(1)+1,ix) - dc(abs(IZ)+1,2,i).*field.V{i}(1,IZ+field.nzN(i)+1,ix);
              field.V{1}(2,IZ+field.nzN(1)+1,ix) = field.V{1}(2,IZ+field.nzN(1)+1,ix) - dc(abs(IZ)+1,3,i).*field.V{i}(2,IZ+field.nzN(i)+1,ix);
              field.V{1}(3,IZ+field.nzN(1)+1,ix) = field.V{1}(3,IZ+field.nzN(1)+1,ix) - dc(abs(IZ)+1,4,i).*field.V{i}(3,IZ+field.nzN(i)+1,ix);
            end
    end
end

fclose(f);