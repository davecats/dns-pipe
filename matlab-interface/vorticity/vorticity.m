function   [field] = vorticity(dns, field, derivatives, dc)
% This function calculates the vorticty for a given field
if ~isfield(field,'dU')
   disp('execute velocity_gradient.m first') 
end
%% Alocate vorticity
field.omega=cell(dns.ny+1,1); 
%% Calculate vortictiy for IY=2:end
for iy=1+(1:dns.ny)
        field.omega{iy}(1,:,:)=squeeze(field.V{iy}(3,:,:))/field.y(iy)+squeeze(field.dU{iy}(3,2,:,:))-squeeze(field.dU{iy}(2,3,:,:))/field.y(iy);
        field.omega{iy}(2,:,:)=squeeze(field.dU{iy}(1,3,:,:)/field.y(iy)-field.dU{iy}(3,1,:,:));
        field.omega{iy}(3,:,:)=squeeze(field.dU{iy}(2,1,:,:))-squeeze(field.dU{iy}(1,2,:,:));
end

%% Use continuity condition at the centerline of the pipe to compute the vorticity at index 1
field.omega{1}=complex(zeros(3,2*field.nzN(1)+1,dns.nx+1),0);
for IX=0:dns.nx
    ix=IX+1;
    for IZ=-field.nzN(1):field.nzN(1); iz=IZ+field.nzN(1)+1;
            field.omega{1}(1:3,iz,ix)=complex(0,0);
            for i=2:3
              jz=IZ+field.nzN(i)+1;
              field.omega{1}(1,iz,ix) = field.omega{1}(1,iz,ix) - dc(abs(IZ)+1,2,i).*field.omega{i}(1,jz,ix);
              field.omega{1}(2,iz,ix) = field.omega{1}(2,iz,ix) - dc(abs(IZ)+1,3,i).*field.omega{i}(2,jz,ix);
              field.omega{1}(3,iz,ix) = field.omega{1}(3,iz,ix) - dc(abs(IZ)+1,4,i).*field.omega{i}(3,jz,ix);
            end
    end
end