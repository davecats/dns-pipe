function   [field] = velocity_gradient(dns, field, derivatives, dc)

%% Allocate velocity gradient
field.dU=cell(dns.ny+1,1);
for IY=0:dns.ny; iy=IY+1; field.dU{iy}=complex(zeros(3,3,2*field.nzN(iy)+1,dns.nx+1),0); end
%% Compute velocity gradient
for IY=0:dns.ny
    iy=IY+1;
    kx = repmat(reshape((0:dns.nx)*dns.alfa0          ,[1,1,dns.nx+1]         ),3,2*field.nzN(iy)+1,1);
    kz = repmat(reshape((-field.nzN(iy):field.nzN(iy)),[1,2*field.nzN(iy)+1,1]),3,1                ,dns.nx+1);
    field.dU{iy}(:,1,:,:) = 1j.*kx.*field.V{iy}(:,:,:);
    field.dU{iy}(:,3,:,:) = 1j.*kz.*field.V{iy}(:,:,:);
end
for IZ=-dns.nz:dns.nz; iz=dns.nz+1+IZ;
    iy0=field.iy0(iz);
    tmpU = complex(zeros(dns.ny-iy0+1,3*(dns.nx+1)));
        % matlab does not allow : along cells of arrays, so copy...
        for IY=iy0:dns.ny; iy=IY+1; jy=IY-iy0+1; jz=field.nzN(iy)+1+IZ; 
            tmpU(jy,:) = reshape(field.V{iy}(:,jz,:),[1,3*(dns.nx+1)]);
        end
        % copute compact derivatives from iy0+1.ny
        tmpU = derivatives.d0{iz}\(derivatives.d1{iz}*tmpU);
        % recopy back to the velocity gradient array
        for IY=iy0:dns.ny; iy=IY+1; jy=IY-iy0+1; jz=field.nzN(iy)+1+IZ; 
            field.dU{iy}(:,2,jz,:) = reshape(tmpU(jy,:),[3,1,1,dns.nx+1]);
        end
end
