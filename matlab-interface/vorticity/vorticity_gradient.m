function [field] = vorticity_gradient(dns, field, derivatives, dc)

%This function computes the gradient of vorticity (jacobi) and the second spatial derivatives
%note that for the second radial derivatives drd=d/dr(r(d/dr(omega))) is used

%% Allocate vorticity jacobian and hessian
field.domega=cell(dns.ny+1,1); 
field.domega2=cell(dns.ny+1,1); % second spatial derivative:
     %order: omega_zz,drd (omega), omega_phiphi
%% Compute gradient of vorticity (jacobi) and the second spatial derivatives (hessian)
% azimuthal and spatial direction
for IY=0:dns.ny
    iy=IY+1;
     % Allocate vorticity gradient(iy)
     field.domega{iy}=complex(zeros(3,3,2*field.nzN(iy)+1,dns.nx+1),0);
     field.domega2{iy}=complex(zeros(3,3,2*field.nzN(iy)+1,dns.nx+1),0); 
        
    kx = repmat(reshape((0:dns.nx)*dns.alfa0          ,[1,1,dns.nx+1]         ),3,2*field.nzN(iy)+1,1);
    kz = repmat(reshape((-field.nzN(iy):field.nzN(iy)),[1,2*field.nzN(iy)+1,1]),3,1                ,dns.nx+1);
    field.domega{iy}(:,1,:,:) = 1j.*kx.*field.omega{iy}(:,:,:); % Derivation into axial direction
    field.domega{iy}(:,3,:,:) = 1j.*kz.*field.omega{iy}(:,:,:); % Derivation into azimuthal direction
    field.domega2{iy}(:,1,:,:) = -1.*kx.*kx.*field.omega{iy}(:,:,:); % 2nd Derivation into axial direction
    field.domega2{iy}(:,3,:,:) = -1.*kz.*kz.*field.omega{iy}(:,:,:); % 2nd Derivation into azimuthal direction

end
% radial direction
for IZ=-dns.nz:dns.nz; iz=dns.nz+1+IZ;
    iy0=field.iy0(IZ+dns.nz+1);
    tmpOm = complex(zeros(dns.ny-iy0+1,3*(dns.nx+1))); % for d1
    tmpOm2 = complex(zeros(dns.ny-iy0+1,3*(dns.nx+1))); %for d2
        % matlab does not allow : along cells of arrays, so copy...
        for IY=iy0:dns.ny; iy=IY+1; jy=IY-iy0+1; jz=field.nzN(iy)+1+IZ; 
            tmpOm(jy,:) = reshape(field.omega{iy}(:,jz,:),[1,3*(dns.nx+1)]);
            tmpOm2(jy,:) = reshape(field.omega{iy}(:,jz,:),[1,3*(dns.nx+1)]);
        end
        % compute compact derivatives from iy0+1.ny
        tmpOm = derivatives.d0{iz}\(derivatives.d1{iz}*tmpOm);
        tmpOm2 = derivatives.d0{iz}\(derivatives.drd{iz}*tmpOm2); % drd! 
        % recopy back to the velocity gradient array
        for IY=iy0:dns.ny; iy=IY+1; jy=IY-iy0+1; jz=field.nzN(iy)+1+IZ; 
            field.domega{iy}(:,2,jz,:) = reshape(tmpOm(jy,:),[3,1,1,dns.nx+1]);
            field.domega2{iy}(:,2,jz,:) = reshape(tmpOm2(jy,:),[3,1,1,dns.nx+1]);
        end
end
