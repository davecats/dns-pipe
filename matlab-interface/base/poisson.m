function   [field] = poisson(dns, field, derivatives, dc)

% 
% This program solves the Poisson equation for the pressure
% 


% Allocate pressure field and velocity gradient
field.p=cell(dns.ny+1,1); field.dU=cell(dns.ny+1,1);
for IY=0:dns.ny
    iy=IY+1; field.p{iy} =complex(zeros(    2*field.nzN(iy)+1,dns.nx+1),0);
             field.dU{iy}=complex(zeros(3,3,2*field.nzN(iy)+1,dns.nx+1),0);
end

% Compute velocity gradient
for IY=0:dns.ny
    iy=IY+1;
    kx = repmat(reshape((0:dns.nx)*dns.alfa0          ,[1,1,dns.nx+1]         ),3,2*field.nzN(iy)+1,1);
    kz = repmat(reshape((-field.nzN(iy):field.nzN(iy)),[1,2*field.nzN(iy)+1,1]),3,1                ,dns.nx+1);
    field.dU{iy}(:,1,:,:) = 1j.*kx.*field.V{iy}(:,:,:);
    field.dU{iy}(:,3,:,:) = 1j.*kz.*field.V{iy}(:,:,:);
end
for IZ=-dns.nz:dns.nz; iz=dns.nz+1+IZ;
    iy0=field.iy0(IZ+dns.nz+1);
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


% Compute right hand side (RHS)
% 
% for iy=1+(1:dns.ny-1)
%     tmpdUd=complex(zeros(3,3,field.nzd(iy),2*dns.nxd),0);
%     tmpVd =complex(zeros(1,3,field.nzd(iy),2*dns.nxd),0);
%     % Backward Fourier transform
%     tmpVd=plane_ift(reshape(field.V{iy},[1,3,2*field.nzN(iy)+1,dns.nx+1]),tmpVd,[4,3]);
%     tmpdUd=plane_ift(field.dU{iy},tmpdUd,[4,3]);
%     %
%     RHS = squeeze( -tmpdUd(2,2,:,:).^2 -(tmpdUd(3,3,:,:)/field.y(iy)).^2 -tmpdUd(1,1,:,:).^2 ...
%                    -2*tmpdUd(2,3,:,:).*tmpdUd(3,2,:,:)/field.y(iy)                           ...
%                    -2*tmpdUd(2,1,:,:).*tmpdUd(1,2,:,:)                                       ...
%                    -2*tmpdUd(3,1,:,:).*tmpdUd(1,3,:,:)/field.y(iy)                           ...
%                    -(tmpVd(1,2,:,:)/field.y(iy)).^2                                            ...
%                    -2*tmpVd(1,2,:,:).*tmpdUd(3,3,:,:)/(field.y(iy)^2)                          ...
%                    +2*tmpVd(1,3,:,:).*tmpdUd(3,2,:,:)/field.y(iy)                              );
%  end
% 
% % Solve Neumann problem for mode (0,0)
% dpdy = calcdpdy(0,0,dns,derivatives,field);
% A = complex(zeros(dns.ny+1,dns.ny+1),0);



% Solve Dirichlet problem for modes but (0,0)

end

% Compute RHS of Poisson equation
function dpdy = calcdpdy(IX,IZ,dns,derivatives,field)
    dpdy=0; ix=IX+1; iz=field.nzN(end)+1+IZ;
    for i=-1:1
      j=dns.ny+i;
      pn = pn -(1/Re/field.y(end))*(derivatives.drd(dns.ny+1,j)*field.V{j}(2,iz,ix));
    end
end

function pn = calcpn(IX,IZ,dns,derivatives,field)
    pn=0; ix=IX+1; iz=field.nzN(end)+1+IZ;
    for i=-1:1
      j=dns.ny+i;
      pn = pn -(1/Re)*(derivatives.drd(dns.ny+1,j)*1j*(iz*field.V{j}(3,iz,ix) + ix*dns.alfa0*field.V{iz}(1,iz,ix)));
    end
end