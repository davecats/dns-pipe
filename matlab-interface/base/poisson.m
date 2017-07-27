function   [field] = poisson(dns, field, derivatives, dc)

% 
% This program solves the Poisson equation for the pressure
% 


%% Allocate pressure field
field.p=cell(dns.ny+1,1);
for IY=0:dns.ny
    iy=IY+1; field.p{iy} =complex(zeros(    2*field.nzN(iy)+1,dns.nx+1),0);
end

%% Compute right hand side (RHS)
for iy=1+(1:dns.ny-1)
    tmpdUd=complex(zeros(3,3,field.nzd(iy),2*dns.nxd),0);
    tmpVd =complex(zeros(1,3,field.nzd(iy),2*dns.nxd),0);
    % Backward Fourier transform
    tmpVd=plane_ift(reshape(field.V{iy},[1,3,2*field.nzN(iy)+1,dns.nx+1]),tmpVd);
    tmpdUd=plane_ift(field.dU{iy},tmpdUd);
    %
    field.p{iy}(:,:)= plane_fft(-tmpdUd(2,2,:,:).^2 -(tmpdUd(3,3,:,:)/field.y(iy)).^2 -tmpdUd(1,1,:,:).^2 ...
                                -2*tmpdUd(2,3,:,:).*tmpdUd(3,2,:,:)/field.y(iy)                           ...
                                -2*tmpdUd(2,1,:,:).*tmpdUd(1,2,:,:)                                       ...
                                -2*tmpdUd(3,1,:,:).*tmpdUd(1,3,:,:)/field.y(iy)                           ...
                                -(tmpVd(1,2,:,:)/field.y(iy)).^2                                          ...
                                -2*tmpVd(1,2,:,:).*tmpdUd(3,3,:,:)/(field.y(iy)^2)                        ...
                                +2*tmpVd(1,3,:,:).*tmpdUd(3,2,:,:)/field.y(iy),                           ...
                                reshape(field.p{iy}(:,:),[1,1,2*field.nzN(iy)+1,dns.nx+1]));                          
end

%% Solve Dirchlet problem for mode (0,0): mean pressure set to 0 at the wall
A = complex(zeros(dns.ny+1,dns.ny+1),0); tmpp=complex(zeros(dns.ny+1,1),0);
% Regularity condition at the axis
A(1,1:3) = dc(1,1,1:3); 
% Drd operator 
A(1+(1:dns.ny-1),:)=derivatives.drd{dns.nz+1}(1+(1:dns.ny-1),:)./repmat(field.y(1+(1:dns.ny-1)),1,dns.ny+1);
% Dirichlet boundary condition at the wall
A(end,end)=1; %A(dns.ny+1,:) = derivatives.d1{dns.nz+1}(dns.ny+1,:);  field.p{dns.ny+1}(dns.nz+1,1) = calcdpdy(0,0,dns,derivatives,field);
% Finally solve
for iy=1:dns.ny+1; tmpp(iy)=field.p{iy}(field.nzN(iy)+1,1); end
tmpp=A\tmpp;
for iy=1:dns.ny+1; field.p{iy}(field.nzN(iy)+1,1)=tmpp(iy); end

%% Solve Dirichlet problem for all other modes
for IZ=-dns.nz:dns.nz; iz=dns.nz+1+IZ;
    iy0=field.iy0(iz); s=dns.ny-iy0+1; y=repmat(field.y(1+(iy0+1:dns.ny-1)),1,s-2); kz2=IZ*IZ./y.^2;
    B=zeros(s-2,s-2);  tmpp=complex(zeros(s-2,dns.nx+1),0);
    B(:,:)=derivatives.drd{iz}(2:s-1,2:s-1)./y -kz2;
    idx=(double(IZ==0):dns.nx);
    % Copy field.p -> tmpp
    for IY=iy0+1:dns.ny-1; iy=IY+1; jz=field.nzN(iy)+IZ+1; tmpp(IY-iy0,idx+1)=field.p{iy}(jz,idx+1); end
    % Solve
    for IX=idx; ix=IX+1;
        %
        kx2=(IX*dns.alfa0)^2;
        % calcpn
        field.p{dns.ny+1}(iz,ix) = calcpn(IX,IZ,dns,derivatives,field);
        % 2nd derivative in x
        A=B-kx2;
        % Regularity condition
        A(1,1:2) = A(1,1:2) - (derivatives.drd{iz}(2,1)./y(1) -kz2(1)-kx2)*reshape(dc(abs(IZ)+1,1,2:3),[1,2]);
        % Dirichlet boundary condition at the wall
        tmpp(:,ix) = tmpp(end,ix)-(derivatives.drd{iz}(end-1,end)./y(end) -kz2(end)-kx2)*field.p{dns.ny+1}(iz,ix);
        % Finally solve
        tmpp(:,ix)=A\tmpp(:,ix);
    end
    % Recover pressure with regularity
    field.p{iy0+1}(field.nzN(iy0+1)+IZ+1,idx+1)=-dc(abs(IZ)+1,1,2)*tmpp(2,idx+1)-dc(abs(IZ)+1,1,3)*tmpp(3,idx+1);
    % Copy field.p <- tmpp
    for IY=iy0+1:dns.ny-1; iy=IY+1; jz=field.nzN(iy)+IZ+1; field.p{iy}(jz,idx+1)=tmpp(IY-iy0,idx+1); end
end

end

function pn = calcpn(IX,IZ,dns,derivatives,field)
    pn=0; ix=IX+1; iz=field.nzN(end)+1+IZ; k2=-IZ*IZ-(IX*dns.alfa0)^2;
    for i=-1:1; j=dns.ny+i;
      pn = pn -(1/dns.Re/k2)*(derivatives.drd{iz}(end,end+i-1)*1j*(IZ*field.V{j}(3,iz,ix) + IX*dns.alfa0*field.V{j}(1,iz,ix)));
    end
end