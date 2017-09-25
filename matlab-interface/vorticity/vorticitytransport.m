function   [field] = vorticitytransport(dns, field, derivatives, dc)
%This function calculates the time independent terms of the vorticity
%transport equation: As they are:
%field.SR: stretching and reorientation
%field.CV: convection
%field.DIF: viscous diffusion
%the signs are set that all terms are on the right hand side of the
%equation
% all terms are multiplied by to avoid singularities
%% Needed spatial derivatives
[field] = velocity_gradient(dns, field, derivatives, dc)
[field] = vorticity(dns, field, derivatives, dc)
[field] = vorticity_gradient(dns, field, derivatives, dc)
%derivatives are. domega_i/dz, domega_i/dr, domega_i/dphi
%Second spatial derivatives are omega_zz,drd (omega), omega_phiphi

%% Alocate terms of the vorticty transport equation
field.SR=cell(dns.ny+1,1);
field.CV=cell(dns.ny+1,1);
field.DIF=cell(dns.ny+1,1);
%% RHS of the vorticty transport equation
for IY=0:dns.ny
    iy=IY+1;
    field.SR{iy}=complex(zeros(3,2*field.nzN(iy)+1,dns.nx+1),2*dns.nxd); %order z,r,phi component
    field.CV{iy}=complex(zeros(3,2*field.nzN(iy)+1,dns.nx+1),2*dns.nxd); %order z,r,phi component
    field.DIF{iy}=complex(zeros(3,2*field.nzN(iy)+1,dns.nx+1),2*dns.nxd); %order z,r,phi component
    
    tmpdUd=complex(zeros(3,3,field.nzd(iy),2*dns.nxd),0);
    tmpV =complex(zeros(1,3,field.nzd(iy),2*dns.nxd),0);
    tmpOmega =complex(zeros(1,3,field.nzd(iy),2*dns.nxd),0);
    tmpdOmega =complex(zeros(3,3,field.nzd(iy),2*dns.nxd),0);

    tmpV=plane_ift(reshape(field.V{iy},[1,3,2*field.nzN(iy)+1,dns.nx+1]),tmpV);
    tmpdUd=plane_ift(field.dU{iy},tmpdUd);
    tmpOmega=plane_ift(reshape(field.omega{iy},[1,3,2*field.nzN(iy)+1,dns.nx+1]),tmpOmega);
    tmpdOmega=plane_ift(reshape(field.domega{iy},[3,3,2*field.nzN(iy)+1,dns.nx+1]),tmpdOmega);
    tmpOmega_fourier=reshape(field.omega{iy},[1,3,2*field.nzN(iy)+1,dns.nx+1]);
   
    % Stretching and reorientation * r^2 
    field.SR{iy}(2,:,:)=plane_fft( ...
                             tmpOmega(1,2,:,:).*tmpdUd(2,2,:,:)*field.y(iy)^2 ...
                            +tmpOmega(1,3,:,:).*tmpdUd(2,3,:,:)*field.y(iy) ...
                            +tmpOmega(1,1,:,:).*tmpdUd(2,1,:,:)*field.y(iy)^2 ...
                            -tmpOmega(1,3,:,:).*tmpV(1,3,:,:)*field.y(iy), ...
                            reshape(field.SR{iy}(2,:,:),[1,1,2*field.nzN(iy)+1,dns.nx+1]));
    field.SR{iy}(3,:,:)=plane_fft( ...
                             tmpOmega(1,2,:,:).*tmpdUd(3,2,:,:)*field.y(iy)^2 ...
                            +tmpOmega(1,3,:,:).*tmpdUd(3,3,:,:)*field.y(iy) ...
                            +tmpOmega(1,1,:,:).*tmpdUd(3,1,:,:)*field.y(iy)^2 ...
                            +tmpOmega(1,3,:,:).*tmpV(1,2,:,:)*field.y(iy), ...
                            reshape(field.SR{iy}(3,:,:),[1,1,2*field.nzN(iy)+1,dns.nx+1]));
    field.SR{iy}(1,:,:)=plane_fft( ...
                             tmpOmega(1,2,:,:).*tmpdUd(1,2,:,:)*field.y(iy)^2 ...
                            +tmpOmega(1,3,:,:).*tmpdUd(1,3,:,:)*field.y(iy) ...
                            +tmpOmega(1,1,:,:).*tmpdUd(1,1,:,:)*field.y(iy)^2, ...
                            reshape(field.SR{iy}(1,:,:),[1,1,2*field.nzN(iy)+1,dns.nx+1]));    
    % Convection * r^2 '*(-1)' to put it on the rhs of the equation
    field.CV{iy}(2,:,:)=(-1)*plane_fft( ...
                             tmpV(1,2,:,:).*tmpdOmega(2,2,:,:)*field.y(iy)^2 ...
                            +tmpV(1,3,:,:).*tmpdOmega(2,3,:,:)*field.y(iy) ...
                            +tmpV(1,1,:,:).*tmpdOmega(2,1,:,:)*field.y(iy)^2 ...
                            -tmpV(1,3,:,:).*tmpOmega(1,3,:,:)*field.y(iy), ...
                            reshape(field.CV{iy}(2,:,:),[1,1,2*field.nzN(iy)+1,dns.nx+1]));
    field.CV{iy}(3,:,:)=(-1)*plane_fft( ...
                             tmpV(1,2,:,:).*tmpdOmega(3,2,:,:)*field.y(iy)^2 ...
                            +tmpV(1,3,:,:).*tmpdOmega(3,3,:,:)*field.y(iy) ...
                            +tmpV(1,1,:,:).*tmpdOmega(3,1,:,:)*field.y(iy)^2 ...
                            +tmpV(1,3,:,:).*tmpOmega(1,2,:,:)*field.y(iy), ...
                            reshape(field.CV{iy}(3,:,:),[1,1,2*field.nzN(iy)+1,dns.nx+1]));
    field.CV{iy}(1,:,:)=(-1)*plane_fft( ...
                             tmpV(1,2,:,:).*tmpdOmega(1,2,:,:)*field.y(iy)^2 ...
                            +tmpV(1,3,:,:).*tmpdOmega(1,3,:,:)*field.y(iy) ...
                            +tmpV(1,1,:,:).*tmpdOmega(1,1,:,:)*field.y(iy)^2, ...
                            reshape(field.CV{iy}(1,:,:),[1,1,2*field.nzN(iy)+1,dns.nx+1]));              
    % Diffusion * r^2 
    field.omega{iy}=reshape(field.omega{iy},[1,3,2*field.nzN(iy)+1,dns.nx+1]);
    field.DIF{iy}(2,:,:)=1/dns.Re*squeeze( ...
                         field.domega2{iy}(2,2,:,:)*field.y(iy) ...
                        +field.domega2{iy}(2,3,:,:) ...
                        +field.domega2{iy}(2,1,:,:)*field.y(iy)^2 ...
                        -tmpOmega_fourier(1,2,:,:) ...
                        -2*field.domega{iy}(3,3,:,:));
    field.DIF{iy}(3,:,:)=1/dns.Re*squeeze( ...
                         field.domega2{iy}(3,2,:,:)*field.y(iy) ...
                        +field.domega2{iy}(3,3,:,:) ...
                        +field.domega2{iy}(3,1,:,:)*field.y(iy)^2 ...
                        -tmpOmega_fourier(1,3,:,:) ...
                        +2*field.domega{iy}(2,3,:,:));                
    field.DIF{iy}(1,:,:)=1/dns.Re*squeeze( ...
                         field.domega2{iy}(1,2,:,:)*field.y(iy) ...
                        +field.domega2{iy}(1,3,:,:) ...
                        +field.domega2{iy}(1,1,:,:)*field.y(iy)^2);

%   Correct by deviding by r.^2
%     field.DIF{iy}=field.DIF{iy}/(field.y(iy)^2);
%     field.SR{iy}=field.SR{iy}/(field.y(iy)^2);
%     field.CV{iy}=field.CV{iy}/(field.y(iy)^2);

end