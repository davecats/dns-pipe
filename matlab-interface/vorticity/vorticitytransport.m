function   [field] = vorticitytransport(dns, field, derivatives, dc)
%This function calculates the time independent terms of the vorticity
%transport equation: As they are:
%field.SR: stretching and reorientation
%field.CV: convection
%field.DIF: viscous diffusion
%the signs are set that all terms are on the right hand side of the
%equation
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
% dada=2;
% for da=1:5
%     field.rhs{da}=complex(zeros(size(field.rhs{da})));
% 
% end
return;

for IX=0:dns.nx
    ix=IX+1;
    for IZ=-field.nzN(1):field.nzN(1); iz=IZ+field.nzN(1)+1;
            field.DIF{1}(1:3,iz,ix)=complex(0,0);
            for i=2:3
              jz=IZ+field.nzN(i)+1;
              field.DIF{1}(1,iz,ix) = field.DIF{1}(1,iz,ix) - dc(abs(IZ)+1,2,i).*field.DIF{i}(1,jz,ix);
              field.DIF{1}(2,iz,ix) = field.DIF{1}(2,iz,ix) - dc(abs(IZ)+1,3,i).*field.DIF{i}(2,jz,ix);
              field.DIF{1}(3,iz,ix) = field.DIF{1}(3,iz,ix) - dc(abs(IZ)+1,4,i).*field.DIF{i}(3,jz,ix);
            end
    end
end
[plotit] = total_ift_nat(dns,field,field.omega, dc);
plotwhat.name{1}='rphi';
plotwhat.plane{1}=20;
plotwhat.lvl{1}=NaN;
plotwhat.dim{1}=[1,3];
plotfield(dns, field, plotit,plotwhat)



return; 
sum(sum(real(field.omega{88}(2,:,:))))
sum(sum(real(field.domega{88}(2,2,:,:))))
sum(sum(real(field.dU{88}(2,:,:))))
sum(sum(real(field.SR{88}(1,2,:,:))))


[field.omegar] = total_ift_nat(dns, field, field.omega, dc)


[field.omegar] = total_ift_full(dns, field, field.omega, dc)
plotwhat.name{1}='rphi';
plotwhat.plane{1}=20;
plotwhat.lvl{1}=2;
plotwhat.dim{1}=[1,1];
plotfield_full(dns, field, field.omegar,plotwhat)

%%
plotwhat.name{1}='rphi';
plotwhat.plane{1}=20;
plotwhat.lvl{1}=NaN;
plotwhat.dim{1}=[1,3];
plotfield(dns, field, field.omegar,plotwhat)

return
%%
figure
subplot(1,3,1)
contourf(field.y(1:end),linspace(0,2*pi/dns.alfa0,2*dns.nxd),squeeze(real( field.omegar(1,1:end,11,:)))',1000,'LineStyle','none')
caxis([-3,3])
colorbar
title('\omega_z')
subplot(1,3,2)
contourf(field.y(1:end),linspace(0,2*pi/dns.alfa0,2*dns.nxd),squeeze(real( field.omegar(2,1:end,11,:)))',1000,'LineStyle','none')
caxis([-3,3])
colorbar
title('\omega_r')
subplot(1,3,3)
contourf(field.y(1:end),linspace(0,2*pi/dns.alfa0,2*dns.nxd),squeeze(real( field.omegar(3,1:end,11,:)))',1000,'LineStyle','none')
caxis([-3,3])
colorbar
title('\omega_phi')

%% Dissipation of enstrophy

field.DOmega=cell(dns.ny+1,1);
field.DEr=zeros(3,dns.ny+1,field.nzd(end),2*dns.nxd);
for iy=1+(0:dns.ny)
        field.dOmega{iy}(2,2,:,:)=-squeeze(field.dU{iy}(1,2,:,:))*field.y(iy)^-2 ... 
        +squeeze(field.dU2{iy}(1,6,:,:))/field.y(iy) ...
        -squeeze(field.dU2{iy}(2,4,:,:)); % domega_r/dr
        field.dOmega{iy}(1,1,:,:)=squeeze(field.dU{iy}(3,1,:,:))/field.y(iy)...
            +squeeze(field.dU2{iy}(3,4,:,:))...
            -squeeze(field.dU2{iy}(2,5,:,:))/field.y(iy); %domega_z/dz
        field.dOmega{iy}(3,3,:,:)=squeeze(field.dU2{iy}(2,5,:,:))...
            -squeeze(field.dU2{iy}(1,6,:,:))/field.y(iy); %domega_phi/dphi
        
end



%%

return;
% Compute the continuity equation in fourier space

field.cont=cell(dns.ny+1,1)
for iy=1+(0:dns.ny)
    %u_r/r+
    meanVm(iy)=field.V{iy}(3,field.nzN(iy)+1,1);
    field.cont{iy}=squeeze(field.V{iy}(2,:,:))/field.y(iy)+squeeze(field.dU{iy}(2,2,:,:))+squeeze(1/field.y(iy)*field.dU{iy}(3,3,:,:))+squeeze(field.dU{iy}(1,1,:,:));
    res_cont(iy)=field.cont{iy}(field.nzN(iy)+1,1);
    
    %Mittelwerte: Sollten 0 sein
    mean_term(iy,1)=field.V{iy}(2,field.nzN(iy)+1,1)/field.y(iy);
    mean_term(iy,2)=field.dU{iy}(2,2,field.nzN(iy)+1,1);
    mean_term(iy,3)=1/field.y(iy)*field.dU{iy}(3,3,field.nzN(iy)+1,1);
    mean_term(iy,4)=field.dU{iy}(1,1,field.nzN(iy)+1,1);
    
    % Einzelne Terme der Kontiunitätsgleichung
    term1cont{iy}=squeeze(field.V{iy}(2,:,:))/field.y(iy);
    term2cont{iy}=squeeze(field.dU{iy}(2,2,:,:));
    term3cont{iy}=squeeze(1/field.y(iy)*field.dU{iy}(3,3,:,:));
    term4cont{iy}=squeeze(field.dU{iy}(1,1,:,:));
    
    %Varianz der Kontiunitätsgleichung über phi und z;
    var_conti(iy)= 2*sum(sum(real(field.cont{iy}(:,2:dns.nx+1).*conj(field.cont{iy}(:,2:dns.nx+1)) )))  + ...
                                       sum(sum(real( field.cont{iy}(1:field.nzN(iy),1).*conj(field.cont{iy}(1:field.nzN(iy),1)) ))) + ...
                                       sum(sum(real( field.cont{iy}(field.nzN(iy)+2:2*field.nzN(iy)+1,1).*conj(field.cont{iy}(field.nzN(iy)+2:2*field.nzN(iy)+1,1)) )));
    var_term1(iy)= 2*sum(sum(real(term1cont{iy}(:,2:dns.nx+1).*conj(term1cont{iy}(:,2:dns.nx+1)) )))  + ...
                                       sum(sum(real( term1cont{iy}(1:field.nzN(iy),1).*conj(term1cont{iy}(1:field.nzN(iy),1)) ))) + ...
                                       sum(sum(real( term1cont{iy}(field.nzN(iy)+2:2*field.nzN(iy)+1,1).*conj(term1cont{iy}(field.nzN(iy)+2:2*field.nzN(iy)+1,1)) )));       
    var_term2(iy)= 2*sum(sum(real(term2cont{iy}(:,2:dns.nx+1).*conj(term2cont{iy}(:,2:dns.nx+1)) )))  + ...
                                       sum(sum(real( term2cont{iy}(1:field.nzN(iy),1).*conj(term2cont{iy}(1:field.nzN(iy),1)) ))) + ...
                                       sum(sum(real( term2cont{iy}(field.nzN(iy)+2:2*field.nzN(iy)+1,1).*conj(term2cont{iy}(field.nzN(iy)+2:2*field.nzN(iy)+1,1)) )));       
    var_term3(iy)= 2*sum(sum(real(term3cont{iy}(:,2:dns.nx+1).*conj(term3cont{iy}(:,2:dns.nx+1)) )))  + ...
                                       sum(sum(real( term3cont{iy}(1:field.nzN(iy),1).*conj(term3cont{iy}(1:field.nzN(iy),1)) ))) + ...
                                       sum(sum(real( term3cont{iy}(field.nzN(iy)+2:2*field.nzN(iy)+1,1).*conj(term3cont{iy}(field.nzN(iy)+2:2*field.nzN(iy)+1,1)) )));       
    var_term4(iy)= 2*sum(sum(real(term4cont{iy}(:,2:dns.nx+1).*conj(term4cont{iy}(:,2:dns.nx+1)) )))  + ...
                                       sum(sum(real( term4cont{iy}(1:field.nzN(iy),1).*conj(term4cont{iy}(1:field.nzN(iy),1)) ))) + ...
                                       sum(sum(real( term4cont{iy}(field.nzN(iy)+2:2*field.nzN(iy)+1,1).*conj(term4cont{iy}(field.nzN(iy)+2:2*field.nzN(iy)+1,1)) )));       
   
end
plane=111;
figure
subplot(2,3,1)
contourf(abs(field.cont{plane}),100,'LineStyle','none')
colorbar
title('Konti')
subplot(2,3,2)
contourf(abs(term1cont{plane}),100,'LineStyle','none')
colorbar
title('u_r/r')
subplot(2,3,3)
contourf(abs(term2cont{plane}),100,'LineStyle','none')
colorbar
title('d(u_r)/dr')
subplot(2,3,4)
contourf(abs(term3cont{plane}),100,'LineStyle','none')
colorbar
title('1/r du_{\phi}/d \phi')
subplot(2,3,5)
contourf(abs(term4cont{plane}),100,'LineStyle','none')
colorbar
title('du_z/dz')

figure

plot(field.y(1:end),var_conti,'LineWidth',3)
hold on
plot(field.y(1:end),var_term1)
hold on
plot(field.y(1:end),var_term2)
hold on
plot(field.y(1:end),var_term3)
hold on
plot(field.y(1:end),var_term4)
title('variance')
legend('conti','u_r/r','du_r/dr','du_phi/dphi','du_z/dz')
ylim([0,1])

figure
plot(field.y(1:end),meanVm)
figure
plot(field.y(1:end),res_cont)
hold on
plot(field.y(1:end),mean_term(:,1))
plot(field.y(1:end),mean_term(:,2))
plot(field.y(1:end),mean_term(:,3))
plot(field.y(1:end),mean_term(:,4))


field.ContReal=zeros(dns.ny+1,field.nzd(end),2*dns.nxd);
field.ContTerm1=zeros(dns.ny+1,field.nzd(end),2*dns.nxd);
field.ContTerm2=zeros(dns.ny+1,field.nzd(end),2*dns.nxd);
field.ContTerm3=zeros(dns.ny+1,field.nzd(end),2*dns.nxd);
field.ContTerm4=zeros(dns.ny+1,field.nzd(end),2*dns.nxd);
field.Vr=zeros(dns.ny+1,field.nzd(end),2*dns.nxd);
field.UdReal=zeros(3,3,dns.ny+1,field.nzd(end),2*dns.nxd);
for iy=1+(1:dns.ny-1)
    tmpCont=complex(zeros(1,1,field.nzd(end),2*dns.nxd),0);
    tmpVd =complex(zeros(1,3,field.nzd(end),2*dns.nxd),0);
    tmpdUd=complex(zeros(3,3,field.nzd(end),2*dns.nxd),0);
    tmpdUd=plane_ift(field.dU{iy},tmpdUd,[4,3]);
    field.UdReal(:,:,iy,:,:)=tmpdUd;
    %     % Backward Fourier transform
       tmpVd=plane_ift(reshape(field.V{iy},[1,3,2*field.nzN(iy)+1,dns.nx+1]),tmpVd,[4,3]);
       field.Vr(iy,:,:)=squeeze(tmpVd(1,3,:,:));
       tmpCont=complex(zeros(1,1,field.nzd(end),2*dns.nxd),0);
       tmpCont=plane_ift(reshape(field.cont{iy},[1,1,2*field.nzN(iy)+1,dns.nx+1]),tmpCont,[4,3]);
       field.ContReal(iy,:,:)=squeeze(tmpCont);
       tmpCont=complex(zeros(1,1,field.nzd(end),2*dns.nxd),0);
       tmpCont=plane_ift(reshape(term1cont{iy},[1,1,2*field.nzN(iy)+1,dns.nx+1]),tmpCont,[4,3]);
       field.ContTerm1(iy,:,:)=squeeze(tmpCont);
       tmpCont=complex(zeros(1,1,field.nzd(end),2*dns.nxd),0);
       tmpCont=plane_ift(reshape(term2cont{iy},[1,1,2*field.nzN(iy)+1,dns.nx+1]),tmpCont,[4,3]);
       field.ContTerm2(iy,:,:)=squeeze(tmpCont);
       tmpCont=complex(zeros(1,1,field.nzd(end),2*dns.nxd),0);
       tmpCont=plane_ift(reshape(term3cont{iy},[1,1,2*field.nzN(iy)+1,dns.nx+1]),tmpCont,[4,3]);
       field.ContTerm3(iy,:,:)=squeeze(tmpCont);
       tmpCont=complex(zeros(1,1,field.nzd(end),2*dns.nxd),0);
       tmpCont=plane_ift(reshape(term4cont{iy},[1,1,2*field.nzN(iy)+1,dns.nx+1]),tmpCont,[4,3]);
       field.ContTerm4(iy,:,:)=squeeze(tmpCont);
       
end
figure
contourf(field.y(1:end),linspace(0,2*pi/dns.alfa0,2*dns.nxd),squeeze(real( field.UdReal(2,2,5:end,11,:)))',100,'LineStyle','none')
figure
contourf(field.y(1:end),linspace(0,2*pi/dns.alfa0,2*dns.nxd),squeeze(real(field.Vr(1:end,11,:)))',100,'LineStyle','none')
figure
subplot(1,5,1)
contourf(field.y(1:end),linspace(0,2*pi/dns.alfa0,2*dns.nxd),squeeze(real(field.ContReal(1:end,11,:)))',100,'LineStyle','none')
colorbar
title('Konti')
subplot(1,5,2)
contourf(field.y(1:end),linspace(0,2*pi/dns.alfa0,2*dns.nxd),squeeze(real(field.ContTerm1(1:end,11,:)))',100,'LineStyle','none')
colorbar
title('u_r/r')
subplot(1,5,3)
contourf(field.y(1:end),linspace(0,2*pi/dns.alfa0,2*dns.nxd),squeeze(real(field.ContTerm2(1:end,11,:)))',100,'LineStyle','none')
colorbar
title('d(u_r)/dr')
subplot(1,5,4)
contourf(field.y(1:end),linspace(0,2*pi/dns.alfa0,2*dns.nxd),squeeze(real(field.ContTerm3(1:end,11,:)))',100,'LineStyle','none')
colorbar
title('1/r du_{\phi}/d \phi')
subplot(1,5,5)
contourf(field.y(1:end),linspace(0,2*pi/dns.alfa0,2*dns.nxd),squeeze(real(field.ContTerm4(1:end,11,:)))',100,'LineStyle','none')
colorbar
title('du_z/dz')

[R PHI] = meshgrid(r,phi);
Z = F(R,PHI); % which assumes your function is vectorized
surf(R.*cos(PHI), R.*sin(PHI), Z);


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
% function dpdy = calcdpdy(IX,IZ,dns,derivatives,field)
%     dpdy=0; ix=IX+1; iz=field.nzN(end)+1+IZ;
%     for i=-1:1
%       j=dns.ny+i;
%       pn = pn -(1/Re/field.y(end))*(derivatives.drd(dns.ny+1,j)*field.V{j}(2,iz,ix));
%     end
% end
% 
% function pn = calcpn(IX,IZ,dns,derivatives,field)
%     pn=0; ix=IX+1; iz=field.nzN(end)+1+IZ;
%     for i=-1:1
%       j=dns.ny+i;
%       pn = pn -(1/Re)*(derivatives.drd(dns.ny+1,j)*1j*(iz*field.V{j}(3,iz,ix) + ix*dns.alfa0*field.V{iz}(1,iz,ix)));
%     end
% end