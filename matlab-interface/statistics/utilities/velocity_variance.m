function [field]=velocity_variance(dns, field, derivatives, dc)
%% Calculates the variance in dependence of the radial position field.y(IY)
field.var = complex(zeros(3,3,1,dns.ny+1),0); %( vel1,vel2, derivative(not needed here, thus 1),:)
for IY=0:dns.ny
    iy=IY+1;
    for iV1=1:3 % ...for all velocity components
        for iV2=iV1:3 % ...for all velocity components
            field.var(iV1,iV2,1,iy)=      2*sum(sum(real( field.V{iy}(iV1,:,2:dns.nx+1).*conj(field.V{iy}(iV2,:,2:dns.nx+1)) )))  + ...
                sum(sum(real( field.V{iy}(iV1,1:field.nzN(iy),1).*conj(field.V{iy}(iV2,1:field.nzN(iy),1)) ))) + ...
                sum(sum(real( field.V{iy}(iV1,field.nzN(iy)+2:2*field.nzN(iy)+1,1).*conj(field.V{iy}(iV2,field.nzN(iy)+2:2*field.nzN(iy)+1,1)) )));
        end
    end
end
