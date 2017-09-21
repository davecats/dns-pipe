function field=velocity_mean(dns, field, derivatives, dc)
%% Calculates the mean velocities and its radial derivatives
field.Vm = complex(zeros(3,dns.ny),0);

for IY= 0:dns.ny
    iy= IY+1;
    for iV=1:3
        field.Vm(iV,iy) = field.V{iy}(iV,field.nzN(iy)+1,1);  
    end    
end

field.dVm=(derivatives.d0{dns.nz+1}\(derivatives.d1{dns.nz+1}*field.Vm(:,:)'))';
