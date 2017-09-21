function field=bl_thinkness(dns, field, derivatives, dc)    
Vm=field.Vm(3,:)';
Vm(1)=0;
derr=(field.y-Vm)./field.y;
derr(derr<0)=0;
field.d99=field.y(max(1,numel(Vm)-sum(derr>0.01)));
field.d99_pos=numel(Vm)-sum(derr>0.01);



%% THIS SHOULD BE EXTENDED BY: d2,d1,etc.

