function [field] = write_vtk_binary(dns, field,fieldstowrite,nameoffield,filepath_short,filename,fraction)
%% This function creates a binary vtk file for paraview
% So far it is only able to store vectors. An extension to store scalars
% can easily be implemented

% Variables
%fieldstowrite: A cell array of string names of fields that should be
%written
%nameoffield: Start of the filename of the produced vtk file
%fraction: gives the fraction of data points in phi and axial direction that
%shall be stored. Example: [1,0.5] only plots half the azimuthal extension
%of the pipe but the whole length. fraction should only numbers that lead
%to integer amount of gridpoints

%Example usage:  [field] = write_vtk_binary(dns, field,{'omega','V'},'Omega',filepath_short,filenamelist{iF},[0.3,1])
%%

mkdir(['../../../Results/',filepath_short,'/VTK_BIN/'])
filename=['../../../Results/',filepath_short,'/VTK_BIN/',nameoffield,'_',filename(1:end-4),'.vtk'];
fid = fopen(filename, 'w');

%ASCII file header
sizegrid = [floor(fraction(1)*2*dns.nxd), dns.ny+1, floor(fraction(2)*dns.nzd)];
fprintf(fid, '# vtk DataFile Version 3.0\n');
fprintf(fid, 'VTK from Matlab\n');
fprintf(fid, 'BINARY\n\n');
fprintf(fid, 'DATASET STRUCTURED_GRID\n');
fprintf(fid, ['DIMENSIONS ' num2str(sizegrid(1)) ' ' num2str(sizegrid(2)) ' ' num2str( sizegrid(3)+1) '\n']);
fprintf(fid, ['POINTS ' num2str(sizegrid(1)*sizegrid(2)*(sizegrid(3)+1)) ' float\n']);
fclose(fid);


fid = fopen(filename, 'a');
for iz=0:sizegrid(3)
    tmpy=ones(sizegrid(1),1)*field.y'*sin(2*pi*iz/(dns.nzd));
    tmpz=ones(sizegrid(1),1)*field.y'*cos(2*pi*iz/(dns.nzd));
    tmpx=2*pi/dns.alfa0*[0:sizegrid(1)-1]'/(2*dns.nxd)*ones(dns.ny+1,1)';
    tmpdata=[tmpx(:),tmpy(:),tmpz(:)];
    tmpdata=tmpdata';
    fwrite(fid,tmpdata(:),'float','b');
end

fprintf(fid, ['\nPOINT_DATA ' num2str(sizegrid(1)*sizegrid(2)*(sizegrid(3)+1)) '\n']);
for iN=1:numel(fieldstowrite)
%append another ASCII sub header
fprintf(fid, ['VECTORS ',fieldstowrite{iN}, ' float \n']);
[fieldtowrite] = total_ift_full(dns, field,field.(fieldstowrite{iN}));
writefield=NaN(3,dns.ny+1,sizegrid(3),sizegrid(1));
for iy=1:dns.ny+1
    writefield(:,iy,:,:)=squeeze(fieldtowrite{iy}(1,:,1:sizegrid(3),1:sizegrid(1)));
end

for iz=1:sizegrid(3)
    tmpx=squeeze(writefield(1,:,iz,:))';
    tmpy=squeeze(writefield(2,:,iz,:))';
    tmpz=squeeze(writefield(3,:,iz,:))';
    tmpdata=[tmpx(:),tmpy(:),tmpz(:)];
    tmpdata=tmpdata';
    fwrite(fid,tmpdata(:),'float','b');
end
% The next lines exists to ensure a closed circle in the paraview output pipe.
% Therefore the 0th azimuthal row is replicated as the dns.nzd+1th one.
if sizegrid(3)==dns.nzd
    tmpx=squeeze(writefield(1,:,1,:))';
    tmpy=squeeze(writefield(2,:,1,:))';
    tmpz=squeeze(writefield(3,:,1,:))';
    tmpdata=[tmpx(:),tmpy(:),tmpz(:)];
    tmpdata=tmpdata';
    fwrite(fid,tmpdata(:),'float','b');
else
    writefield=NaN(3,dns.ny+1,sizegrid(3)+1,sizegrid(1));
    for iy=1:dns.ny+1
        writefield(:,iy,:,:)=squeeze(fieldtowrite{iy}(1,:,1:sizegrid(3)+1,1:sizegrid(1)));
    end
    tmpx=squeeze(writefield(1,:,sizegrid(3)+1,:))';
    tmpy=squeeze(writefield(2,:,sizegrid(3)+1,:))';
    tmpz=squeeze(writefield(3,:,sizegrid(3)+1,:))';
    tmpdata=[tmpx(:),tmpy(:),tmpz(:)];
    tmpdata=tmpdata';
    fwrite(fid,tmpdata(:),'float','b');
end


end
fclose(fid);


