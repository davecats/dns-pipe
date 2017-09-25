function plotfield(dns, field, plotfield,plotwhat)
%Plots the fields given as argument: plotfield
%Hereby plotwhat is a struct containing the wished plottypes and its needed
%details.

%% Example Settings
% plotwhat.name{1}='rphi';
% plotwhat.plane{1}=20;
% plotwhat.lvl{1}=NaN;
% plotwhat.dim{1}=[1,1];
% plotwhat.name{2}='rz';
% plotwhat.plane{2}=20;
% plotwhat.lvl{2}=NaN;
% plotwhat.dim{2}=[1,1];
% plotwhat.name{3}='3D';
% plotwhat.plane{3}=20;
% plotwhat.lvl{3}=NaN;
% plotwhat.dim{3}=[1,1];

%%

for pNo=1:numel(plotwhat)
    switch plotwhat.name{pNo}
        case 'rphi'
            scatteredData=[];
            for iy=1:dns.ny(end)+1;
                xtmp=field.y(iy)*sin(2*pi/field.nzd(iy)*[1:field.nzd(iy)]);
                ytmp=field.y(iy)*cos(2*pi/field.nzd(iy)*[1:field.nzd(iy)]);
                ztmp=squeeze(plotfield{iy}(plotwhat.dim{pNo}(1),plotwhat.dim{pNo}(2),:,plotwhat.plane{pNo}));
                scatteredData=[ scatteredData;xtmp(:),ytmp(:),ztmp];
            end
            
            %%
            tri = delaunay(scatteredData(:,1),scatteredData(:,2));
            %             [r,c] = size(tri);
            figure
            
            zz=real(scatteredData(:,3));
            
            %             limit=10;
            %             zz(zz>limit)=limit;
            %             zz(zz<-limit)=-limit;
            %             zz(1)=limit;
            %             zz(2)=-limit;
            h = trisurf(tri, scatteredData(:,1),scatteredData(:,2),zz,'LineStyle','none');
            view([0,0,180])
            
            axis equal
            colorbar
            cmap=colormap;
            cmap=vorcolorbar(zz);
            colormap(cmap)
            
            %%
            figure
            
            h = trisurf(tri, scatteredData(:,1),scatteredData(:,2),zz,'LineStyle','none');
            l = light('Position',[-50 -15 29])
            set(gca,'CameraPosition',[ -5.2075  -11.1087  117.2614])
            lighting phong%gouraud
            shading interp
            %     colormapeditor
            colorbar EastOutside
            axis off
            colormap(cmap)
            
            
        case 'rz'
            % Here always the phi=0 plane is plotted
            ztmp=zeros(dns.ny(end)+1,2*dns.nxd);
            for iy=2:dns.ny(end)+1;
                ztmp(iy,:)=squeeze(plotfield{iy}(plotwhat.dim{pNo}(1),plotwhat.dim{pNo}(2),1,:));
            end
            zz=real(ztmp);
            %             limit=10;
            %             zz(zz>limit)=limit;
            %             zz(zz<-limit)=-limit;
            %             zz(1)=limit;
            %             zz(2)=-limit;
            figure
            contourf(field.y(end-20:end),linspace(0,2*pi/dns.alfa0,2*dns.nxd),zz(end-20:end,:)',1000,'LineStyle','none')
            cmap=vorcolorbar(real(zz));
            colormap(cmap)
            axis equal
            colorbar

            
        case '3D'
            disp('not yet implemented')
            
    end
    
end
