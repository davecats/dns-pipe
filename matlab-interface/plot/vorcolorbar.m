function [cmap]=vorcolorbar(datatoplot)
datatoplot=datatoplot(:);
% Create colorbar where 0 is white, -2 is blue, 2 is red
            mindata=min(real(datatoplot));
            maxdata=max(real(datatoplot));
            delta_color=maxdata-mindata;
            delta_bar=delta_color/300;
            k_zero=300-round(max(real(datatoplot))/delta_bar);
            if max(real(datatoplot))>2
            k_2=300-round((max(real(datatoplot))-2)/delta_bar);
            else
                k_2=300;
            end
            if min(real(datatoplot))<-2
            k_min2=300-round((max(real(datatoplot))+2)/delta_bar);
            else 
               k_min2=0; 
            end
            dblue=[0,0,0.5];
            blue = [0, 0, 1];
            lblue=[0.8, 0.8, 1]
            red = [1, 0, 0];
            lred=[1, 0.8, 0.8];
            white=[1,1,1];
            dred=[0.5,0,0];
            a=[linspace(dblue(1),blue(1),k_min2);linspace(dblue(2),blue(2),k_min2);linspace(dblue(3),blue(3),k_min2)];
            b=[linspace(blue(1),lblue(1),floor((k_zero-k_min2-21/2)));linspace(blue(2),lblue(2),floor((k_zero-k_min2-21/2)));linspace(blue(3),lblue(3),floor((k_zero-k_min2-21/2)))];
            c=[linspace(lblue(1),white(1),7);linspace(lblue(2),white(2),7);linspace(lblue(3),white(3),7)];
            d=ones(3,7);
            e=[linspace(white(1), lred(1),7);linspace(white(2), lred(2),7);linspace(white(3), lred(3),7)];
            f=[linspace(lred(1),red(1),floor((k_2-k_zero-21/2)));linspace(lred(2),red(2),floor((k_2-k_zero-21/2)));linspace(lred(3),red(3),floor((k_2-k_zero-21/2)))];
            g=[linspace(red(1),dred(1),300-k_2);linspace(red(2),dred(2),300-k_2);linspace(red(3),dred(3),300-k_2)];
            cmap=[a,b,c,d,e,f,g];
            cmap=cmap';