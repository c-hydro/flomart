%{

    Trova coordinate matrice di una serie di punti coerentemente con 
    l'ascii raster indicato tragli input e non considerando tutte le
    celle che sono NaN nella matrice 'masknan'


IT-Tor prato: 45.844349, 7.578145 
IT-Trf bosco larice: 45.823761, 7.560866

    *** Se, per esempio, masknan è da ricavare da un raster choice
    mettendo NaN in tutte le celle che NON sono canale:
    	tmp = arcgridread(nomefilechoice);
    	masknan = nan(size(tmp));
    	masknan(tmp==1) = 333;


    x_points = [long1, long2, long3];
    y_points = [lat1, lat2, lat3];
    namfascii = '/home/pippo/nomefileraster.txt';
    masknan = matrice matlab con tante righe e tante colonne quante
    "namfascii" e avere NaN in tutti i pixel da escludere, e un valore
    finito qualsiasi nelle celle candidate***

%}


function [rr_points,cc_points] = Get_MatrixCoordinates_FromPoints_NoNaN(x_points,y_points,namfascii,masknan)



%% DOING
[map,nr,nc,xll,yll,dx] = ReadAsciiRaster(namfascii);

xcenter = (xll+dx/2)+(0:nc-1)*dx;
ycenter = (yll+dx/2)+(0:nr-1)*dx;

[xx,yy] = meshgrid(xcenter,ycenter);
yy = flipud(yy);
 
xx(isnan(masknan)) = NaN;
yy(isnan(masknan)) = NaN;

[kmin,rr_points,cc_points] = deal(nan(length(x_points),1));
for i=1:length(x_points)
    dist2 = (xx-x_points(i)).^2 + (yy-y_points(i)).^2;
    [valmin,kmin(i)] = min(dist2(:));
    [rr_points(i),cc_points(i)] = ind2sub([nr,nc],kmin(i));
end
 
 
 
