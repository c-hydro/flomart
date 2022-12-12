function [map,nr,nc,xllcorner,yllcorner,dx] = ReadAsciiRaster(namf)
%{

    map = matrix with NaN where values match "NODATA_value"

%}

map = arcgridread(namf);
[nr,nc] = size(map);

fid = fopen(namf);
for iline=1:6
str = fgetl(fid);
[fname,val] = strtok(str);
switch lower(fname)
    case 'xllcorner'
        xllcorner = str2double(val);
    case 'yllcorner'
        yllcorner = str2double(val);
    case 'cellsize'
        dx = str2double(val);
end
end
fclose(fid);
