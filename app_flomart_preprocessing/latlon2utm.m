function [x,y] = latlon2utm(lat, lon)

% Converte gradi lat/lon in coordinate metriche utm 
%
% INPUT:
%    lat e lon
% OUTPUT:
%    coordinate utm in metri

% get the utm zone from the coordinates :
p1 = [lat,lon];
z1 = utmzone(p1);
disp(z1)
% get the geoid of this zone and  construct the projection 
% structure using the following functions:
[ellipsoid,estr] = utmgeoid(z1);
utmstruct = defaultm('utm');
utmstruct.zone = z1;
utmstruct.geoid = ellipsoid;
utmstruct = defaultm(utmstruct);
% use mfwdtran to convert coordinates lat lon to utm meters: 
[x,y] = mfwdtran(utmstruct,lat,lon);
