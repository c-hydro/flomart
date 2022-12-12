clc; clear;


%%
%% COMPUTE Qindex if not available:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Period with available discharges must be >= 10 years !!


%% INPUTS:
%%%%%%%%%%
% path with the choice grid layer:
sPathGridGeoData      = '/home/matteo/Documents/CIMA_projects/RT_FloodMapping/data/data_Marche_Chienti/data_static/PREPARATION/gridded_marche';
% path containing discharge values for each cell of the grid for a time:
sPathQmaps            = [path_preparation_data, '/discharge_2010_2020_NEW'];
sDateFrom             = '201101012300';
sDateTo               = '202012312300';
dt                    =  1;  % time step of files (days)
path_preparation_data = '/home/matteo/Documents/CIMA_projects/RT_FloodMapping/data/data_Marche_Chienti/data_static/PREPARATION';
domain_name           = 'Chienti';




%% START:
%%%%%%%%%
sFileName_choice = [sPathGridGeoData, '/marche.choice.txt'];
[a2dMap_choice, a2dCoord_choice] = arcgridread(sFileName_choice);
a2iChoice = a2dMap_choice;
a2iChoice(isnan(a2iChoice)) = -1;  %replace all NaN with -1
[iNRows,iNCols]= size(a2iChoice);
iNoData        = -9999; % valore per dati mancanti nelle mappe Netcdf
% define startuing and ending time of simulations:
nDateFrom      = datenum(sDateFrom,'yyyymmddHHMM');
nDateTo        = datenum(sDateTo,'yyyymmddHHMM');
%initialise current time:
nNow           = nDateFrom;
%initialise matrixes with discharge:
a3dMapQ        = zeros([iNRows,iNCols,10]);
maxQannual     = zeros([iNRows,iNCols]);
% initialise counter:
iCountDay      = 0;
iCountYear     = 1;

while nNow<=nDateTo
    iCountDay = iCountDay + 1;
    sDate = datestr(nNow,'yyyymmddHHMM');
    disp(sDate);
    % extract month, day, hour from file name:  
    iYear   = str2double(sDate(1:4));
    iMonth  = str2double(sDate(5:6));
    iDay    = str2double(sDate(7:8));
    iHour   = str2double(sDate(9:10));
    sPathNow = [sPathQmaps,'/',sDate(1:4),'/',sDate(5:6),'/',sDate(7:8),'/'];
    try
        sFileNameMap = ['hmc.output-grid.',sDate,'.nc.gz'];   
        a2dMap = Continuum_getMap_NC(sPathNow, sFileNameMap,'Discharge');  %istantaneous Discharge (m)
        a2dMap(a2dMap==iNoData) = NaN;
        % compute maximum between actual Q(t) and all previous ones of the
        % year:
        a3dMapQ(:,:,iCountYear) = max(a3dMapQ(:,:,iCountYear), a2dMap);
        %a3dMapQ(:,:,iCountYear) = arrayfun(@(x,y) max(x(:),y(:)), a3dMapQ(:,:,iCountYear), a2dMap);
    catch
        disp('problem with format of netcdf file!! skip!');
    end

    if iMonth==12 & iDay ==31 & iHour == 23
            maxQannual(:,:,iCountYear) =  a3dMapQ(:,:,iCountYear);
            %display(maxQannual);
            display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
            iCountDay = 0;
            % go to next year:
            iCountYear = iCountYear + 1;
    end   
    % pass to next time step:
    nNow = datenum(sDate,'yyyymmddHHMM')+dt;
end



%% Calcolo mappa di portata massima annuale media su deici anni:
maxQannual_mean(:,:) = zeros([iNRows,iNCols]);
for iCountYear=1:10
    tmp = maxQannual_mean(:,:);
    %tmp(OutOfCN) = NaN;
    maxQannual_mean(:,:) = tmp + maxQannual(:,:,iCountYear);
end
maxQannual_mean(:,:) = maxQannual_mean(:,:)./10;
a2dQindice = maxQannual_mean;




%% Save obtained Qindex layer:
save([path_preparation_data, '/Qindex_bis_',domain_name,'.mat'], 'a2dQindice');