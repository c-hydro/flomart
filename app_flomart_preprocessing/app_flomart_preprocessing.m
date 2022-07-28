clc; clear;

%%
%%%%%%%%%%%%%%
%   INPUTS   %
%%%%%%%%%%%%%%
% please, if a new case study need to be analysed, then insert the name of the 
% new domain in the list below. Then execute next command to select the
% case study from prompt pop-up:
list_domains = {'Foglia', 'Entella', 'EntellaCompleto', 'Scrivia', 'Caccamo-Grazie', 'Grazie-Foce', 'Polverina-Caccamo', 'Chienti'}
domain_name = strjoin(list_domains([listdlg('ListString', ['Choose case study:', ' ', list_domains])]-2))





if strcmp(domain_name,'Foglia')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % percorso dove si trova questo script matlab (e i suoi moduli) per il preprocessing in Flomart:
    path_code = '/home/matteo/Documents/CIMA_projects/RT_FloodMapping/flomart-2.0.0_test/app_flomart_preprocessing';
    % percorso dove preparare i files statici:
    path_preparation_data = '/home/matteo/Documents/CIMA_projects/RT_FloodMapping/data/data_Marche_Foglia/data_static/PREPARATION';
    % percorso dove cercare dati statici geografici con griglia: 
    %sPathGridGeoData = '/home/matteo/Documents/CIMA_projects/RT_FloodMapping/data/data_Marche_Foglia/data_static/PREPARATION/fp_marche';
    sPathGridGeoData = '/home/matteo/Documents/CIMA_projects/RT_FloodMapping/data/data_Marche_Foglia/data_static/PREPARATION/gridded_marche';
    % Percorso dove salvare il file mat output 'info_{domain_name}.mat':
    sPathOutputInfoMat = '/home/matteo/Documents/CIMA_projects/RT_FloodMapping/data/data_Marche_Foglia/data_static/geo_data';
    % percorso dove si trovano le simulazioni:
    sPathHazardData = '/home/matteo/Documents/CIMA_projects/RT_FloodMapping/data/data_Marche_Foglia/data_static/hazard_data';
    %name_hazardmaps = [domain_name,'_WD_max_Q.tif'];% nome delle hazard maps di partenza senza le tre cifre finali (es:  "Foglia_WD_max_Q001.tif" )
    name_hazardmaps = ['Foglia_WD_max_T'];
    type_hazardmaps ='.tif';  % format of the hazard map files
    % tempo di ritorno al di sotto (strettamente) del quale la mappa di hazard si annulla
    TR_min=1;
    % tempo di ritorno massimo (eventuali valori superiori vengono saturati a questo valore)
    TR_max=500;
    modeSelectionSections = 'hardcoded';  % mode for the selection of control sections ('hardcoded' or 'byprompt')
    % coordinate system EPSG utm: 
    EPSG_domain = '32633';  % Marche Region zone 33 utm
    
    
    
elseif strcmp(domain_name,'Entella')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % percorso dove si trova questo script matlab (e i suoi moduli) per il preprocessing in Flomart:
    path_code = '/home/matteo/Documents/CIMA_projects/RT_FloodMapping/flomart-2.0.0_test/app_flomart_preprocessing';
    % percorso dove preparare i files statici:
    path_preparation_data = '/home/matteo/Documents/CIMA_projects/RT_FloodMapping/data/data_Liguria_Entella/data_static/PREPARATION';
    % percorso dove cercare dati statici geografici con griglia: 
    sPathGridGeoData = '/home/matteo/Documents/CIMA_projects/RT_FloodMapping/data/data_Liguria_Entella/data_static/PREPARATION/grid';
    % Percorso dove salvare il file mat output 'info_{domain_name}.mat':
    sPathOutputInfoMat = '/home/matteo/Documents/CIMA_projects/RT_FloodMapping/data/data_Liguria_Entella/data_static/geo_data';
    % percorso dove si trovano le simulazioni:
    sPathHazardData = '/home/matteo/Documents/CIMA_projects/RT_FloodMapping/data/data_Liguria_Entella/data_static/hazard_data';
    %name_hazardmaps = [domain_name,'_WD_max_Q.tif'];% nome delle hazard maps di partenza senza le tre cifre finali (es:  "Foglia_WD_max_Q001.tif" )
    name_hazardmaps = ['EntellaCompleto_WD_max_T'];
    type_hazardmaps ='.tif';  % format of the hazard map files
    % tempo di ritorno al di sotto (strettamente) del quale la mappa di hazard si annulla
    TR_min=1;
    % tempo di ritorno massimo (eventuali valori superiori vengono saturati a questo valore)
    TR_max=500;
    modeSelectionSections = 'hardcoded';  % mode for the selection of control sections ('hardcoded' or 'byprompt')
    % coordinate system EPSG utm: 
    EPSG_domain = '32632'; % Liguria Region zone 32 utm

    
    
elseif strcmp(domain_name,'Scrivia')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % percorso dove si trova questo script matlab (e i suoi moduli) per il preprocessing in Flomart:
    path_code = '/home/matteo/Documents/CIMA_projects/RT_FloodMapping/flomart-2.0.0_test/app_flomart_preprocessing';
    % percorso dove preparare i files statici:
    path_preparation_data = '/home/matteo/Documents/CIMA_projects/RT_FloodMapping/data/data_Liguria_Scrivia/data_static/PREPARATION';
    % percorso dove cercare dati statici geografici con griglia: 
    sPathGridGeoData = '/home/matteo/Documents/CIMA_projects/RT_FloodMapping/data/data_Liguria_Scrivia/data_static/PREPARATION/grid';
    % Percorso dove salvare il file mat output 'info_{domain_name}.mat':
    sPathOutputInfoMat = '/home/matteo/Documents/CIMA_projects/RT_FloodMapping/data/data_Liguria_Scrivia/data_static/geo_data';
    % percorso dove si trovano le simulazioni:
    sPathHazardData = '/home/matteo/Documents/CIMA_projects/RT_FloodMapping/data/data_Liguria_Scrivia/data_static/hazard_data';
    %name_hazardmaps = [domain_name,'_WD_max_Q.tif'];% nome delle hazard maps di partenza senza le tre cifre finali (es:  "Foglia_WD_max_Q001.tif" )
    name_hazardmaps = ['Scrivia_WD_max_T'];
    type_hazardmaps ='.tif';  % format of the hazard map files
    % tempo di ritorno al di sotto (strettamente) del quale la mappa di hazard si annulla
    TR_min=1;
    % tempo di ritorno massimo (eventuali valori superiori vengono saturati a questo valore)
    TR_max=500;
    modeSelectionSections = 'hardcoded';  % mode for the selection of control sections ('hardcoded' or 'byprompt')
    % coordinate system EPSG utm: 
    EPSG_domain = '32632'; % Liguria Region zone 32 utm

    
    
elseif strcmp(domain_name,'Grazie-Foce')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % percorso dove si trova questo script matlab (e i suoi moduli) per il preprocessing in Flomart:
    path_code = '/home/matteo/Documents/CIMA_projects/RT_FloodMapping/flomart-2.0.0_test/app_flomart_preprocessing';
    % percorso dove preparare i files statici:
    path_preparation_data = '/home/matteo/Documents/CIMA_projects/RT_FloodMapping/data/data_Marche_Chienti/data_static/PREPARATION';
    % percorso dove cercare dati statici geografici con griglia: 
    sPathGridGeoData = '/home/matteo/Documents/CIMA_projects/RT_FloodMapping/data/data_Marche_Chienti/data_static/PREPARATION/gridded_marche';
    % Percorso dove salvare il file mat output 'info_{domain_name}.mat':
    sPathOutputInfoMat = '/home/matteo/Documents/CIMA_projects/RT_FloodMapping/data/data_Marche_Chienti/data_static/geo_data';
    % percorso dove si trovano le simulazioni:
    sPathHazardData = '/home/matteo/Documents/CIMA_projects/RT_FloodMapping/data/data_Marche_Chienti/data_static/hazard_data';
    %name_hazardmaps = [domain_name,'_WD_max_Q.tif'];% nome delle hazard maps di partenza senza le tre cifre finali (es:  "Grazie-Foce_WD_max_T001.tif" )
    name_hazardmaps = ['Grazie-Foce_WD_max_T'];
    type_hazardmaps ='.tif';  % format of the hazard map files
    % tempo di ritorno al di sotto (strettamente) del quale la mappa di hazard si annulla
    TR_min=1;
    % tempo di ritorno massimo (eventuali valori superiori vengono saturati a questo valore)
    TR_max=500;
    modeSelectionSections = 'hardcoded';  % mode for the selection of control sections ('hardcoded' or 'byprompt')
    % coordinate system EPSG utm: 
    EPSG_domain = '32633'; % Marche Region zone 33 utm

    
    
    
elseif strcmp(domain_name,'Caccamo-Grazie')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % percorso dove si trova questo script matlab (e i suoi moduli) per il preprocessing in Flomart:
    path_code = '/home/matteo/Documents/CIMA_projects/RT_FloodMapping/flomart-2.0.0_test/app_flomart_preprocessing';
    % percorso dove preparare i files statici:
    path_preparation_data = '/home/matteo/Documents/CIMA_projects/RT_FloodMapping/data/data_Marche_Chienti/data_static/PREPARATION';
    % percorso dove cercare dati statici geografici con griglia: 
    sPathGridGeoData = '/home/matteo/Documents/CIMA_projects/RT_FloodMapping/data/data_Marche_Chienti/data_static/PREPARATION/gridded_marche';
    % Percorso dove salvare il file mat output 'info_{domain_name}.mat':
    sPathOutputInfoMat = '/home/matteo/Documents/CIMA_projects/RT_FloodMapping/data/data_Marche_Chienti/data_static/geo_data';
    % percorso dove si trovano le simulazioni:
    sPathHazardData = '/home/matteo/Documents/CIMA_projects/RT_FloodMapping/data/data_Marche_Chienti/data_static/hazard_data';
    %name_hazardmaps = [domain_name,'_WD_max_Q.tif'];% nome delle hazard maps di partenza senza le tre cifre finali (es:  "Grazie-Foce_WD_max_T001.tif" )
    name_hazardmaps = ['Caccamo-Grazie_WD_max_T'];
    type_hazardmaps ='.tif';  % format of the hazard map files
    % tempo di ritorno al di sotto (strettamente) del quale la mappa di hazard si annulla
    TR_min=1;
    % tempo di ritorno massimo (eventuali valori superiori vengono saturati a questo valore)
    TR_max=500;
    modeSelectionSections = 'hardcoded';  % mode for the selection of control sections ('hardcoded' or 'byprompt')
    % coordinate system EPSG utm: 
    EPSG_domain = '32633'; % Marche Region zone 33 utm

    
    
elseif strcmp(domain_name,'Polverina-Caccamo')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % percorso dove si trova questo script matlab (e i suoi moduli) per il preprocessing in Flomart:
    path_code = '/home/matteo/Documents/CIMA_projects/RT_FloodMapping/flomart-2.0.0_test/app_flomart_preprocessing';
    % percorso dove preparare i files statici:
    path_preparation_data = '/home/matteo/Documents/CIMA_projects/RT_FloodMapping/data/data_Marche_Chienti/data_static/PREPARATION';
    % percorso dove cercare dati statici geografici con griglia: 
    sPathGridGeoData = '/home/matteo/Documents/CIMA_projects/RT_FloodMapping/data/data_Marche_Chienti/data_static/PREPARATION/gridded_marche';
    % Percorso dove salvare il file mat output 'info_{domain_name}.mat':
    sPathOutputInfoMat = '/home/matteo/Documents/CIMA_projects/RT_FloodMapping/data/data_Marche_Chienti/data_static/geo_data';
    % percorso dove si trovano le simulazioni:
    sPathHazardData = '/home/matteo/Documents/CIMA_projects/RT_FloodMapping/data/data_Marche_Chienti/data_static/hazard_data';
    %name_hazardmaps = [domain_name,'_WD_max_Q.tif'];% nome delle hazard maps di partenza senza le tre cifre finali (es:  "Grazie-Foce_WD_max_T001.tif" )
    name_hazardmaps = ['Polverina-Caccamo_WD_max_T'];
    type_hazardmaps ='.tif';  % format of the hazard map files
    % tempo di ritorno al di sotto (strettamente) del quale la mappa di hazard si annulla
    TR_min=1;
    % tempo di ritorno massimo (eventuali valori superiori vengono saturati a questo valore)
    TR_max=500;
    modeSelectionSections = 'hardcoded';  % mode for the selection of control sections ('hardcoded' or 'byprompt')
    % coordinate system EPSG utm: 
    EPSG_domain = '32633'; % Marche Region zone 33 utm
    
elseif strcmp(domain_name,'Chienti')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % percorso dove si trova questo script matlab (e i suoi moduli) per il preprocessing in Flomart:
    path_code = '/home/matteo/Documents/CIMA_projects/RT_FloodMapping/flomart-2.0.0_test/app_flomart_preprocessing';
    % percorso dove preparare i files statici:
    path_preparation_data = '/home/matteo/Documents/CIMA_projects/RT_FloodMapping/data/data_Marche_Chienti/data_static/PREPARATION';
    % percorso dove cercare dati statici geografici con griglia: 
    sPathGridGeoData = '/home/matteo/Documents/CIMA_projects/RT_FloodMapping/data/data_Marche_Chienti/data_static/PREPARATION/gridded_marche';
    % Percorso dove salvare il file mat output 'info_{domain_name}.mat':
    sPathOutputInfoMat = '/home/matteo/Documents/CIMA_projects/RT_FloodMapping/data/data_Marche_Chienti/data_static/geo_data';
    % percorso dove si trovano le simulazioni:
    sPathHazardData = '/home/matteo/Documents/CIMA_projects/RT_FloodMapping/data/data_Marche_Chienti/data_static/hazard_data';
    %name_hazardmaps = [domain_name,'_WD_max_Q.tif'];% nome delle hazard maps di partenza senza le tre cifre finali (es:  "Grazie-Foce_WD_max_T001.tif" )
    name_hazardmaps = ['Chienti_WD_max_T'];
    type_hazardmaps ='.tif';  % format of the hazard map files
    % tempo di ritorno al di sotto (strettamente) del quale la mappa di hazard si annulla
    TR_min=1;
    % tempo di ritorno massimo (eventuali valori superiori vengono saturati a questo valore)
    TR_max=500;
    modeSelectionSections = 'hardcoded';  % mode for the selection of control sections ('hardcoded' or 'byprompt')
    % coordinate system EPSG utm: 
    EPSG_domain = '32633'; % Marche Region zone 33 utm

end













%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read geographical gridded information:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if strcmp(domain_name,'Foglia')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% names of files with gridded info (choice, dem, ...):
    %MARCHE:
%     sFileName_choice = [sPathGridGeoData, '/MarcheDomain.choice.txt'];
%     sFileName_area = [sPathGridGeoData,'/MarcheDomain.area.txt'];
%     sFileName_cell = [sPathGridGeoData,'/MarcheDomain.areacell.txt'];
%     sFileName_dem = [sPathGridGeoData,'/MarcheDomain.dem.txt'];
%     sFileName_lon = [sPathGridGeoData,'/MarcheDomain.lon.txt'];
%     sFileName_lat = [sPathGridGeoData,'/MarcheDomain.lat.txt'];
%     sFileName_pnt = [sPathGridGeoData,'/MarcheDomain.pnt.txt'];
    sFileName_choice = [sPathGridGeoData, '/marche.choice.txt'];
    sFileName_area = [sPathGridGeoData, '/marche.area.txt'];
    sFileName_cell = [sPathGridGeoData, '/marche.areacell.txt'];
    sFileName_dem = [sPathGridGeoData, '/marche.dem.txt'];
    sFileName_lon = [sPathGridGeoData, '/marche.lon.txt'];
    sFileName_lat = [sPathGridGeoData, '/marche.lat.txt'];
    sFileName_pnt = [sPathGridGeoData, '/marche.pnt.txt'];
    
    %% DEFINING VARIABLES OF FILE .MAT:
    % Choice:
    [a2dMap_choice, a2dCoord_choice] = arcgridread(sFileName_choice);
    a2iChoice = a2dMap_choice;
    a2iChoice(isnan(a2iChoice)) = -1;  %replace all NaN with -1
    % Area:
    [a2dMap_area, a2dCoord_area] = arcgridread(sFileName_area);
    a2dArea = a2dMap_area;
    % Cell:
    [a2dMap_cell, a2dCoord_cell] = arcgridread(sFileName_cell);
    a2dCelle = a2dMap_cell;
    a2dCelle(isnan(a2dCelle)) = 1; %replace all NaN with 1
    % Dem:
    [a2dMap_dem, a2dCoord_dem] = arcgridread(sFileName_dem);
    a2dDem = a2dMap_dem;
    a2dDem(isnan(a2dDem)) = -99;  %replace all NaN with -99
    % Pointers:
    [a2dMap_pnt, a2dCoord_pnt] = arcgridread(sFileName_pnt);
    a2iPunt = a2dMap_pnt;
    a2iPunt(isnan(a2iPunt)) = 0;   %replace all NaN with 0
    % Lon:
    [a2dMap_lon, a2dCoord_lon] = arcgridread(sFileName_lon);
    Londem = a2dMap_lon;
    Londem(isnan(Londem)) = -99;
    % Lat:
    [a2dMap_lat, a2dCoord_lat] = arcgridread(sFileName_lat);
    Latdem = a2dMap_lat;
    Latdem(isnan(Latdem)) = -99;

    % % important: if the imported "lat" and "lon" have -9999 values you need to 
    % % modify them in order to avoid any -9999. Create the grid manually:
    % [lon, lat] = meshgrid(12.0717-0.005006:0.005006:14.0591 + 0.005006*2,  ...
    %                       42.5264:0.005006:44.0733+0.005006);
    % lat = flip(lat,1);
    % %lat = lat.';
    % Latdem = lat;
    % Londem = lon;

    %% Qindex:
    % COMPUTE Qindex if not available:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % path containing discharge values for each cell of the grid for a time
%     % period >= 10 years:
%     sPathQmaps = [path_preparation_data, '/discharge_2010_2020_NEW'];
%     [iNRows,iNCols]= size(a2iChoice);
%     sDateFrom      = '201101012300';
%     sDateTo        = '202012312300';
%     %define time step of files:
%     dt             =  1;  %days
%     iNoData        = -9999; % valore per dati mancanti nelle mappe Netcdf
%     % define startuing and ending time of simulations:
%     nDateFrom      = datenum(sDateFrom,'yyyymmddHHMM');
%     nDateTo        = datenum(sDateTo,'yyyymmddHHMM');
%     %initialise current time:
%     nNow           = nDateFrom;
%     %initialise matrixes with discharge:
%     a3dMapQ        = zeros([iNRows,iNCols,10]);
%     maxQannual     = zeros([iNRows,iNCols]);
%     % initialise counter:
%     iCountDay      = 0;
%     iCountYear     = 1;
% 
%     while nNow<=nDateTo
%         iCountDay = iCountDay + 1;
%         sDate = datestr(nNow,'yyyymmddHHMM');
%         disp(sDate);
%         % extract month, day, hour from file name:  
%         iYear   = str2double(sDate(1:4));
%         iMonth  = str2double(sDate(5:6));
%         iDay    = str2double(sDate(7:8));
%         iHour   = str2double(sDate(9:10));
%         sPathNow = [sPathQmaps,'/',sDate(1:4),'/',sDate(5:6),'/',sDate(7:8),'/'];
%         try
%             sFileNameMap = ['hmc.output-grid.',sDate,'.nc.gz'];   
%             a2dMap = Continuum_getMap_NC(sPathNow, sFileNameMap,'Discharge');  %istantaneous Discharge (m)
%             a2dMap(a2dMap==iNoData) = NaN;
%             % compute maximum between actual Q(t) and all previous ones of the
%             % year:
%             a3dMapQ(:,:,iCountYear) = max(a3dMapQ(:,:,iCountYear), a2dMap);
%             %a3dMapQ(:,:,iCountYear) = arrayfun(@(x,y) max(x(:),y(:)), a3dMapQ(:,:,iCountYear), a2dMap);
%         catch
%             disp('problem with format of netcdf file!! skip!');
%         end
% 
%         if iMonth==12 & iDay ==31 & iHour == 23
%                 maxQannual(:,:,iCountYear) =  a3dMapQ(:,:,iCountYear);
%                 %display(maxQannual);
%                 display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
%                 iCountDay = 0;
%                 % go to next year:
%                 iCountYear = iCountYear + 1;
%         end   
%         % pass to next time step:
%         nNow = datenum(sDate,'yyyymmddHHMM')+dt;
%     end
% 
%     % Calcolo mappa di portata massima annuale media su deici anni:
%     maxQannual_mean(:,:) = zeros([iNRows,iNCols]);
%     for iCountYear=1:10
%         tmp = maxQannual_mean(:,:);
%         %tmp(OutOfCN) = NaN;
%         maxQannual_mean(:,:) = tmp + maxQannual(:,:,iCountYear);
%     end
%     maxQannual_mean(:,:) = maxQannual_mean(:,:)./10;
%     a2dQindice = maxQannual_mean;
% 
%     % Save obtained Qindex layer for future applications:
%     save([path_preparation_data, '/Qindex_bis_',domain_name,'.mat'], 'a2dQindice');
%     
    % or just load a previously computed .mat file with "a2dQindice":
    load([path_preparation_data, '/Qindex_bis_',domain_name,'.mat'])
    
    
    
elseif strcmp(domain_name,'Entella')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %already available import the Data_*_*Domain.mat file:
    sFileNameOutput =   [sPathGridGeoData, '/Data_LiguriaDomain.mat'];
    load(sFileNameOutput);
    
    
elseif strcmp(domain_name,'Scrivia')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %already available import the Data_*_*Domain.mat file:
    sFileNameOutput =   [sPathGridGeoData, '/Data_LiguriaDomain.mat'];
    load(sFileNameOutput);

    
elseif strcmp(domain_name,'Caccamo-Grazie') | strcmp(domain_name,'Grazie-Foce') | strcmp(domain_name,'Polverina-Caccamo') | strcmp(domain_name,'Chienti')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sFileName_choice = [sPathGridGeoData, '/marche.choice.txt'];
    sFileName_area = [sPathGridGeoData, '/marche.area.txt'];
    sFileName_cell = [sPathGridGeoData, '/marche.areacell.txt'];
    sFileName_dem = [sPathGridGeoData, '/marche.dem.txt'];
    sFileName_lon = [sPathGridGeoData, '/marche.lon.txt'];
    sFileName_lat = [sPathGridGeoData, '/marche.lat.txt'];
    sFileName_pnt = [sPathGridGeoData, '/marche.pnt.txt'];
    % Choice:
    [a2dMap_choice, a2dCoord_choice] = arcgridread(sFileName_choice);
    a2iChoice = a2dMap_choice;
    a2iChoice(isnan(a2iChoice)) = -1;  %replace all NaN with -1
    % Area:
    [a2dMap_area, a2dCoord_area] = arcgridread(sFileName_area);
    a2dArea = a2dMap_area;
    % Cell:
    [a2dMap_cell, a2dCoord_cell] = arcgridread(sFileName_cell);
    a2dCelle = a2dMap_cell;
    a2dCelle(isnan(a2dCelle)) = 1; %replace all NaN with 1
    % Dem:
    [a2dMap_dem, a2dCoord_dem] = arcgridread(sFileName_dem);
    a2dDem = a2dMap_dem;
    a2dDem(isnan(a2dDem)) = -99;  %replace all NaN with -99
    % Pointers:
    [a2dMap_pnt, a2dCoord_pnt] = arcgridread(sFileName_pnt);
    a2iPunt = a2dMap_pnt;
    a2iPunt(isnan(a2iPunt)) = 0;   %replace all NaN with 0
    % Lon:
    [a2dMap_lon, a2dCoord_lon] = arcgridread(sFileName_lon);
    Londem = a2dMap_lon;
    Londem(isnan(Londem)) = -99;
    % Lat:
    [a2dMap_lat, a2dCoord_lat] = arcgridread(sFileName_lat);
    Latdem = a2dMap_lat;
    Latdem(isnan(Latdem)) = -99;
    %Qindex:
    load([path_preparation_data, '/Qindex_bis_Chienti.mat']) 
end













%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Creazione mat file "Aree_finali_Foglia.mat":
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%wm = webmap('World Imagery'); % satellite world map

%% define the domain grid (corners): 
if strcmp(domain_name,'Foglia')
        Lat_min=43.776;
        Lat_max=43.926;
        Lon_min=12.485;
        Lon_max=12.929;

elseif strcmp(domain_name,'Scrivia')
        Lat_min=44.48;
        Lat_max=44.68;
        Lon_min=8.93;
        Lon_max=9.1;

elseif strcmp(domain_name,'Entella')
        Lat_min=44.3;
        Lat_max=44.44;
        Lon_min=9.18;
        Lon_max=9.5;  
%         Lat_min=44.328;
%         Lat_max=44.350;
%         Lon_min=9.36;
%         Lon_max=9.44;
elseif strcmp(domain_name,'EntellaCompleto')
        Lat_min=44.3;
        Lat_max=44.44;
        Lon_min=9.18;
        Lon_max=9.5;

elseif strcmp(domain_name,'Grazie-Foce')
%         Lat_min=43.03;
%         Lat_max=43.32;
%         Lon_min=13.23;
%         Lon_max=13.8;
        Lat_min = 43.18; 
        Lat_max = 43.3;
        Lon_min = 13.25; 
        Lon_max = 13.76;

elseif strcmp(domain_name,'Polverina-Caccamo')
        Lat_min=43.142;
        Lat_max=43.190;
        Lon_min=13.205;
        Lon_max=13.275;

elseif strcmp(domain_name,'Caccamo-Grazie')
        Lat_min=43.085;
        Lat_max=43.148;
        Lon_min=13.109;
        Lon_max=13.207;

elseif strcmp(domain_name,'Chienti')
        Lat_min=43.0845;
        Lat_max=43.301;
        Lon_min=13.103;
        Lon_max=13.7497;
   
end




% define corners indexes:
temp=find(((Londem(1,:)-Lon_min)<0)&(Londem(1,:)-Lon_min)>-100);
indice_y_min=temp(end);clear temp
[x,y] =find(((Londem(:,:)-Lon_max)<0)&(Londem(:,:)-Lon_max)>-100);
indice_y_max=y(end);clear temp      
temp =find(((Latdem(:,1)-Lat_min)>0)&(Latdem(:,1)-Lat_min)>-100);
indice_x_max=temp(end);clear temp
temp=find(((Latdem(:,1)-Lat_max)>0)&(Latdem(:,1)-Lat_max)>-100);
indice_x_min=temp(end);clear temp
% save([path_preparation_data, '/Aree_competenza_',domain_name,'.mat'],'Lat_min','Lat_max',...
%      'Lon_min','Lon_max', 'indice_y_min','indice_y_max','indice_x_max','indice_x_min')




        
       


%% Define the new layers according to the selected domain grid
% clip maps:
Area_dominio    = a2dArea(indice_x_min:indice_x_max, indice_y_min:indice_y_max);
Choice_dominio  = a2iChoice(indice_x_min:indice_x_max,indice_y_min:indice_y_max);
Punt_dominio    = a2iPunt(indice_x_min:indice_x_max,indice_y_min:indice_y_max);
Lat_dominio     = Latdem(indice_x_min:indice_x_max,indice_y_min:indice_y_max);
Lon_dominio     = Londem(indice_x_min:indice_x_max,indice_y_min:indice_y_max);
Qindice_dominio = a2dQindice(indice_x_min:indice_x_max,indice_y_min:indice_y_max);

% plot to verify:
imagesc(Area_dominio);
caxis([2 30]);








%% SECTIONS: select the control sections

if strcmp(modeSelectionSections,'byprompt')
    %% number of sections (USER SELECTION BY PROMPT):
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    quante_sez=input('Quante sezioni?')
    figure 
    imagesc(Area_dominio)
    %imagesc(a2dArea)
    caxis([2 20])
    % scelgo dove sono le sezioni
    [pointclick_x,pointclick_y] = ginput(quante_sez);
    close all
    % costruisco file delle sezioni
    sezioni_indici_relativi=[round(pointclick_y),round(pointclick_x)];
    % nomi sezioni
    for ind=1:quante_sez
        nomi_sezioni{ind}='NoName';
    end
    sezioni_indici_relativi_corr2 = sezioni_indici_relativi;
    
    
    
    
    
    
  
%% Instead, here chosen sections are hard-coded:
elseif strcmp(modeSelectionSections,'hardcoded')
    
       
    if strcmp(domain_name,'Foglia')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % define the number of control sections:
        quante_sez=3;
        Ymax= length(Latdem(:,1));
        Xmax= length(Latdem(1,:));
        % define indexes of the selected sections [Y X]:
        sezioni_indici_relativi=[106 155;...  % Bronzo   
                                 92 241;... % Montecchio
                                 72 286];    % OUT   
%         sezioni_indici_relativi =[41 7;...  % Bronzo        % GRID MARCHE   
%                                   27 93;... % Montecchio
%                                   8 138];   % OUT
        % define the names of each selected section:                     
        nomi_sezioni{1}='Bronzo';
        nomi_sezioni{2}='Montecchio';
        nomi_sezioni{3}='Foglia3';
        % define name of catchment and section:
        bacino_sezione{1} = 'Foglia_Bronzo'; 
        bacino_sezione{2} = 'Foglia_Montecchio'; 
        bacino_sezione{3} = 'Foglia_Foce'; 

        % verifica aree e calcolo Qindex)
        for i=1:size(sezioni_indici_relativi,1)
            AreaBas(i)=a2dArea(sezioni_indici_relativi(i,1),sezioni_indici_relativi(i,2));
            a1dQindex(i)=a2dQindice(sezioni_indici_relativi(i,1),sezioni_indici_relativi(i,2));
        end
        %Ricavo indici relativi al ritaglio
        sezioni_indici_relativi_corr2(:,1)=sezioni_indici_relativi(:,1)-indice_x_min+1;
        sezioni_indici_relativi_corr2(:,2)=sezioni_indici_relativi(:,2)-indice_y_min+1;
    
%         quante_sez=3;
%         % define indexes of the selected sections [Y X]:
%         sezioni_indici_relativi=[25 6;...   % Bronzo        % GRID FLOODPROOF   
%                                  17 57;...  % Montecchio
%                                  5 84];     % OUT
% 
%         % define the names of each selected section:                     
%         nomi_sezioni{1}='Bronzo';
%         nomi_sezioni{2}='Montecchio';
%         nomi_sezioni{3}='Foglia3';
%         % define name of catchment and section:
%         bacino_sezione{1} = 'Foglia_Bronerence system code in the inputs!  EPSG_domain  
% salvo geotiff e .mat aree competenza (controllare):
geotiffwrite([path_preparation_data,'/prova_area', domain_name,'.tif'], double(AreeCompetenza),R, 'CoordRefSysCode',['EPSG:',EPSG_domain])



% Save obtained information into mat file (this file is used as static input by Flomart application):
save([sPathOutputInfoMat, '/info_',domain_name,'.mat'], ...
    'LonLL','LonUR','LatLL','LatUR',...
    'AreeCompetenza','mappa_aree_allargata','mappa_aree','Lat_dominio_UTM','Lon_dominio_UTM','a1dQindex','nomi_sezioni_sort', ...
     'indici_sort','bacino_sezione_sort', 'EPSG_domain', 'drainage_area_section_km2');
  







%%zo'; 
%         bacino_sezione{2} = 'Foglia_Montecchio'; 
%         bacino_sezione{3} = 'Foglia_Foce'; 
%         %% INDEXES CORRECTION (added by M.Darienzo):
%         % Ricavo indici corretti
%         sezioni_indici_relativi_corr=sezioni_indici_relativi;
%         sezioni_indici_relativi_corr2=sezioni_indici_relativi_corr;
%         % verifica aree e calcolo Qindex)
%         for i=1:size(sezioni_indici_relativi_corr,1)
%             AreaBas(i)=Area_dominio(sezioni_indici_relativi_corr(i,1),sezioni_indici_relativi_corr(i,2));
%             a1dQindex(i)=Qindice_dominio(sezioni_indici_relativi_corr(i,1),sezioni_indici_relativi_corr(i,2));
%         end    

        
               
        
        
    elseif strcmp(domain_name,'Scrivia')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% define the number of control sections:
        quante_sez=4;
        Ymax=442;
        Xmax=1052;
        % define indexes of the selected sections [Y X]:
        sezioni_indici_relativi=[388	576;...  % Busalla   
                                 365	595;... % Casella
                                 418	582; ... %Isola C
                                 356	615; ... %Mont
                                 %415	590; ...%Vobbia
                                 ];    % OUT 
        % define the names of each selected section:                     
        nomi_sezioni{1}='SCRBUS';
        nomi_sezioni{2}='SCRI01';
        nomi_sezioni{3}='SCRISO';
        nomi_sezioni{4}='SCRMON';
        %nomi_sezioni{5}='SCRISO';
        % AGGIUNTO FRANCESCO
        bacino_sezione{1}='Scrivia_Busalla';
        bacino_sezione{2}='Scrivia_Casella';
        bacino_sezione{3}='Scrivia_IsolaDelCantone';
        bacino_sezione{4}='Scrivia_Montoggio';
                
        % INDEXES CORRECTION (added by F. Silvestro):
        % Ricavo indici corretti
        sezioni_indici_relativi_corr(:,1)=Ymax-sezioni_indici_relativi(:,1)+1;
        sezioni_indici_relativi_corr(:,2)=sezioni_indici_relativi(:,2)+0;

        % verifica aree e calcolo Qindex)
        for i=1:size(sezioni_indici_relativi_corr,1)
            AreaBas(i)=a2dArea(sezioni_indici_relativi_corr(i,1),sezioni_indici_relativi_corr(i,2));
            a1dQindex(i)=a2dQindice(sezioni_indici_relativi_corr(i,1),sezioni_indici_relativi_corr(i,2));
        end

        %Ricavo indici relativi al ritaglio
        sezioni_indici_relativi_corr2(:,1)=sezioni_indici_relativi_corr(:,1)-indice_x_min+1;
        sezioni_indici_relativi_corr2(:,2)=sezioni_indici_relativi_corr(:,2)-indice_y_min+1;
    
      
        
        
    elseif strcmp(domain_name,'Entella')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% define the number of control sections:
        quante_sez=6;
        Ymax=442;
        Xmax=1052;
        % define indexes of the selected sections [Y X]:
        sezioni_indici_relativi=[271	745;...  % Caminata   
                                 275	723;... % Carasco
                                 281    721; ... %S. Martino
                                 268	723; ... %Panesi
                                 292	732; ...%Vignolo
                                 257    717];    % OUT 
        % define the names of each selected section:                     
        nomi_sezioni{1}='CAMINA';
        nomi_sezioni{2}='CARASC';
        nomi_sezioni{3}='LAVA01';
        nomi_sezioni{4}='PANESI';
        nomi_sezioni{5}='STIVIG';
        nomi_sezioni{6}='OUT';
        
        % AGGIUNTO FRANCESCO
        bacino_sezione{1}='Graveglia_Caminata';
        bacino_sezione{2}='Lavagna_Carasco';
        bacino_sezione{3}='Lavagna_SMartino';
        bacino_sezione{4}='Entella_Panesi';
        bacino_sezione{5}='Sturla_Vignolo';    
        bacino_sezione{6}='Entella_Foce';

        % INDEXES CORRECTION (added by F. Silvestro):
        % Ricavo indici corretti
        sezioni_indici_relativi_corr(:,1)=Ymax-sezioni_indici_relativi(:,1)+1;
        sezioni_indici_relativi_corr(:,2)=sezioni_indici_relativi(:,2)+0;
        % verifica aree e calcolo Qindex)
        for i=1:size(sezioni_indici_relativi_corr,1)
            AreaBas(i)=a2dArea(sezioni_indici_relativi_corr(i,1),sezioni_indici_relativi_corr(i,2));
            a1dQindex(i)=a2dQindice(sezioni_indici_relativi_corr(i,1),sezioni_indici_relativi_corr(i,2));
        end
        %Ricavo indici relativi al ritaglio
        sezioni_indici_relativi_corr2(:,1)=sezioni_indici_relativi_corr(:,1)-indice_x_min+1;
        sezioni_indici_relativi_corr2(:,2)=sezioni_indici_relativi_corr(:,2)-indice_y_min+1;
        
%         % add Foce section:
%         quante_sez=6;
%         bacino_sezione{6}='Entella_Foce';
%         nomi_sezioni{6}='OUT';
%         sezioni_indici_relativi_corr2(6,:) = [65 57]; 
        

        
     elseif strcmp(domain_name,'EntellaCompleto')
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% define the number of control sections:
        quante_sez=5;
        Ymax= length(Latdem(:,sezioni_indici_relativi_corr21));
        Xmax= length(Latdem(1,:));
        % define indexes of the selected sections [Y X]:
        sezioni_indici_relativi=[271	745;...  % Caminata   
                                 275	723;... % Carasco
                                 281    721; ... %S. Martino
                                 268	723; ... %Panesi
                                 292	732; ...%Vignolo
                                 ];    % OUT 
        % define the names of each selected section:                     
        nomi_sezioni{1}='CAMINA';
        nomi_sezioni{2}='CARASC';
        nomi_sezioni{3}='LAVA01';
        nomi_sezioni{4}='PANESI';
        nomi_sezioni{5}='STIVIG';
        
        
        % AGGIUNTO FRANCESCO
        bacino_sezione{1}='Graveglia_Caminata';
        bacino_sezione{2}='Lavagna_Carasco';
        bacino_sezione{3}='Lavagna_SMartino';
        bacino_sezione{4}='Entella_Panesi';
        bacino_sezione{5}='Sturla_Vignolo';    
  

        % INDEXES CORRECTION (added by F. Silvestro):
        % Ricavo indici corretti
        sezioni_indici_relativi_corr(:,1)=Ymax-sezioni_indici_relativi(:,1)+1;
        sezioni_indici_relativi_corr(:,2)=sezioni_indici_relativi(:,2)+0;

        % verifica aree e calcolo Qindex)
        for i=1:size(sezioni_indici_relativi_corr,1)
            AreaBas(i)=a2dArea(sezioni_indici_relativi_corr(i,1),sezioni_indici_relativi_corr(i,2));
            a1dQindex(i)=a2dQindice(sezioni_indici_relativi_corr(i,1),sezioni_indici_relativi_corr(i,2));
        end

        %Ricavo indici relativi al ritaglio
        sezioni_indici_relativi_corr2(:,1)=sezioni_indici_relativi_corr(:,1)-indice_x_min+1;
        sezioni_indici_relativi_corr2(:,2)=sezioni_indici_relativi_corr(:,2)-indice_y_min+1;
        
        % sezioni_indici_relativi_corr2 = sezioni_indici_relativi;
        
        
        
    elseif strcmp(domain_name,'Grazie-Foce')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % define the number of control sections:
%         quante_sez=5;
%         Ymax= length(Latdem(:,1));
%         Xmax= length(Latdem(1,:));
%         % define indexes of the selected sections [Y X]:
%         sezioni_indici_relativi=[299 454;...  % Fiastra
%                                  297 444;...  % Chienti1  
%                                  288 515;...  % Chienti2
%                                  296 543;...  % EteMorto
%                                  277 565];    % OUT   
%         % define the names of each selected section:
%         nomi_sezioni{1}='Fiastra';
%         nomi_sezioni{2}='Chienti1';
%         nomi_sezioni{3}='Chienti2';
%         nomi_sezioni{4}='EteMorto';
%         nomi_sezioni{5}='FoceChienti';
%         % define name of catchment and section:
% %         bacino_sezione{1} = 'Chienti_Fiastra'; 
% %         bacino_sezione{2} = 'Chienti_Chienti1'; 
% %         bacino_sezione{3} = 'Chienti_Chienti2'; 
% %         bacino_sezione{4} = 'Chienti_EteMorto'; 
% %         bacino_sezione{5} = 'Chienti_Foce'; 
%         bacino_sezione{1} = 'Fiastra_Chienti'; 
%         bacino_sezione{2} = 'Chienti1_Chienti'; 
%         bacino_sezione{3} = 'Chienti2_Chienti'; 
%         bacino_sezione{4} = 'EteMorto_EteMorto'; 
%         bacino_sezione{5} = 'FoceChienti_Chienti'; 
        % define the number of control sections:
        quante_sez=4;
        Ymax= length(Latdem(:,1));
        Xmax= length(Latdem(1,:));
        % define indexes of the selected sections [Y X]:
        sezioni_indici_relativi=[299 454;...  % Fiastra
                                 297 444;...  % Chienti1  
                                 288 515;...  % Chienti2
                                 277 565];    % OUT   
        % define the names of each selected section:
        nomi_sezioni{1}='Fiastra';
        nomi_sezioni{2}='Chienti1';
        nomi_sezioni{3}='Chienti2';
        nomi_sezioni{4}='FoceChienti';
        % define name of catchment and section:
        bacino_sezione{1} = 'Fiastra_Chienti'; 
        bacino_sezione{2} = 'Chienti1_Chienti'; 
        bacino_sezione{3} = 'Chienti2_Chienti'; 
        bacino_sezione{4} = 'FoceChienti_Chienti';         
        

        % verifica aree e calcolo Qindex)
        for i=1:size(sezioni_indici_relativi,1)
            AreaBas(i)=a2dArea(sezioni_indici_relativi(i,1),sezioni_indici_relativi(i,2));
            a1dQindex(i)=a2dQindice(sezioni_indici_relativi(i,1),sezioni_indici_relativi(i,2));
        end
        %Ricavo indici relativi al ritaglio
        sezioni_indici_relativi_corr2(:,1)=sezioni_indici_relativi(:,1)-indice_x_min+1;
        sezioni_indici_relativi_corr2(:,2)=sezioni_indici_relativi(:,2)-indice_y_min+1;
   
       
        
        
        
        
    elseif strcmp(domain_name,'Polverina-Caccamo')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        quante_sez=1;
        sezioni_indici_relativi_corr2 = [3 32];
        nomi_sezioni{1}='DigaCaccamo';
        % define name of catchment and section:
        bacino_sezione{1} = 'DigaCaccamo_Chienti'; 

        % verifica aree e calcolo Qindex)
        for i=1:size(sezioni_indici_relativi_corr2,1)
            AreaBas(i)  =Area_dominio(sezioni_indici_relativi_corr2(i,1),sezioni_indici_relativi_corr2(i,2));
            a1dQindex(i)=a2dQindice(sezioni_indici_relativi_corr2(i,1),sezioni_indici_relativi_corr2(i,2));
        end
     
        
    elseif strcmp(domain_name,'Caccamo-Grazie')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        quante_sez=1;
        sezioni_indici_relativi_corr2 = [3 22];
        nomi_sezioni{1}='DigaGrazie';
        % define name of catchment and section:
        bacino_sezione{1} = 'DigaGrazie_Chienti'; 

        % verifica aree e calcolo Qindex)
        for i=1:size(sezioni_indici_relativi_corr2,1)
            AreaBas(i)  =Area_dominio(sezioni_indici_relativi_corr2(i,1),sezioni_indici_relativi_corr2(i,2));
            a1dQindex(i)=a2dQindice(sezioni_indici_relativi_corr2(i,1),sezioni_indici_relativi_corr2(i,2));
        end
        
        
        
        
       

    elseif strcmp(domain_name,'Chienti')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % define the number of control sections:
        quante_sez=3;
        Ymax= length(Latdem(:,1));
        Xmax= length(Latdem(1,:));
        % define indexes of the selected sections [Y X]:
        sezioni_indici_relativi=[%299 454;...  % Fiastra
                                 297 444;...  % Chienti1  
                                 288 515;...  % Chienti2
                                 277 565];    % OUT   
        % define the names of each selected section:
        %nomi_sezioni{1}='Fiastra';
        nomi_sezioni{1}='Chienti1';
        nomi_sezioni{2}='Chienti2';
        nomi_sezioni{3}='FoceChienti';
        % define name of catchment and section:
        %bacino_sezione{1} = 'Fiastra_Chienti'; 
        bacino_sezione{1} = 'Chienti1_Chienti'; 
        bacino_sezione{2} = 'Chienti2_Chienti'; 
        bacino_sezione{3} = 'FoceChienti_Chienti';         

        % verifica aree e calcolo Qindex)
        for i=1:size(sezioni_indici_relativi,1)
            AreaBas(i)=a2dArea(sezioni_indici_relativi(i,1),sezioni_indici_relativi(i,2));
            a1dQindex(i)=a2dQindice(sezioni_indici_relativi(i,1),sezioni_indici_relativi(i,2));
        end
        %Ricavo indici relativi al ritaglio
        sezioni_indici_relativi_corr2(:,1)=sezioni_indici_relativi(:,1)-indice_x_min+1;
        sezioni_indici_relativi_corr2(:,2)=sezioni_indici_relativi(:,2)-indice_y_min+1;
        
        
    end
    
else
    error('ERROR: user-defined mode of sections selection not available');     
end  
      



% % plot to verify sections on starting grid:
% figure 
% imagesc(a2dArea)
% caxis([2 50])
% hold on
% for indicew=1:length(sezioni_indici_relativi)
%     plot(sezioni_indici_relativi(indicew,2), sezioni_indici_relativi(indicew,1),'or','markersize',13, 'LineWidth',5)
% end


% plot to verify sections on clipped grid:
figure 
imagesc(Area_dominio)
caxis([2 50])
hold on
for indicew=1:length(sezioni_indici_relativi_corr2)
    plot(sezioni_indici_relativi_corr2(indicew,2), sezioni_indici_relativi_corr2(indicew,1),'or','markersize',13, 'LineWidth',5)
end



     


     

%% AREE DRENATE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the drainage area for each selected section:
cd(path_code)
aree = aree_drenate(Punt_dominio, sezioni_indici_relativi_corr2);

% Delete superposed areas:
mappa_aree=zeros(size(Punt_dominio));
L=cellfun(@length,aree);
[ordinati,indici_sort]=sort(L);
for i=length(L):-1:1
        valori=unique(mappa_aree(indici_sort(i)));
        mappa_aree(aree{indici_sort(i)})=i;
        nomi_sezioni_sort{i}=nomi_sezioni{indici_sort(i)};     %Aggiunto FRANCESCO
        bacino_sezione_sort{i}=bacino_sezione{indici_sort(i)}; %Aggiunto FRANCESCO
        drainage_area_sort{i}=AreaBas(indici_sort(i));         %Aggiunto MATTEO
end

    

% Verify areas with a PLOT:
figure
imagesc(mappa_aree)
hold on
[canalix,canaliy]=find(Choice_dominio==1);
for indicew=1:length(canalix)
    plot(canaliy(indicew),canalix(indicew),'.k','markersize',6)
end
for indicew=1:length(sezioni_indici_relativi_corr2)
    plot(sezioni_indici_relativi_corr2(indicew,2),sezioni_indici_relativi_corr2(indicew,1),'or','markersize',6)
end
    




%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% modification in order to increase/adjust the dominium boundaries
[righe,colonne] = size(mappa_aree);
mappa_aree_allargata = mappa_aree;
for indice=1:colonne
        temp=find(mappa_aree(:,indice)>0,1,'first');
        if isempty(temp)
        else
            if temp>1
                mappa_aree_allargata(temp-1,indice)=mappa_aree(temp,indice);
            end
        end
        clear temp
        temp=find(mappa_aree(:,indice)>0,1,'last');
        if isempty(temp)
        else
            if temp<righe
                mappa_aree_allargata(temp+1,indice)=mappa_aree(temp,indice);
            end
        end
        clear temp
end
for indice=1:righe
        temp=find(mappa_aree(indice,:)>0,1,'first');
        if isempty(temp)
        else
            if temp>1
                mappa_aree_allargata(indice,temp-1)=mappa_aree(indice,temp);
            end
        end
        clear temp
        temp=find(mappa_aree(indice,:)>0,1,'last');
        if isempty(temp)
        else
            if temp<colonne
                mappa_aree_allargata(indice,temp+1)=mappa_aree(indice,temp);
            end
        end
        clear temp
end

% verify PLOT:
figure
imagesc(mappa_aree_allargata)
hold on
[canalix,canaliy]=find(Choice_dominio==1);
for indicew=1:length(canalix)
        plot(canaliy(indicew),canalix(indicew),'.k','markersize',6)
end
for indicew=1:length(sezioni_indici_relativi_corr2)
        plot(sezioni_indici_relativi_corr2(indicew,2),sezioni_indici_relativi_corr2(indicew,1),'or','markersize',6)
end





    
% IF NECESSARY:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% manually modify "mappa_aree_allargata" in order to include the river 
% delta or other zones known to be flooded that are not included in the domain:
if strcmp(domain_name,'Foglia')
%       mappa_aree_allargata(2:6,83:89) = quante_sez;        %foce  GRID FP
%        mappa_aree_allargata(5:18,161:171) = quante_sez;    %foce GRID MARCHE
        mappa_aree_allargata(3:10,135:145) = quante_sez;     %foce GRID MARCHE
          
elseif strcmp(domain_name,'Graveglia')
        disp('No allargamento !');
        
        
elseif strcmp(domain_name,'Scrivia')
        disp('No allargamento !');
        
elseif strcmp(domain_name,'Entella')
        mappa_aree_allargata(54:62,59:66) = quante_sez; %foce
        mappa_aree_allargata(60:66,49:63) = quante_sez; %foce
        % mappa_aree_allargata(2:4,81:83) = quante_sez;  

elseif strcmp(domain_name,'EntellaCompleto')
        mappa_aree_allargata(54:62,59:66) = quante_sez; %foce
        mappa_aree_allargata(60:66,49:63) = quante_sez; %foce
        % mappa_aree_allargata(2:4,81:83) = quante_sez;  
        
elseif strcmp(domain_name,'Grazie-Foce')
        %mappa_aree_allargata(8:11,168:172) = quante_sez; %foce
        mappa_aree_allargata(3:8,160:170) = quante_sez; %foce

elseif strcmp(domain_name,'Polverina-Caccamo')
        disp('No allargamento !');

elseif strcmp(domain_name,'Caccamo-Grazie')
        disp('No allargamento !');

elseif strcmp(domain_name,'Chienti')
        %mappa_aree_allargata(8:11,168:172) = quante_sez; %foce
        mappa_aree_allargata(2:9,210:218) = quante_sez; %foce
end   

    
% verify PLOT:
figure
imagesc(mappa_aree_allargata)
hold on
[canalix,canaliy]=find(Choice_dominio==1);
for indicew=1:length(canalix)
        plot(canaliy(indicew),canalix(indicew),'.k','markersize',6)
end
for indicew=1:length(sezioni_indici_relativi_corr2)
        plot(sezioni_indici_relativi_corr2(indicew,2),sezioni_indici_relativi_corr2(indicew,1),'or','markersize',6)
end
%%%%%%%%%% .mat aree competenza (controllare):
save([path_preparation_data,'/mappa_aree_allargata_', domain_name,'.mat'], 'mappa_aree_allargata')
    
    



    
   
    
%%
% recover info of metric grid UTM from the hazard maps:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
name_file_hazard_read=[sPathHazardData, '/', domain_name, '/', name_hazardmaps, sprintf('%03.0f',TR_min), type_hazardmaps];
% name_file_hazard_read=[sPathHazardData, '/', domain_name, '/', name_hazardmaps, sprintf('%03.0f',13), type_hazardmaps];

%name_file_read=[sPathHazardData, '/', domain_name,'/', name_hazardmaps, num2str(TR_max),nome_hazmaps1];
%name_file_hazard_read=[path_preparation_data, '/telemac_data/', domain_name,'_WD_max_Q100.tif'];
[A, R]=geotiffread(name_file_hazard_read);



% % ONLY IF NECESSARY:
% % it may be possible that the available hazard maps are located into 3
% % separated abacus (for different branch of the river).
% % Here, you can merge the hazard maps in one unique map for the whole basin:
% name_file_hazard_read2 =[sPathHazardData, '/Polverina-Caccamo/Polverina-Caccamo_WD_max_T', sprintf('%03.0f',TR_min), type_hazardmaps];
% [A2, R2] = geotiffread(name_file_hazard_read2);
% name_file_hazard_read3 =[sPathHazardData, '/Caccamo-Grazie/Caccamo-Grazie_WD_max_T', sprintf('%03.0f',TR_min), type_hazardmaps];
% [A3, R3]= geotiffread(name_file_hazard_read3);
% [Z_out, R_out] = imfuse(A2,R2,A3,R3);  % requires image processing toolbox







%%
% Rigriglio su griglia metrica UTM:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function 'meshgrid(x,y)' transforms the domain specified by vectors x and y into arrays X and Y, 
% which can be used to evaluate functions of two variables and three-dimensional mesh/surface plots. 
% The rows of the output array X are copies of the vector x;
% The columns of the output array Y are copies of the vector y.
[new_x,new_y]= meshgrid(R.XWorldLimits(1):R.CellExtentInWorldX:R.XWorldLimits(2) - R.CellExtentInWorldX, ...
                        R.YWorldLimits(1):R.CellExtentInWorldY:R.YWorldLimits(2) - R.CellExtentInWorldY);
%AGGIUNTO Francesco S.
LonLL=R.XWorldLimits(1);
LonUR=R.XWorldLimits(2);
LatLL=R.YWorldLimits(1);
LatUR=R.YWorldLimits(2);

% convert coord lat/lon to coordinates utm32: 
[coord_left, coord_bottom] = latlon2utm(Lat_min, Lon_min);
[coord_right, coord_top]   = latlon2utm(Lat_max, Lon_max);



% get the resolution of the domain.mat in degrees and meters:
% res_Lon = ceil((coord_right - coord_left)/size(mappa_aree_allargata,2));
% res_Lat = ceil((coord_top - coord_bottom)/size(mappa_aree_allargata,1));
res_Lon = ((coord_right - coord_left)/(size(mappa_aree_allargata,2)-1));
res_Lat = ((coord_top - coord_bottom)/(size(mappa_aree_allargata,1)-1));
[Lon_dominio_UTM, Lat_dominio_UTM] = meshgrid(coord_left:res_Lon:coord_right, ...
                                              coord_bottom:res_Lat:coord_top);

                                          
% Compute drainage area in km2 for each selected section:                                  
if strcmp(domain_name,'Foglia') | strcmp(domain_name,'Grazie-Foce') | strcmp(domain_name,'Caccamo-Grazie') | strcmp(domain_name,'Polverina-Caccamo')| strcmp(domain_name,'Chienti') 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for  i=1:length(drainage_area_sort)
        drainage_area_section_km2{i} = drainage_area_sort{i}/1000000*res_Lon*res_Lat; 
    end
    
elseif strcmp(domain_name,'Entella') | strcmp(domain_name,'EntellaCompleto') | strcmp(domain_name,'Scrivia')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for  i=1:length(drainage_area_sort)
        drainage_area_section_km2{i} = drainage_area_sort{i};   %/1000000*res_Lon*res_Lat; 
    end
end





% Interpolate the sample grid data using the griddata() and the 'nearest' method:
% griddata(x, y, v, xq, yq) fits a surface of the form v = f(x,y) to the scattered 
% data in the vectors (x,y,v). 
% The griddata function interpolates the surface at the query points specified by (xq,yq) 
% and returns the interpolated values. The surface always passes through the data points defined by x and y.

% check if Lat_dominio_UTM(1,1) < Lat_dominio_UTM(length(Lat_dominio_UTM)), then flip:
if new_y(1,1) > new_y(length(new_y))
    if new_x(1,1) > new_x(length(new_x))
        display('Flipping both y and x is needed! Please, wait ...');
        % flipping both y and x:
        AreeCompetenza = griddata(Lat_dominio_UTM, Lon_dominio_UTM, mappa_aree_allargata, ...
                                  flipud(new_y), flipud(new_x), 'nearest');
    else
        display('Flipping y is needed! Please, wait ...');
        % flipping y only:
        AreeCompetenza = griddata(Lat_dominio_UTM, Lon_dominio_UTM, mappa_aree_allargata, ...
                                  flipud(new_y), new_x, 'nearest');
    end
        
else
    
    if new_x(1,1) > new_x(length(new_x))
        %  flipping x only:     
        display('Flipping x is needed! Please, wait ...');
        AreeCompetenza = griddata(Lat_dominio_UTM, Lon_dominio_UTM, mappa_aree_allargata, ...
                                  new_y, flipud(new_x), 'nearest');
    else
        % No flipping needed:     
        display('No flipping needed, Please, wait ...');
        AreeCompetenza = griddata(Lon_dominio_UTM, Lat_dominio_UTM, mappa_aree_allargata, ...
                                  new_x, new_y, 'nearest');
    end
end
figure; imagesc(AreeCompetenza);



% AreeCompetenzaBIG = AreeCompetenza;
% [Ny,Nx,~] = size(AreeCompetenzaBIG);
% idx1 = 1:4:Ny;
% idx2 = 1:4:Nx;
% B = AreeCompetenzaBIG(idx1,idx2,:);
% AreeCompetenza = B;

%% ! please, verify the geo reference system code in the inputs!  EPSG_domain  
% salvo geotiff e .mat aree competenza (controllare):
geotiffwrite([path_preparation_data,'/prova_area', domain_name,'.tif'], double(AreeCompetenza),R, 'CoordRefSysCode',['EPSG:',EPSG_domain])



% Save obtained information into mat file (this file is used as static input by Flomart application):
save([sPathOutputInfoMat, '/info_',domain_name,'.mat'], ...
    'LonLL','LonUR','LatLL','LatUR',...
    'AreeCompetenza','mappa_aree_allargata','mappa_aree','Lat_dominio_UTM','Lon_dominio_UTM','a1dQindex','nomi_sezioni_sort', ...
     'indici_sort','bacino_sezione_sort', 'EPSG_domain', 'drainage_area_section_km2');
  







%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    EXAMPLE OF OPERATIONAL FLOOD MAP GENERATION:                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% insert values of peak discharge for each selected section:
discharge_values_sections =[800, 2000, 3050]; %Foglia
discharge_values_sections =[300, 200, 300, 400, 500, 1000]; %Entella
discharge_values_sections =[200, 300, 400, 500]; %Scrivia
discharge_values_sections =[200, 300, 400, 500, 1000]; %Grazie-Foce


examp_map = example_generation_flood_map(domain_name, sPathHazardData, name_hazardmaps, path_preparation_data, ...
                                         discharge_values_sections, quante_sez, Qindice_dominio, ...
                                         sezioni_indici_relativi_corr2, nomi_sezioni, ...
                                         AreeCompetenza, R, TR_max, TR_min, EPSG_domain);



