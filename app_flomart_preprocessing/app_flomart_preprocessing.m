clc; clear;

%%
%%%%%%%%%%%%%%
%   INPUTS   %
%%%%%%%%%%%%%%
% please, if a new case study need to be analysed, then insert the name of the 
% new domain in the list below. Then execute next command to select the
% case study from prompt pop-up:
list_domains = {'Foglia', 'Entella', 'Scrivia', 'Chienti'}
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
    name_hazardmaps = ['compressed/Chienti_WD_max_T'];
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
    sFileName_AreaDrenata = [path_preparation_data, '/Mappe_Regionalizzazione_Q_v2/Mappa_aree_drenate.asc'];
    
    
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
    % load a previously computed .mat file with "a2dQindice":
    % load([path_preparation_data, '/Qindex_bis_',domain_name,'.mat'])

    %% Regionalizzazione Marche:
    % Choice:
    [a2dMap_areadrenata, a2dCoord_areadrenata] = arcgridread(sFileName_AreaDrenata);
    a2iAreaDrenata = a2dMap_areadrenata;
    % a2iAreaDrenata(isnan(a2iAreaDrenata)) = -1;  %replace all NaN with -1
    a2dQindice = 1.6119*a2iAreaDrenata(:,:).^0.9735;
    
    % plot to verify:
    imagesc(a2dMap_area);
    caxis([2 30]);

    
    
elseif strcmp(domain_name,'Entella') | strcmp(domain_name,'Scrivia')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   if already available import the Data_*_*Domain.mat file:
%   sFileNameOutput =   [sPathGridGeoData, '/Data_LiguriaDomain.mat'];
%   load(sFileNameOutput);

%   otherwhise upload ascii data:
    sFileName_choice = [sPathGridGeoData, '/LiguriaDomain.choice.txt'];
    sFileName_area = [sPathGridGeoData, '/LiguriaDomain.area.txt'];
    sFileName_pnt = [sPathGridGeoData, '/LiguriaDomain.pnt.txt'];
    % Choice:
    [a2dMap_choice, a2dCoord_choice] = arcgridread(sFileName_choice);
    a2iChoice = a2dMap_choice;
    a2iChoice(isnan(a2iChoice)) = -1;  %replace all NaN with -1
    % Area:
    [a2dMap_area, a2dCoord_area] = arcgridread(sFileName_area);
    a2dArea = a2dMap_area;
    % Pointers:
    [a2dMap_pnt, a2dCoord_pnt] = arcgridread(sFileName_pnt);
    a2iPunt = a2dMap_pnt;
    a2iPunt(isnan(a2iPunt)) = 0;   %replace all NaN with 0


    %upload Qindex, Latdem, Londem from .mat file:
    load([sPathGridGeoData, '/Data_LiguriaDomain.mat']);
    a2dQindex  = a2dQindice;
    Latdem = Latdem; 
    Londem = Londem; 
    
  

    
elseif strcmp(domain_name,'Chienti')
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define the new domain grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% define the new domain grid (corners): 
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





%% extra tool (Create new sections indexes)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % transform lat/lon coordinates into matricial continuum indixes:
% tmp = a2dMap_choice;
% masknan = nan(size(tmp));
% masknan(tmp==1) = 333;
% c = readtable([path_preparation_data, '/nuove_sezioni_flomart_Chienti_Foglia.csv']);
% x_points = c.Lon;
% y_points = c.Lat;
% namfascii = sFileName_choice;
% [rr, cc] = Get_MatrixCoordinates_FromPoints_NoNaN(x_points, y_points, namfascii, masknan);
% %Quindi, salvare le coordinate matrice [rr, cc] nel file “rigacolonna.txt”.
% fid = fopen([path_preparation_data, '/rigacolonna.txt'], 'wt');
% for i = 1:length(rr)
%     fprintf(fid,'%.0f', rr(i));
%     fprintf(fid,'%s', " ");
%     fprintf(fid,'%.0f', cc(i));
%     fprintf(fid,'\n');
% end
% fclose(fid);




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
        quante_sez=6;
        Ymax= length(Latdem(:,1));
        Xmax= length(Latdem(1,:));
        % define indexes of the selected sections [Y X]:
        sezioni_indici_relativi=[106 155;...  % Bronzo 
                                 108 190;...  % FlomartTorreCotogna
                                 108 207;...  % FlomartLaBadia
                                 93 228;...   % FlomartCasellaMontecchio
                                 92 241;...   % Montecchio
                                 72 286];     % Foglia3 (Foce)   
        % define the names of each selected section:                     
        nomi_sezioni{1}='Bronzo_Foglia'; 
        nomi_sezioni{2}='FlomartTorreCotogna_Foglia';
        nomi_sezioni{3}='FlomartLaBadia_Foglia';
        nomi_sezioni{4}='FlomartCasellaMontecchio_Foglia';
        nomi_sezioni{5}='Montecchio_Foglia';
        nomi_sezioni{6}='Foglia3_Foglia';
        
        % define name of catchment and section:
        bacino_sezione{1} = 'Foglia_Bronzo';
        bacino_sezione{2} = 'Foglia_FlomartTorreCotogna';
        bacino_sezione{3} = 'Foglia_FlomartLaBadia';
        bacino_sezione{4} = 'Foglia_FlomartCasellaMontecchio';
        bacino_sezione{5} = 'Foglia_Montecchio'; 
        bacino_sezione{6} = 'Foglia_Foce'; 
        
        % verifica aree e calcolo Qindex)
        for i=1:size(sezioni_indici_relativi,1)
            AreaBas(i)=a2dArea(sezioni_indici_relativi(i,1),sezioni_indici_relativi(i,2));
            a1dQindex(i)=a2dQindice(sezioni_indici_relativi(i,1),sezioni_indici_relativi(i,2));
        end
        %Ricavo indici relativi al ritaglio
        sezioni_indici_relativi_corr2(:,1)=sezioni_indici_relativi(:,1)-indice_x_min+1;
        sezioni_indici_relativi_corr2(:,2)=sezioni_indici_relativi(:,2)-indice_y_min+1;

                  
        
        
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
        
        % AGGIUNTO MATTEO:
        section_name_part1{1}='Graveglia';
        section_name_part1{2}='Lavagna';
        section_name_part1{3}='Lavagna';
        section_name_part1{4}='Entella';
        section_name_part1{5}='Sturla';    
        section_name_part1{6}='Entella';
        
        section_name_part2{1}='Caminata';
        section_name_part2{2}='Carasco';
        section_name_part2{3}='SMartino';
        section_name_part2{4}='Panesi';
        section_name_part2{5}='Vignolo';    
        section_name_part2{6}='Foce';
        
        
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
        
        

    elseif strcmp(domain_name,'Chienti')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % define the number of control sections:
        quante_sez=11;
        Ymax= length(Latdem(:,1));
        Xmax= length(Latdem(1,:));
        % define indexes of the selected sections [Y X]:
        sezioni_indici_relativi=[337 362;...  % FlomartPolverina
                                 326 387;...  % BorgianoDiga
                                 % 319 402;...% FlomartBorgiano (REMOVED)
                                 312 409;...  % LeGrazieDiga
                                 306 415;...  % FlomartLeGrazie
                                 297 444;...  % Chienti1  
                                 288 472;...  % FlomartPrimaFiastra
                                 284 482;...  % FlomartDopoFiastra
                                 288 515;...  % Chienti2
                                 285 522;...  % FlomartTorrione
                                 285 532;...  % FlomartMonteCosaro
                                 277 565];    % FoceChienti   
        % define the names of each selected section:
        nomi_sezioni{1}='FlomartPolverina_Chienti';
        nomi_sezioni{2}='Diga_Borgiano';
        nomi_sezioni{3}='Diga_LeGrazie';
        nomi_sezioni{4}='FlomartLeGrazie_Chienti';
        nomi_sezioni{5}='Chienti1_Chienti';
        nomi_sezioni{6}='FlomartPrimaFiastra_Chienti';
        nomi_sezioni{7}='FlomartDopoFiastra_Chienti';
        nomi_sezioni{8}='Chienti2_Chienti';
        nomi_sezioni{9}='FlomartTorrione_Chienti';
        nomi_sezioni{10}='FlomartMonteCosaro_Chienti';
        nomi_sezioni{11}='FoceChienti_Chienti';
        
        % define name of catchment and section:
        bacino_sezione{1} = 'Chienti_FlomartPolverina'; 
        bacino_sezione{2} = 'Chienti_Diga_Borgiano'; 
        bacino_sezione{3} = 'Chienti_Diga_LeGrazie'; 
        bacino_sezione{4} = 'Chienti_FlomartLeGrazie';
        bacino_sezione{5} = 'Chienti_Chienti1'; 
        bacino_sezione{6} = 'Chienti_FlomartPrimaFiastra';
        bacino_sezione{7} = 'Chienti_FlomartDopoFiastra';
        bacino_sezione{8} = 'Chienti_Chienti2'; 
        bacino_sezione{9} = 'Chienti_FlomartTorrione';
        bacino_sezione{10} = 'Chienti_FlomartMonteCosaro';
        bacino_sezione{11} = 'Chienti_FoceChienti';
        
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
if strcmp(domain_name,'Foglia') | strcmp(domain_name,'Chienti') 
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




%% ! please, verify the geo reference system code in the inputs!  EPSG_domain  
% salvo geotiff e .mat aree competenza (controllare):
save([path_preparation_data,'/AreeCompetenza_', domain_name,'.mat'], 'AreeCompetenza')
geotiffwrite([path_preparation_data,'/prova_area', domain_name,'.tif'], double(AreeCompetenza),R, 'CoordRefSysCode',['EPSG:',EPSG_domain])


    
% Save obtained information into mat file (this file is used as static input by Flomart application):
save([sPathOutputInfoMat, '/info_',domain_name,'.mat'], ...
    'LonLL','LonUR','LatLL','LatUR',...
    'AreeCompetenza','mappa_aree_allargata','mappa_aree','Lat_dominio_UTM','Lon_dominio_UTM','a1dQindex','nomi_sezioni_sort', ...
     'indici_sort','bacino_sezione_sort', 'EPSG_domain', 'drainage_area_section_km2');



 
% % if an error occurs while saving mat file (mainly because of size of variabile AreeCompetenza then save into h5df format ('-v7.3'):
% save([sPathOutputInfoMat, '/info_',domain_name,'.mat'], ...
%     'LonLL','LonUR','LatLL','LatUR',...
%     'AreeCompetenza','mappa_aree_allargata','mappa_aree','Lat_dominio_UTM','Lon_dominio_UTM','a1dQindex','nomi_sezioni_sort', ...
%      'indici_sort','bacino_sezione_sort', 'EPSG_domain', 'drainage_area_section_km2', '-v7.3');




%% create and save flood map for tr=0 (with zeros):
HazardMap0 = zeros(size(A));
geotiffwrite([sPathHazardData, '/', domain_name, '/', name_hazardmaps, sprintf('%03.0f',0), type_hazardmaps], double(HazardMap0),R, 'CoordRefSysCode',['EPSG:',EPSG_domain])



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 EXAMPLE OF OPERATIONAL FLOOD MAP GENERATION:            % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% insert values of peak discharge for each selected section:
discharge_values_sections =[800, 2000, 3050]; %Foglia
discharge_values_sections =[300, 200, 300, 400, 500, 1000]; %Entella
discharge_values_sections =[200, 300, 400, 500]; %Scrivia
discharge_values_sections =[200, 300, 400, 500, 1000]; % Chienti


examp_map = example_generation_flood_map(domain_name, sPathHazardData, name_hazardmaps, path_preparation_data, ...
                                         discharge_values_sections, quante_sez, Qindice_dominio, ...
                                         sezioni_indici_relativi_corr2, nomi_sezioni, ...
                                         AreeCompetenza, R, TR_max, TR_min, EPSG_domain);



                                     
                                     

                                     
                                     

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIND RELATION Q vs. T:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sFileName_QT2 = [path_preparation_data, '/Mappe_Regionalizzazione_Q_v2/Mappa_quantili_T2.asc'];
sFileName_QT5 = [path_preparation_data, '/Mappe_Regionalizzazione_Q_v2/Mappa_quantili_T5.asc'];
sFileName_QT10 = [path_preparation_data, '/Mappe_Regionalizzazione_Q_v2/Mappa_quantili_T10.asc'];
sFileName_QT20 = [path_preparation_data, '/Mappe_Regionalizzazione_Q_v2/Mappa_quantili_T20.asc'];
sFileName_QT50 = [path_preparation_data, '/Mappe_Regionalizzazione_Q_v2/Mappa_quantili_T50.asc'];
sFileName_QT100 = [path_preparation_data, '/Mappe_Regionalizzazione_Q_v2/Mappa_quantili_T100.asc'];
sFileName_QT150 = [path_preparation_data, '/Mappe_Regionalizzazione_Q_v2/Mappa_quantili_T150.asc'];
sFileName_QT200 = [path_preparation_data, '/Mappe_Regionalizzazione_Q_v2/Mappa_quantili_T200.asc'];
sFileName_QT500 = [path_preparation_data, '/Mappe_Regionalizzazione_Q_v2/Mappa_quantili_T500.asc'];
% Q(T=2):
[a2dMap_QT2, a2dCoord_QT2] = arcgridread(sFileName_QT2);
a2iQT2 = a2dMap_QT2;
% Q(T=5):
[a2dMap_QT5, a2dCoord_QT5] = arcgridread(sFileName_QT5);
a2iQT5 = a2dMap_QT5;
% Q(T=10):
[a2dMap_QT10, a2dCoord_QT10] = arcgridread(sFileName_QT10);
a2iQT10 = a2dMap_QT10;
% Q(T=20):
[a2dMap_QT20, a2dCoord_QT20] = arcgridread(sFileName_QT20);
a2iQT20 = a2dMap_QT20;
% Q(T=50):
[a2dMap_QT50, a2dCoord_QT50] = arcgridread(sFileName_QT50);
a2iQT50 = a2dMap_QT50;
% Q(T=100):
[a2dMap_QT100, a2dCoord_QT100] = arcgridread(sFileName_QT100);
a2iQT100 = a2dMap_QT100;
% Q(T=150):
[a2dMap_QT150, a2dCoord_QT150] = arcgridread(sFileName_QT150);
a2iQT150 = a2dMap_QT150;
% Q(T=200):
[a2dMap_QT200, a2dCoord_QT200] = arcgridread(sFileName_QT200);
a2iQT200 = a2dMap_QT200;
% Q(T=500):
[a2dMap_QT500, a2dCoord_QT500] = arcgridread(sFileName_QT500);
a2iQT500 = a2dMap_QT500;
    
%
for i=1:size(sezioni_indici_relativi,1)
        a1dQT2(i)   = a2iQT2(sezioni_indici_relativi(i,1),sezioni_indici_relativi(i,2));
        a1dQT5(i)   = a2iQT5(sezioni_indici_relativi(i,1),sezioni_indici_relativi(i,2));
        a1dQT10(i)  = a2iQT10(sezioni_indici_relativi(i,1),sezioni_indici_relativi(i,2));
        a1dQT20(i)  = a2iQT20(sezioni_indici_relativi(i,1),sezioni_indici_relativi(i,2));
        a1dQT50(i)  = a2iQT50(sezioni_indici_relativi(i,1),sezioni_indici_relativi(i,2));
        a1dQT100(i) = a2iQT100(sezioni_indici_relativi(i,1),sezioni_indici_relativi(i,2));
        a1dQT150(i) = a2iQT150(sezioni_indici_relativi(i,1),sezioni_indici_relativi(i,2));
        a1dQT200(i) = a2iQT200(sezioni_indici_relativi(i,1),sezioni_indici_relativi(i,2));
        a1dQT500(i) = a2iQT500(sezioni_indici_relativi(i,1),sezioni_indici_relativi(i,2));
end

figure
RT= [0 2 5 10 20 50 100 150 200 500];
for i=length(L):-1:1
        a1dQT2_sort(i) = a1dQT2(indici_sort(i));               
        a1dQT5_sort(i) = a1dQT5(indici_sort(i));  
        a1dQT10_sort(i) = a1dQT10(indici_sort(i));    
        a1dQT20_sort(i) = a1dQT20(indici_sort(i));         
        a1dQT50_sort(i) = a1dQT50(indici_sort(i));        
        a1dQT100_sort(i) = a1dQT100(indici_sort(i));              
        a1dQT150_sort(i) = a1dQT150(indici_sort(i));              
        a1dQT200_sort(i) = a1dQT200(indici_sort(i));          
        a1dQT500_sort(i) = a1dQT500(indici_sort(i));             
        
        QQ = [0 a1dQT2_sort(i) a1dQT5_sort(i) a1dQT10_sort(i) a1dQT20_sort(i) ... 
              a1dQT50_sort(i) a1dQT100_sort(i) a1dQT150_sort(i) a1dQT200_sort(i) a1dQT500_sort(i)];
        plot(RT, QQ, 'o-', 'markersize',6, 'MarkerFaceColor', 'k');
        
        x = log(RT);
        y = QQ;
        hold on;
        a=[];
        for k =1:length(x)
            a=[a ; x(k) 1];
        end
        c = a\y';
        hold on;
        x =min(RT):0.01:max(RT);
        y=c(1)*log(x)+c(2);
        %plot(x,y, 'k');

        txt2 = nomi_sezioni_sort{i};
        txt2 = strrep(txt2,'_','\_');
        if c(2) < 0
            txt = [txt2, '   \rightarrow     y = ', num2str(c(1)), ' log(x) - ' num2str(abs(c(2)))];
        else
            txt = [txt2, '   \rightarrow     y = ', num2str(c(1)), ' log(x) + ' num2str(c(2))];
        end    

        %text(200, c(1)*log(200) + c(2) + 50, txt);
        display(txt);
        display(QQ);
        %legend({txt}, 'Location','southeast');

        
        Qlimit_sort(i) = c(1)*log(1) + c(2);   % lower limit of the relation
        hold on;
end
grid off;
title('FOGLIA')
%title('CHIENTI')
xlabel('Return Time T') 
ylabel('Discharge Q_T [m^3/s]')                                     
     


x = exp((y - c(2))/c(1));
