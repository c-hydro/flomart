
%%%%%%%%%%
%% INPUTs
%%%%%%%%%%
path_code = '/home/matteo/Documents/CIMA_projects/RT_FloodMapping/flomart_testing/code/flomart_preprocessing';
sPathLandData = '/home/matteo/Documents/CIMA_projects/RT_FloodMapping/flomart_testing/data/Marche/data_static/domain_data';
nome_dominio = 'Foglia';  % define the name of the domain.
sPathGeoData = '/home/matteo/Documents/CIMA_projects/RT_FloodMapping/flomart_testing/data/Marche/data_static/geo_data';
% percorso dove si trovano le simulazioni:
sPathHazardData='/home/matteo/Documents/CIMA_projects/RT_FloodMapping/flomart_testing/data/Marche/data_static/hazard_data';
sPathTelemacData='/home/matteo/Documents/CIMA_projects/RT_FloodMapping/flomart_testing/data/Marche/data_static/telemac_data';
sFileNameOutput = [sPathGeoData, '/Data_MarcheDomain.mat'];
% percorso dove si trovano gli scenari:
path_mappe_input=[pwd];
path_mappe_output=[pwd];
% composizione nome file scenari:
nome_hazmaps=[nome_dominio,'_WD_max_Q'];% nome delle hazard maps di partenza senza le tre cifre finali (es:  "Foglia_WD_max_Q001.tif" )
nome_hazmaps1='.tif';  % format of the hazard map files
% nome della mappa risultato
nome_floodmap_out='FloodMap';
% tempo di ritorno al di sotto (strettamente) del quale la mappa di hazard si annulla
Tsoglia=1;
% tempo di ritorno massimo (eventuali valori superiori vengono saturati a questo valore)
TR_max=500;
% coordinate system EPSG: 
EPSG_domain = '32633';









%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Creazione file Data_MarcheDomain.mat:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%
% If already available import the Data_*_*Domain.mat file:
load(sFileNameOutput);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Otherwhise create a new Data_*_*Domain.mat file:

%% names of files with gridded info (choice, dem, ...):
%MARCHE:
sFileName_choice = '/fp_marche/MarcheDomain.choice.txt';
sFileName_area = '/fp_marche/MarcheDomain.area.txt';
sFileName_cell = '/fp_marche/MarcheDomain.areacell.txt';
sFileName_dem = '/fp_marche/MarcheDomain.dem.txt';
sFileName_lon = '/fp_marche/MarcheDomain.lon.txt';
sFileName_lat = '/fp_marche/MarcheDomain.lat.txt';
sFileName_pnt = '/fp_marche/MarcheDomain.pnt.txt';
sFileName_lon = '/gridded_marche/marche.lon.txt';
sFileName_lat = '/gridded_marche/marche.lat.txt';

% file .mat containing variable 'a2dQindex':
namfmat2load_bis = [sPathLandData, '/fp_marche/MarcheDomainStatistica_FPI_classes.mat'];







%% LOADING
%load(namfmat2load)


%% DEFINING VARIABLES OF FILE .MAT:
% Choice:
[a2dMap_choice, a2dCoord_choice] = arcgridread([sPathLandData,sFileName_choice]);
a2iChoice = a2dMap_choice;
a2iChoice(isnan(a2iChoice)) = -1;  %replace all NaN with -1

% Area:
[a2dMap_area, a2dCoord_area] = arcgridread([sPathLandData,sFileName_area]);
a2dArea = a2dMap_area;

% Cell:
[a2dMap_cell, a2dCoord_cell] = arcgridread([sPathLandData,sFileName_cell]);
a2dCelle = a2dMap_cell;
a2dCelle(isnan(a2dCelle)) = 1; %replace all NaN with 1

% Dem:
[a2dMap_dem, a2dCoord_dem] = arcgridread([sPathLandData,sFileName_dem]);
a2dDem = a2dMap_dem;
a2dDem(isnan(a2dDem)) = -99;  %replace all NaN with -99

% Pointers:
[a2dMap_pnt, a2dCoord_pnt] = arcgridread([sPathLandData,sFileName_pnt]);
a2iPunt = a2dMap_pnt;
a2iPunt(isnan(a2iPunt)) = 0;   %replace all NaN with 0

% Lon:
% [a2dMap_lon, a2dCoord_lon] = arcgridread([sPathLandData,sFileName_lon]);
% Londem = a2dMap_lon;
% Londem(isnan(Londem)) = -99;

% Lat:
% [a2dMap_lat, a2dCoord_lat] = arcgridread([sPathLandData,sFileName_lat]);
% Latdem = a2dMap_lat;
% Latdem(isnan(Latdem)) = -99;

% important: if the imported "lat" and "lon" have -9999 values you need to 
% modify them in order to avoid any -9999. Create the grid manually:
[lon, lat] = meshgrid(12.0717-0.005006:0.005006:14.0591 + 0.005006*2,  ...
                      42.5264:0.005006:44.0733+0.005006);
lat = flip(lat,1);
%lat = lat.';
Latdem = lat;
Londem = lon;





%% Qindex:
% [a2dMapQindice, a2dCoord_Qindice] = arcgridread([sPathLandData,sFileName_Qindice]);
% a2dQindice = a2dMapQindice;
% a2dQindice(isnan(a2dQindice)) = 0;

% or if you have it in a .mat file:
load(namfmat2load_bis)
a2dQindice = a2dQindex;


% if not known:
%%%  create Qindex layer from historic hydrographs (> 10 years).
%%%  ...
%%%  ... TO DO


%%
% save final .mat file:
save(sFileNameOutput,'a2dDem','a2dCelle','a2dArea','a2iChoice', 'Latdem', 'Londem', 'a2iPunt', 'a2dQindice')















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Creazione mat file "Aree_finali_Foglia.mat":
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% define the domain grid: 
if strcmp(nome_dominio,'Foglia')
        % define the corners:
        Lat_min=43.776;
        Lat_max=43.926;
        Lon_min=12.485;
        Lon_max=12.929;
        
        temp=find(((Londem(1,:)-Lon_min)<0)&(Londem(1,:)-Lon_min)>-100);
        indice_y_min=temp(end);clear temp
        [x,y] =find(((Londem(:,:)-Lon_max)<0)&(Londem(:,:)-Lon_max)>-100);
        indice_y_max=y(end);clear temp
       
        temp =find(((Latdem(:,1)-Lat_min)>0)&(Latdem(:,1)-Lat_min)>-100);
        indice_x_max=temp(end);clear temp
        temp=find(((Latdem(:,1)-Lat_max)>0)&(Latdem(:,1)-Lat_max)>-100);
        indice_x_min=temp(end);clear temp

        save([sPathGeoData, '/Aree_competenza_',nome_dominio,'.mat'],'Lat_min','Lat_max',...
            'Lon_min','Lon_max', 'indice_y_min','indice_y_max','indice_x_max','indice_x_min')  
end



        


        



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


% % number of sections (USER ADDS BY PROMPT):
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% quante_sez=input('Quante sezioni?')
% figure 
% imagesc(Area_dominio)
% caxis([2 20])
% % scelgo dove sono le sezioni
% [pointclick_x,pointclick_y] = ginput(quante_sez);
% close all
% % costruisco file delle sezioni
% sezioni_indici_relativi=[round(pointclick_y),round(pointclick_x)];
% % nomi sezioni
% for ind=1:quante_sez
%         nomi_sezioni{ind}='NoName';
% end



% Instead, here chosen sections are hard-coded:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(nome_dominio,'Foglia')
        % define the number of control sections:
        quante_sez=3;
        % define indexes of the selected sections [Y X]:
        sezioni_indici_relativi=[25 6;...  % Bronzo   
                                 17 57;... % Montecchio
                                 5 84];    % OUT    
        % define the names of each selected section:                     
        nomi_sezioni{1}='Bronzo';
        nomi_sezioni{2}='Montecchio';
        nomi_sezioni{3}='OUT';
        
        % plot to verify:
        figure 
        imagesc(Area_dominio)
        caxis([2 50])
        hold on
        for indicew=1:length(sezioni_indici_relativi)
            plot(sezioni_indici_relativi(indicew,2), sezioni_indici_relativi(indicew,1),'o','markersize',12)
        end
end

     
     

%% AREE DRENATE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the drainage area for each selected section:
cd(path_code)
aree = aree_drenate(Punt_dominio, sezioni_indici_relativi);

% Delete superposed areas:
mappa_aree=zeros(size(Punt_dominio));
L=cellfun(@length,aree);
[ordinati,indici_sort]=sort(L);
for i=length(L):-1:1
        valori=unique(mappa_aree(indici_sort(i)));
        mappa_aree(aree{indici_sort(i)})=i;
end
    
% Verify areas with a PLOT:
figure
imagesc(mappa_aree)
hold on
[canalix,canaliy]=find(Choice_dominio==1);
for indicew=1:length(canalix)
    plot(canaliy(indicew),canalix(indicew),'.k','markersize',6)
end
for indicew=1:length(sezioni_indici_relativi)
    plot(sezioni_indici_relativi(indicew,2),sezioni_indici_relativi(indicew,1),'or','markersize',6)
end
    




%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% modification in order to increase the dominium (borders)
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
for indicew=1:length(sezioni_indici_relativi)
        plot(sezioni_indici_relativi(indicew,2),sezioni_indici_relativi(indicew,1),'or','markersize',6)
end



    
% IF NECESSARY:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% manually modify "mappa_aree_allargata" in order to include the river 
% delta or other zones known to be flooded that are not included in the domain:
if nome_dominio=='Foglia'
        mappa_aree_allargata(2:6,86:89) = quante_sez;
        mappa_aree_allargata(2:4,84:85) = quante_sez;
        mappa_aree_allargata(2:4,81:83) = quante_sez; 
          
elseif nome_dominio=='Graveglia'
        disp('No allargamento !');
end
    
% verify PLOT:
figure
imagesc(mappa_aree_allargata)
hold on
[canalix,canaliy]=find(Choice_dominio==1);
for indicew=1:length(canalix)
        plot(canaliy(indicew),canalix(indicew),'.k','markersize',6)
end
for indicew=1:length(sezioni_indici_relativi)
        plot(sezioni_indici_relativi(indicew,2),sezioni_indici_relativi(indicew,1),'or','markersize',6)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    
    
%%
% recover info of metric grid UTM33N from the hazard maps:
name_file_read=[sPathTelemacData, '/', nome_dominio,'/',nome_hazmaps,num2str(TR_max),nome_hazmaps1];
[A, R]=geotiffread(name_file_read);


    


% Rigriglio su griglia metrica UTM33N:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function 'meshgrid(x,y)' transforms the domain specified by vectors x and y into arrays X and Y, 
% which can be used to evaluate functions of two variables and three-dimensional mesh/surface plots. 
% The rows of the output array X are copies of the vector x;
% The columns of the output array Y are copies of the vector y.
[new_x,new_y]= meshgrid(R.XWorldLimits(1):R.CellExtentInWorldX:R.XWorldLimits(2) - R.CellExtentInWorldX, ...
                        R.YWorldLimits(1):R.CellExtentInWorldY:R.YWorldLimits(2) - R.CellExtentInWorldY);

% convert coord lat/lon to coordinates utm32: 
[coord_left, coord_bottom] = latlon2utm(Lat_min, Lon_min);
[coord_right, coord_top]   = latlon2utm(Lat_max, Lon_max);


% get the resolution of the domain.mat in degrees and meters:
res_Lon = ceil((coord_right - coord_left)/size(mappa_aree_allargata,2));
res_Lat = ceil((coord_top - coord_bottom)/size(mappa_aree_allargata,1));
[Lon_dominio_UTM, Lat_dominio_UTM] = meshgrid(coord_left:res_Lon:coord_right, ...
                                              coord_bottom:res_Lat:coord_top);
        
% Interpolate the sample grid data using the griddata() and the 'nearest' method:
% griddata(x, y, v, xq, yq) fits a surface of the form v = f(x,y) to the scattered data in the vectors (x,y,v). 
% The griddata function interpolates the surface at the query points specified by (xq,yq) 
% and returns the interpolated values. The surface always passes through the data points defined by x and y.
mappa_aree_new = griddata(Lat_dominio_UTM, Lon_dominio_UTM, ...
                          mappa_aree_allargata, ...
                          new_y, new_x, ...
                          'nearest');

% verify PLOT:
figure
imagesc(mappa_aree_new)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
 



%% remember to modify the geo reference system code!  CoordRefSysCode  
%%
% salvo geotiff e .mat aree competenza (controllare):
geotiffwrite([sPathGeoData,'/prova_area', nome_dominio,'.tif'], double(mappa_aree_new),R, 'CoordRefSysCode',['EPSG:',EPSG_domain])
% salvo il file 'Aree_finali_',nome_dominio,'.mat'
save([sPathGeoData, '/Aree_finali_',nome_dominio,'.mat'], 'new_x','new_y',...
    'mappa_aree_new','R','mappa_aree_allargata','mappa_aree','Lat_dominio_UTM','Lon_dominio_UTM')
      


















%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    ESEMPIO GENERAZIONE MAPPE:                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% insert values of peak discharge for each selected section:
portate=[500, 800, 850];

% calcolo scenari:
for uuu=1:quante_sez
    TT(uuu)=round(exp(((Qindice_dominio(sezioni_indici_relativi(uuu,1),sezioni_indici_relativi(uuu,2)).*0.5239) + ...
                  portate(uuu))./(Qindice_dominio(sezioni_indici_relativi(uuu,1),sezioni_indici_relativi(uuu,2)).*1.0433)));
              
    disp(['sezione "',char(nomi_sezioni(uuu)),'": Portata=',num2str(portate(uuu)),' --> scenario=',num2str(TT(uuu))])
end

%% GENERAZIONE MAPPE: crezione degli scenari di unione
% dimensioni della griglia idraulica
[nrows,ncols]=size(mappa_aree_new);   
% cell array con gli indici matrice delle aree di competenza
for i=1:quante_sez
    indici_aree{i}=find(mappa_aree_new==i);                             
end



% costruzione scenario:
mappa_flood=uint16(zeros(nrows,ncols));
TT=max(1,min(TR_max,TT));
T_unici=unique(TT);
for i=1:length(T_unici)
    T=T_unici(i);
    if T<Tsoglia
        continue
    end
    indici_T=find(TT==T);
    load([sPathHazardData, '/', nome_dominio,'/',nome_dominio,'_WD_max_T',sprintf('%03.0f',T),'.mat'],'mappa_h');
    mappa_flood(indici_aree{indici_T})=mappa_h((indici_aree{indici_T}));
    clear mappa_h
end



% definizione mappa finale
mappa_flood=double(mappa_flood)./1000;




% % verifica visiva (sconsigliata per domini grandi)
figure
imagesc(mappa_flood)



% pcolor(new_x,new_y,flipud(mappa_flood))
% shading flat
% caxis([0 7])
% hold on
% colormap([1 1 1;...
%           .8 1 1;...
%           .6 .94 1;...
%           .04 .96 .96;...
%           .04 .67 .96;...
%            0 .45 .74;...
%            0 0 1])
% for i=1:quante_sez
%     plot(sezioni{i,5},sezioni{i,4},'.k','markersize',10)
% end


% scrittura geotiff finale
geotiffwrite([path_mappe_output,'/scenario_operativo-2-',nome_dominio,'.tif'],mappa_flood, R, 'CoordRefSysCode',['EPSG:',EPSG_domain])
disp('SCENARIO GENERATO!!!')
