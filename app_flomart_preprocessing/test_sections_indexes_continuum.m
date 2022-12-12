%%
clc; clear;


%%
%%%%%%%%%%%%%%
%   INPUTS   %
%%%%%%%%%%%%%%

% percorso dove si trova questo script matlab (e i suoi moduli) per il preprocessing in Flomart:
path_code = '/home/matteo/Documents/CIMA_projects/RT_FloodMapping/flomart-2.0.0_test/app_flomart_preprocessing';
% percorso dove salvare i risultati:
path_preparation_data ="/home/matteo/Documents/CIMA_projects/Marche_Qgis/shapefile_marche_REPAIR";

% percorso dove cercare dati statici geografici con griglia: 
%sPathGridGeoData = '/home/matteo/Documents/CIMA_projects/RT_FloodMapping/data/data_Marche_Foglia/data_static/PREPARATION/fp_marche';
sPathGridGeoData = '/home/matteo/Documents/CIMA_projects/Marche_Qgis/StaticiContinuum_Marche';
sPathShapefile = '/home/matteo/Documents/CIMA_projects/Marche_Qgis/20221019_shapefile_marche';




%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read geographical gridnformation:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% names of files with gridded info (choice, dem, ...):
% MARCHE:
% sFileName_choice = [sPathGridGeoData, '/MarcheDomain.choice.txt'];      from FP ITALIA
% sFileName_area = [sPathGridGeoData,'/MarcheDomain.area.txt'];
sFileName_choice = [sPathGridGeoData, '/marche.choice.txt'];             % from REGIONE MARCHE
sFileName_area = [sPathGridGeoData, '/marche.area.txt'];
sFileName_shapefile =  [sPathShapefile, '/fp_sections_marche.shp']; 


%% DEFINING VARIABLES OF FILE .MAT:
% Choice:
[a2dMap_choice, a2dCoord_choice] = arcgridread(sFileName_choice);
a2iChoice = a2dMap_choice;
a2iChoice(isnan(a2iChoice)) = -1;  %replace all NaN with -1

% Area:
[a2dMap_area, a2dCoord_area] = arcgridread(sFileName_area);
a2dArea = a2dMap_area;

% % Plot to verify:
% imagesc(Choice_dominio);
% caxis([2 30]);

SHAPE = shaperead(sFileName_shapefile);
indexes_X = [SHAPE.HMC_X];
indexes_Y = [SHAPE.HMC_Y];
indexes_names = {SHAPE.SEC_NAME};
indexes_basins = {SHAPE.BASIN};
indexes_domains = {SHAPE.DOMAIN};
indexes_area = [SHAPE.AREA];
indexes_Qthr1 = [SHAPE.Q_THR1];
indexes_Qthr2 = [SHAPE.Q_THR2];
indexes_Qthr3 =[SHAPE.Q_THR3];
indexes_baseflow = [SHAPE.BASEFLOW];
indexes_secrs = [SHAPE.SEC_RS];


sect(:, 1) = indexes_X;
sect(:, 2) = indexes_Y;
sezioni_indici_relativi = sect;
          


% plot to verify sections on starting grid:
figure 
imagesc(a2dArea)
caxis([2 50])
hold on
for indicew=1:length(sezioni_indici_relativi)
    plot(sezioni_indici_relativi(indicew,2), sezioni_indici_relativi(indicew,1),'or','markersize',13, 'LineWidth',5)
    text(sezioni_indici_relativi(indicew,2), sezioni_indici_relativi(indicew,1), indexes_names{indicew}, 'FontSize', 11, 'Color', 'y')
end






%% 
%%%%%%%%%%%%%%%%%%%% SAVE
%Quindi, salvare in file section .txt.
fid = fopen([sPathShapefile, '/marche.info_section_NEW.txt'], 'wt');
for i = 1:length(indexes_X)
    if strcmp(indexes_domains{i},'Nera')
        display('Section of Nera basin, removed!');
    else
        fprintf(fid,'%.0f', indexes_X(i));
        fprintf(fid,'%s', " ");
        fprintf(fid,'%.0f', indexes_Y(i));
        fprintf(fid,'%s', " ");
        fprintf(fid,indexes_basins{i});
        fprintf(fid,'%s', " ");
        fprintf(fid,indexes_names{i});
        fprintf(fid,'%s', " ");
        fprintf(fid,num2str(indexes_secrs(i)));
        fprintf(fid,'%s', " ");
        fprintf(fid,'%.1f',indexes_area(i));
        fprintf(fid,'%s', " ");
        fprintf(fid,'%.1f',indexes_Qthr1(i));
        fprintf(fid,'%s', " ");
        fprintf(fid,'%.1f', indexes_Qthr2(i));
        fprintf(fid,'%s', " ");
        fprintf(fid,indexes_domains{i});
        fprintf(fid,'%s', " ");
        fprintf(fid,'%.1f',(indexes_baseflow(i)));
        fprintf(fid,'\n');
    end
end
fclose(fid);






%%
% indice_extra = 81;
% plot(sezioni_indici_relativi(indice_extra,2), sezioni_indici_relativi(indice_extra,1),'ow','markersize',15, 'LineWidth',5)
% indice_extra2 = 85;
% plot(sezioni_indici_relativi(indice_extra2,2), sezioni_indici_relativi(indice_extra2,1),'ow','markersize',15, 'LineWidth',5)
% 
% indice_extra3 = 83;
% plot(sezioni_indici_relativi(indice_extra3,2), sezioni_indici_relativi(indice_extra3,1),'ow','markersize',15, 'LineWidth',5)
% 
% indice_extra4 = 45;
% plot(sezioni_indici_relativi(indice_extra4,2), sezioni_indici_relativi(indice_extra4,1),'ow','markersize',15, 'LineWidth',5)
% 

% % verifica aree
% for i=1:size(sezioni_indici_relativi,1)
%     AreaBas(i)=a2dArea(sezioni_indici_relativi(i,1),sezioni_indici_relativi(i,2));
% end










%% extra tool (Create sections indexes) from lat lon
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % transform lat/lon coordinates into matricial continuum indixes:
% tmp = a2dMap_choice;
% masknan = nan(size(tmp));
% masknan(tmp==1) = 333;
% c = readtable([sPathShapefile, '/sections_chienti_higher30cells_selected2.csv']);
% x_points = c.lon;
% y_points = c.lat;
% namfascii = sFileName_choice;
% [rr, cc] = Get_MatrixCoordinates_FromPoints_NoNaN(x_points, y_points, namfascii, masknan);
% indexes_names = [c.info];
% 
% 
% % plot to verify sections on starting grid:
% figure 
% imagesc(a2dArea)
% caxis([2 50])
% hold on
% for indicew=1:length(rr)
%     plot(cc(indicew), rr(indicew),'or','markersize',13, 'LineWidth',5)
%     text(cc(indicew), rr(indicew), indexes_names{indicew}, 'FontSize', 7, 'Color', 'y')
% end



%%
% %Quindi, salvare le coordinate matrice [rr, cc] nel file “rigacolonna.txt”.
% fid = fopen([path_preparation_data, '/rigacolonna.txt'], 'wt');
% for i = 1:length(rr)
%     fprintf(fid,'%.0f', rr(i));
%     fprintf(fid,'%s', " ");
%     fprintf(fid,'%.0f', cc(i));
%     fprintf(fid,'\n');
% end
% fclose(fid);

