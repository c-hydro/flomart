function x = example_generation_flood_map(domain_name, sPathHazardData, name_hazardmaps, path_mappe_output, discharge_values_sections, quante_sez, ...
                                          Qindice_dominio, sezioni_indici_relativi_corr2, nomi_sezioni, ...
                                          AreeCompetenza, R, TR_max, TR_min, EPSG_domain)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcolo scenari:
for uuu=1:quante_sez
    TT(uuu)=round(exp(((Qindice_dominio(sezioni_indici_relativi_corr2(uuu,1), sezioni_indici_relativi_corr2(uuu,2)).*0.5239) + ...
                  discharge_values_sections(uuu))./(Qindice_dominio(sezioni_indici_relativi_corr2(uuu,1),sezioni_indici_relativi_corr2(uuu,2)).*1.0433)));
              
    disp(['sezione "',char(nomi_sezioni(uuu)),'": Portata=',num2str(discharge_values_sections(uuu)),' --> scenario=',num2str(TT(uuu))])
end

%% GENERAZIONE MAPPE: crezione degli scenari di unione
% dimensioni della griglia idraulica
[nrows,ncols]=size(AreeCompetenza);   
% cell array con gli indici matrice delle aree di competenza
for i=1:quante_sez
    indici_aree{i}=find(AreeCompetenza==i);                             
end

% costruzione scenario:
mappa_flood=uint16(zeros(nrows,ncols));
TT=max(1,min(TR_max,TT));
T_unici=unique(TT);
for i=1:length(T_unici)
    T=T_unici(i);
    if T<TR_min
        continue
    end
    indici_T=find(TT==T);

    % .mat
    name_file_hazard_read = [sPathHazardData, '/', domain_name,'/', name_hazardmaps, sprintf('%03.0f',T),'.mat'];
    load(name_file_hazard_read, 'mappa_h');
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

x = mappa_flood;
% scrittura geotiff finale
geotiffwrite([path_mappe_output,'/scenario_operativo-2-',domain_name,'.tif'],mappa_flood, R, 'CoordRefSysCode',['EPSG:',EPSG_domain])
disp('DONE! SCENARIO GENERATED!!!')
