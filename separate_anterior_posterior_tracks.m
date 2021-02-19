function [dataTracks_rotated,dataTracks_rotated_afterFurrow,dataTracks_rotated_beforeFurrow] = separate_anterior_posterior_tracks...
    (dataTracks_rotated,dataTracks_rotated_afterFurrow,dataTracks_rotated_beforeFurrow,regionXlength,regionXlimit,furrow_position,...
    reference,polarity_metaphase_50,global_analysis,vector_datum)

% to separate tracks that are in the anterior balstomere or in the( furrow_position.image_start_detection - param.delta_furrowDetection_anaphaseOnset *param.sp6 ) > 100
% posterior one

% parameter reference used to choose with temporal windows are studied;
% according to that
% dataTracks_rotated_afterFurrow,dataTracks_rotated_beforeFurrow, input
% args could be metaphase, anaphase, late, prometaphase

global pathMainDirectory
global param
global general_param

if nargin < 8
    polarity_metaphase_50 = 0;
end

% if global analysis of residency time analysis, no saving of the different
% datatrcks for timing periods
if nargin < 9
    global_analysis = 0;
end
if nargin < 10
    vector_datum = [];
end

numTracks_anterior = 0;
numTracks_posterior = 0;
numTracks_outside = 0;


if reference == 1 % reference is furrow onset
    reference_time = furrow_position.image_start_detection;
    reference_time_bis = param.sp3+1;
elseif reference == 2 % reference is m2a
    reference_time = furrow_position.image_start_detection - param.delta_furrowDetection_anaphaseOnset *param.sp6 ;
    dataTracks_rotated_anaphase = dataTracks_rotated_afterFurrow;
    dataTracks_rotated_prometaphase = dataTracks_rotated_beforeFurrow;
    if furrow_position.image_start_detection - param.delta_furrowDetection_anaphaseOnset *param.sp6 +  param.delta_late_anaphase *param.sp6 < param.sp3
        reference_time_bis = furrow_position.image_start_detection - param.delta_furrowDetection_anaphaseOnset *param.sp6 +  ...
            param.delta_late_anaphase *param.sp6 ;
    else
        reference_time_bis = param.sp3+1;
    end
elseif reference == 3 % 2 min after m2a
    reference_time = furrow_position.image_start_detection - param.delta_furrowDetection_anaphaseOnset *param.sp6 + param.delta_late_anaphase *param.sp6;
    dataTracks_rotated_late = dataTracks_rotated_afterFurrow;
    dataTracks_rotated_anaphase = dataTracks_rotated_beforeFurrow;
    reference_time_bis = param.sp3+1;
elseif reference == 4 % 2 min before m2a
    reference_time = furrow_position.image_start_detection - param.delta_furrowDetection_anaphaseOnset *param.sp6 - param.delta_early_metaphase *param.sp6;
    dataTracks_rotated_metaphase = dataTracks_rotated_afterFurrow;
    dataTracks_rotated_prophase = dataTracks_rotated_beforeFurrow;
    reference_time_bis = furrow_position.image_start_detection - param.delta_furrowDetection_anaphaseOnset *param.sp6 ;
end

if polarity_metaphase_50 == 1 % to set the boundary between anteror and posterior to 50% before the furrow ingression
    furrow_position.xCoordinate.timeDependence(1:furrow_position.image_start_detection-1) = ...
        furrow_position.xCoordinate.mean * 50/ furrow_position.percent_length.mean;
end


%%

% full recording
% disp('ID:')
% disp(param.stem_name);
% sp2_bkp = param.sp2;
for i = 1 : dataTracks_rotated.entireEmbryo.numTracks
    if dataTracks_rotated.entireEmbryo.indexXStart(1,i) - param.sp2 <=0
        param.sp2 = 0;
    end
    %     if ( sp2_bkp < 26 && sp2_bkp > 1 ) && ...
    %             ( vector_datum(1,1) <= 2018 || ...
    %             ( vector_datum(1,1) == 2019 && ( vector_datum(1,2) < 5 || (vector_datum(1,2) == 6 && vector_datum(1,3) < 24 ) ) ) )
    %         param.sp2 = sp2_bkp + 26;
    %     end
    %     disp(param.stem_name);
    
    if dataTracks_rotated.entireEmbryo.lengthTracks(i) >= general_param.cortex_analysis.minLength * param.sp6
        
        if ( dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i) < ...
                furrow_position.xCoordinate.timeDependence(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) ) ...
                && ( regionXlimit.entireEmbryo.nbR10(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,1) <= ...
                dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i)  )
            numTracks_anterior = numTracks_anterior +1;
            lengthTracks_anterior(numTracks_anterior) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
            indexXEnd_anterior(numTracks_anterior) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
            indexXStart_anterior(numTracks_anterior) = dataTracks_rotated.entireEmbryo.indexXStart(i);
            tracksX_anterior(:,numTracks_anterior) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
            tracksY_anterior(:,numTracks_anterior) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
        elseif ( dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i) < ...
                regionXlimit.entireEmbryo.nbR10(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,1) + ...
                10*regionXlength.entireEmbryo.nbR10(1,dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) ) ...
                && ( furrow_position.xCoordinate.timeDependence(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) ...
                <= dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i) )
            numTracks_posterior = numTracks_posterior +1;
            lengthTracks_posterior(numTracks_posterior) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
            indexXEnd_posterior(numTracks_posterior) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
            indexXStart_posterior(numTracks_posterior) = dataTracks_rotated.entireEmbryo.indexXStart(i);
            tracksX_posterior(:,numTracks_posterior) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
            tracksY_posterior(:,numTracks_posterior) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
        else
            numTracks_outside = numTracks_outside +1;
            
        end
        
    end
end


%save data from tracking in structure file
if numTracks_anterior > 0
    dataTracks_rotated.anterior.tracksX = tracksX_anterior;
    dataTracks_rotated.anterior.tracksY = tracksY_anterior;
    dataTracks_rotated.anterior.lengthTracks = lengthTracks_anterior;
    dataTracks_rotated.anterior.indexXStart = indexXStart_anterior;
    dataTracks_rotated.anterior.indexXEnd = indexXEnd_anterior;
    dataTracks_rotated.anterior.numTracks = numTracks_anterior;
else
    dataTracks_rotated.anterior.numTracks = 0;
end

if numTracks_posterior > 0
    dataTracks_rotated.posterior.tracksX = tracksX_posterior;
    dataTracks_rotated.posterior.tracksY = tracksY_posterior;
    dataTracks_rotated.posterior.lengthTracks = lengthTracks_posterior;
    dataTracks_rotated.posterior.indexXStart = indexXStart_posterior;
    dataTracks_rotated.posterior.indexXEnd = indexXEnd_posterior;
    dataTracks_rotated.posterior.numTracks = numTracks_posterior;
    dataTracks_rotated.numTracks_outside = numTracks_outside;
else
    dataTracks_rotated.posterior.numTracks = 0;
end

if global_analysis == 0
    name = strcat('dataTracks_rotated-', short_name, '.mat');
    save(fullfile(pathMainDirectory,name), '-struct', 'dataTracks_rotated','-v7.3');
    %save(fullfile(pathMainDirectory,name), '-struct', 'dataTracks_rotated');
end

clear tracksX_anterior
clear tracksY_anterior
clear lengthTracks_anterior
clear indexXStart_anterior
clear indexXEnd_anterior

clear tracksX_posterior
clear tracksY_posterior
clear lengthTracks_posterior
clear indexXStart_posterior
clear indexXEnd_posterior


%%

if reference ~= 0
    
    % after apparition furrow at cortex or metaphase/anaphase
    
    numTracks_anterior = 0;
    numTracks_posterior = 0;
    numTracks_outside = 0;
    numTracks_anterior2 = 0;
    numTracks_posterior2 = 0;
    numTracks_outside2 = 0;
    if reference == 1
        dataTracks_rotated_beforeFurrow.anterior.numTracks = NaN;
        dataTracks_rotated_beforeFurrow.posterior.numTracks = NaN;
        dataTracks_rotated_afterFurrow.anterior.numTracks = NaN;
        dataTracks_rotated_afterFurrow.posterior.numTracks = NaN;
    elseif reference ==2
        dataTracks_rotated_prometaphase.anterior.numTracks = NaN;
        dataTracks_rotated_prometaphase.posterior.numTracks = NaN;
        dataTracks_rotated_anaphase.anterior.numTracks = NaN;
        dataTracks_rotated_anaphase.posterior.numTracks = NaN;
    elseif reference ==3
        dataTracks_rotated_anaphase.anterior.numTracks = NaN;
        dataTracks_rotated_anaphase.posterior.numTracks = NaN;
        dataTracks_rotated_late.anterior.numTracks = NaN;
        dataTracks_rotated_late.posterior.numTracks = NaN;
    elseif reference ==4
        dataTracks_rotated_metaphase.anterior.numTracks = NaN;
        dataTracks_rotated_metaphase.posterior.numTracks = NaN;
        dataTracks_rotated_prophase.anterior.numTracks = NaN;
        dataTracks_rotated_prophase.posterior.numTracks = NaN;        
    end
    
    if ~isnan(furrow_position.image_start_detection)
        
        for i = 1 : dataTracks_rotated.entireEmbryo.numTracks
            
            if dataTracks_rotated.entireEmbryo.lengthTracks(i) >= general_param.cortex_analysis.minLength * param.sp6
                
                if dataTracks_rotated.entireEmbryo.indexXStart(1,i) > reference_time && ... % above refernce time being m2a or furron ingression or start late-anaphase
                        dataTracks_rotated.entireEmbryo.indexXStart(1,i) < reference_time_bis  % below last image or late anaphase or m2a
                    if ( dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i) < ...
                            furrow_position.xCoordinate.timeDependence(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) ) ...
                            && ( regionXlimit.entireEmbryo.nbR10(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,1) <= ...
                            dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i)  )
                        numTracks_anterior = numTracks_anterior +1;
                        lengthTracks_anterior(numTracks_anterior) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
                        indexXEnd_anterior(numTracks_anterior) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
                        indexXStart_anterior(numTracks_anterior) = dataTracks_rotated.entireEmbryo.indexXStart(i);
                        tracksX_anterior(:,numTracks_anterior) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
                        tracksY_anterior(:,numTracks_anterior) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
                    elseif ( dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i) < ...
                            regionXlimit.entireEmbryo.nbR10(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,1) + ...
                            10*regionXlength.entireEmbryo.nbR10(1,dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) ) ...
                            && ( furrow_position.xCoordinate.timeDependence(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) ...
                            <= dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i) )
                        numTracks_posterior = numTracks_posterior +1;
                        lengthTracks_posterior(numTracks_posterior) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
                        indexXEnd_posterior(numTracks_posterior) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
                        indexXStart_posterior(numTracks_posterior) = dataTracks_rotated.entireEmbryo.indexXStart(i);
                        tracksX_posterior(:,numTracks_posterior) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
                        tracksY_posterior(:,numTracks_posterior) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
                    else
                        numTracks_outside = numTracks_outside +1;
                    end
                elseif dataTracks_rotated.entireEmbryo.indexXStart(1,i) < reference_time
                    if ( dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i) < ...
                            furrow_position.xCoordinate.timeDependence(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) ) ...
                            && ( regionXlimit.entireEmbryo.nbR10(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,1) <= ...
                            dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i)  )
                        numTracks_anterior2 = numTracks_anterior2 +1;
                        lengthTracks_anterior2(numTracks_anterior2) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
                        indexXEnd_anterior2(numTracks_anterior2) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
                        indexXStart_anterior2(numTracks_anterior2) = dataTracks_rotated.entireEmbryo.indexXStart(i);
                        tracksX_anterior2(:,numTracks_anterior2) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
                        tracksY_anterior2(:,numTracks_anterior2) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
                    elseif ( dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i) < ...
                            regionXlimit.entireEmbryo.nbR10(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,1) + ...
                            10*regionXlength.entireEmbryo.nbR10(1,dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) ) ...
                            && ( furrow_position.xCoordinate.timeDependence(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) ...
                            <= dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i) )
                        numTracks_posterior2 = numTracks_posterior2 +1;
                        lengthTracks_posterior2(numTracks_posterior2) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
                        indexXEnd_posterior2(numTracks_posterior2) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
                        indexXStart_posterior2(numTracks_posterior2) = dataTracks_rotated.entireEmbryo.indexXStart(i);
                        tracksX_posterior2(:,numTracks_posterior2) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
                        tracksY_posterior2(:,numTracks_posterior2) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
                    else
                        numTracks_outside2 = numTracks_outside2 +1;
                    end
                end
                
            end
        end
        
        if reference == 1
            
            %save data from tracking in structure file
            dataTracks_rotated_afterFurrow.anterior.tracksX = tracksX_anterior;
            dataTracks_rotated_afterFurrow.anterior.tracksY = tracksY_anterior;
            dataTracks_rotated_afterFurrow.anterior.lengthTracks = lengthTracks_anterior;
            dataTracks_rotated_afterFurrow.anterior.indexXStart = indexXStart_anterior;
            dataTracks_rotated_afterFurrow.anterior.indexXEnd = indexXEnd_anterior;
            dataTracks_rotated_afterFurrow.anterior.numTracks = numTracks_anterior;
            
            dataTracks_rotated_afterFurrow.posterior.tracksX = tracksX_posterior;
            dataTracks_rotated_afterFurrow.posterior.tracksY = tracksY_posterior;
            dataTracks_rotated_afterFurrow.posterior.lengthTracks = lengthTracks_posterior;
            dataTracks_rotated_afterFurrow.posterior.indexXStart = indexXStart_posterior;
            dataTracks_rotated_afterFurrow.posterior.indexXEnd = indexXEnd_posterior;
            dataTracks_rotated_afterFurrow.posterior.numTracks = numTracks_posterior;
            dataTracks_rotated_afterFurrow.numTracks_outside = numTracks_outside;
            
            if global_analysis == 0
                name2 = strcat('dataTracks_rotated_afterFurrowDetection-', short_name, '.mat');
                save(fullfile(pathMainDirectory,name2), '-struct', 'dataTracks_rotated_afterFurrow');
            end
            
            dataTracks_rotated_beforeFurrow.anterior.tracksX = tracksX_anterior2;
            dataTracks_rotated_beforeFurrow.anterior.tracksY = tracksY_anterior2;
            dataTracks_rotated_beforeFurrow.anterior.lengthTracks = lengthTracks_anterior2;
            dataTracks_rotated_beforeFurrow.anterior.indexXStart = indexXStart_anterior2;
            dataTracks_rotated_beforeFurrow.anterior.indexXEnd = indexXEnd_anterior2;
            dataTracks_rotated_beforeFurrow.anterior.numTracks = numTracks_anterior2;
            
            dataTracks_rotated_beforeFurrow.posterior.tracksX = tracksX_posterior2;
            dataTracks_rotated_beforeFurrow.posterior.tracksY = tracksY_posterior2;
            dataTracks_rotated_beforeFurrow.posterior.lengthTracks = lengthTracks_posterior2;
            dataTracks_rotated_beforeFurrow.posterior.indexXStart = indexXStart_posterior2;
            dataTracks_rotated_beforeFurrow.posterior.indexXEnd = indexXEnd_posterior2;
            dataTracks_rotated_beforeFurrow.posterior.numTracks = numTracks_posterior2;
            dataTracks_rotated_beforeFurrow.numTracks_outside = numTracks_outside2;
            
            if global_analysis == 0
                name2 = strcat('dataTracks_rotated_beforeFurrowDetection-', short_name, '.mat');
                save(fullfile(pathMainDirectory,name2), '-struct', 'dataTracks_rotated_beforeFurrow');
            end
            
        elseif reference == 2
            
            %save data from tracking in structure file
            if numTracks_anterior > 0
                dataTracks_rotated_anaphase.anterior.tracksX = tracksX_anterior;
                dataTracks_rotated_anaphase.anterior.tracksY = tracksY_anterior;
                dataTracks_rotated_anaphase.anterior.lengthTracks = lengthTracks_anterior;
                dataTracks_rotated_anaphase.anterior.indexXStart = indexXStart_anterior;
                dataTracks_rotated_anaphase.anterior.indexXEnd = indexXEnd_anterior;
                dataTracks_rotated_anaphase.anterior.numTracks = numTracks_anterior;
            end
            
            if numTracks_posterior > 0
                dataTracks_rotated_anaphase.posterior.tracksX = tracksX_posterior;
                dataTracks_rotated_anaphase.posterior.tracksY = tracksY_posterior;
                dataTracks_rotated_anaphase.posterior.lengthTracks = lengthTracks_posterior;
                dataTracks_rotated_anaphase.posterior.indexXStart = indexXStart_posterior;
                dataTracks_rotated_anaphase.posterior.indexXEnd = indexXEnd_posterior;
                dataTracks_rotated_anaphase.posterior.numTracks = numTracks_posterior;
            end
            dataTracks_rotated_anaphase.numTracks_outside = numTracks_outside;
            
            if global_analysis == 0
                name2 = strcat('dataTracks_rotated_anaphase-', short_name, '.mat');
                save(fullfile(pathMainDirectory,name2), '-struct', 'dataTracks_rotated_anaphase');
            end
            
            if numTracks_anterior2 > 0
                dataTracks_rotated_prometaphase.anterior.tracksX = tracksX_anterior2;
                dataTracks_rotated_prometaphase.anterior.tracksY = tracksY_anterior2;
                dataTracks_rotated_prometaphase.anterior.lengthTracks = lengthTracks_anterior2;
                dataTracks_rotated_prometaphase.anterior.indexXStart = indexXStart_anterior2;
                dataTracks_rotated_prometaphase.anterior.indexXEnd = indexXEnd_anterior2;
                dataTracks_rotated_prometaphase.anterior.numTracks = numTracks_anterior2;
            end
            
            if numTracks_posterior2 > 0
                dataTracks_rotated_prometaphase.posterior.tracksX = tracksX_posterior2;
                dataTracks_rotated_prometaphase.posterior.tracksY = tracksY_posterior2;
                dataTracks_rotated_prometaphase.posterior.lengthTracks = lengthTracks_posterior2;
                dataTracks_rotated_prometaphase.posterior.indexXStart = indexXStart_posterior2;
                dataTracks_rotated_prometaphase.posterior.indexXEnd = indexXEnd_posterior2;
                dataTracks_rotated_prometaphase.posterior.numTracks = numTracks_posterior2;
            end
            dataTracks_rotated_prometaphase.numTracks_outside = numTracks_outside2;
            
            if global_analysis == 0
                name2 = strcat('dataTracks_rotated_prometaphase-', short_name, '.mat');
                save(fullfile(pathMainDirectory,name2), '-struct', 'dataTracks_rotated_prometaphase');
            end
            
            dataTracks_rotated_afterFurrow = dataTracks_rotated_anaphase;
            dataTracks_rotated_beforeFurrow = dataTracks_rotated_prometaphase;
            
        elseif reference == 3
            
            %save data from tracking in structure file
            if numTracks_anterior
                dataTracks_rotated_late.anterior.tracksX = tracksX_anterior;
                dataTracks_rotated_late.anterior.tracksY = tracksY_anterior;
                dataTracks_rotated_late.anterior.lengthTracks = lengthTracks_anterior;
                dataTracks_rotated_late.anterior.indexXStart = indexXStart_anterior;
                dataTracks_rotated_late.anterior.indexXEnd = indexXEnd_anterior;
                dataTracks_rotated_late.anterior.numTracks = numTracks_anterior;
            end
            
            if numTracks_posterior
                dataTracks_rotated_late.posterior.tracksX = tracksX_posterior;
                dataTracks_rotated_late.posterior.tracksY = tracksY_posterior;
                dataTracks_rotated_late.posterior.lengthTracks = lengthTracks_posterior;
                dataTracks_rotated_late.posterior.indexXStart = indexXStart_posterior;
                dataTracks_rotated_late.posterior.indexXEnd = indexXEnd_posterior;
                dataTracks_rotated_late.posterior.numTracks = numTracks_posterior;
            end
            dataTracks_rotated_late.numTracks_outside = numTracks_outside;
            
            if global_analysis == 0
                name2 = strcat('dataTracks_rotated_late-', short_name, '.mat');
                save(fullfile(pathMainDirectory,name2), '-struct', 'dataTracks_rotated_late');
            end
            
            %save data from tracking in structure file
            if numTracks_anterior2 > 0
                dataTracks_rotated_anaphase.anterior.tracksX = tracksX_anterior2;
                dataTracks_rotated_anaphase.anterior.tracksY = tracksY_anterior2;
                dataTracks_rotated_anaphase.anterior.lengthTracks = lengthTracks_anterior2;
                dataTracks_rotated_anaphase.anterior.indexXStart = indexXStart_anterior2;
                dataTracks_rotated_anaphase.anterior.indexXEnd = indexXEnd_anterior2;
                dataTracks_rotated_anaphase.anterior.numTracks = numTracks_anterior2;
            end
            
            if numTracks_posterior2 > 0
                dataTracks_rotated_anaphase.posterior.tracksX = tracksX_posterior2;
                dataTracks_rotated_anaphase.posterior.tracksY = tracksY_posterior2;
                dataTracks_rotated_anaphase.posterior.lengthTracks = lengthTracks_posterior2;
                dataTracks_rotated_anaphase.posterior.indexXStart = indexXStart_posterior2;
                dataTracks_rotated_anaphase.posterior.indexXEnd = indexXEnd_posterior2;
                dataTracks_rotated_anaphase.posterior.numTracks = numTracks_posterior2;
            end
            dataTracks_rotated_anaphase.numTracks_outside = numTracks_outside;
            
            if global_analysis == 0
                name2 = strcat('dataTracks_rotated_anaphase-', short_name, '.mat');
                save(fullfile(pathMainDirectory,name2), '-struct', 'dataTracks_rotated_anaphase');
            end
            
            dataTracks_rotated_afterFurrow = dataTracks_rotated_late;
            dataTracks_rotated_beforeFurrow = dataTracks_rotated_anaphase;
            
        elseif reference == 4
            
            %save data from tracking in structure file
            if numTracks_anterior > 0
                dataTracks_rotated_metaphase.anterior.tracksX = tracksX_anterior;
                dataTracks_rotated_metaphase.anterior.tracksY = tracksY_anterior;
                dataTracks_rotated_metaphase.anterior.lengthTracks = lengthTracks_anterior;
                dataTracks_rotated_metaphase.anterior.indexXStart = indexXStart_anterior;
                dataTracks_rotated_metaphase.anterior.indexXEnd = indexXEnd_anterior;
                dataTracks_rotated_metaphase.anterior.numTracks = numTracks_anterior;
            end
            
            if numTracks_posterior > 0
                dataTracks_rotated_metaphase.posterior.tracksX = tracksX_posterior;
                dataTracks_rotated_metaphase.posterior.tracksY = tracksY_posterior;
                dataTracks_rotated_metaphase.posterior.lengthTracks = lengthTracks_posterior;
                dataTracks_rotated_metaphase.posterior.indexXStart = indexXStart_posterior;
                dataTracks_rotated_metaphase.posterior.indexXEnd = indexXEnd_posterior;
                dataTracks_rotated_metaphase.posterior.numTracks = numTracks_posterior;
            end
            dataTracks_rotated_metaphase.numTracks_outside = numTracks_outside;
            
            if global_analysis == 0
                name2 = strcat('dataTracks_rotated_metaphase-', short_name, '.mat');
                save(fullfile(pathMainDirectory,name2), '-struct', 'dataTracks_rotated_metaphase');
            end
            
            if numTracks_anterior2 > 0
                dataTracks_rotated_prophase.anterior.tracksX = tracksX_anterior2;
                dataTracks_rotated_prophase.anterior.tracksY = tracksY_anterior2;
                dataTracks_rotated_prophase.anterior.lengthTracks = lengthTracks_anterior2;
                dataTracks_rotated_prophase.anterior.indexXStart = indexXStart_anterior2;
                dataTracks_rotated_prophase.anterior.indexXEnd = indexXEnd_anterior2;
                dataTracks_rotated_prophase.anterior.numTracks = numTracks_anterior2;
            end
            
            if numTracks_posterior2 > 0
                dataTracks_rotated_prophase.posterior.tracksX = tracksX_posterior2;
                dataTracks_rotated_prophase.posterior.tracksY = tracksY_posterior2;
                dataTracks_rotated_prophase.posterior.lengthTracks = lengthTracks_posterior2;
                dataTracks_rotated_prophase.posterior.indexXStart = indexXStart_posterior2;
                dataTracks_rotated_prophase.posterior.indexXEnd = indexXEnd_posterior2;
                dataTracks_rotated_prophase.posterior.numTracks = numTracks_posterior2;
            end
            dataTracks_rotated_prophase.numTracks_outside = numTracks_outside2;
            
            if global_analysis == 0
                name2 = strcat('dataTracks_rotated_prophase-', short_name, '.mat');
                save(fullfile(pathMainDirectory,name2), '-struct', 'dataTracks_rotated_prophase');
            end
            
            dataTracks_rotated_afterFurrow = dataTracks_rotated_metaphase;
            dataTracks_rotated_beforeFurrow = dataTracks_rotated_prophase;           
            
        end
        
        clear tracksX_anterior
        clear tracksY_anterior
        clear lengthTracks_anterior
        clear indexXStart_anterior
        clear indexXEnd_anterior
        
        clear tracksX_posterior
        clear tracksY_posterior
        clear lengthTracks_posterior
        clear indexXStart_posterior
        clear indexXEnd_posterior
        
        
        clear tracksX_anterior2
        clear tracksY_anterior2
        clear lengthTracks_anterior2
        clear indexXStart_anterior2
        clear indexXEnd_anterior2
        
        clear tracksX_posterior2
        clear tracksY_posterior2
        clear lengthTracks_posterior2
        clear indexXStart_posterior2
        clear indexXEnd_posterior2
        
    end
    
end


end

