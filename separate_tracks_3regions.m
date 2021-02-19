function [dataTracks_rotated,dataTracks_rotated_afterFurrow,dataTracks_rotated_beforeFurrow,dataTracks_rotated_anaphase,...
    dataTracks_rotated_metaphase,dataTracks_rotated_late,dataTracks_rotated_prometaphase, dataTracks_rotated_prophase] ...
    = separate_tracks_3regions(dataTracks_rotated,xStart_3regions,xEnd_3regions,furrow_position,late_separation,early_separation,vector_datum,time_reference_choice)

if nargin < 5
    late_separation = 0;
end
if nargin < 6
    early_separation = 0;
end
if nargin < 7
    vector_datum = [];
end
if nargin < 8
    time_reference_choice = 0;
end

% separate tracks into three regions in first part
% then if furrow detected, do same before and after furrow detection
% image

global pathMainDirectory
global param

numTracks_region1 = 0;
numTracks_region2 = 0;
numTracks_region3 = 0;
numTracks_outside = 0;


%%

% full recording
disp(param.stem_name)
sp2_bkp = param.sp2;
for i = 1 : dataTracks_rotated.entireEmbryo.numTracks
    if dataTracks_rotated.entireEmbryo.indexXStart(1,i) - param.sp2 <=0
        param.sp2 = 0;
    end
%     if ( sp2_bkp < 26 && sp2_bkp > 1 ) && ...
%             ( vector_datum(1,1) <= 2018 || ...
%             ( vector_datum(1,1) == 2019 && ( vector_datum(1,2) < 5 || (vector_datum(1,2) == 6 && vector_datum(1,3) < 24 ) ) ) ) 
%        param.sp2 = sp2_bkp + 26; 
%     end
%     if sp2_bkp > 25
%        param.sp2 = sp2_bkp - 26; 
%     end

    if ( dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i) < ...
            xEnd_3regions.region1(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) ) ...
            && ( xStart_3regions.region1(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) <= ...
            dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i)  )
        numTracks_region1 = numTracks_region1 +1;
        lengthTracks_region1(numTracks_region1) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
        indexXEnd_region1(numTracks_region1) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
        indexXStart_region1(numTracks_region1) = dataTracks_rotated.entireEmbryo.indexXStart(i);
        tracksX_region1(:,numTracks_region1) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
        tracksY_region1(:,numTracks_region1) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
    elseif ( dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i) < ...
            xEnd_3regions.region2(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) ) ...
            && ( xStart_3regions.region2(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) <= ...
            dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i)  )
        numTracks_region2 = numTracks_region2 +1;
        lengthTracks_region2(numTracks_region2) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
        indexXEnd_region2(numTracks_region2) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
        indexXStart_region2(numTracks_region2) = dataTracks_rotated.entireEmbryo.indexXStart(i);
        tracksX_region2(:,numTracks_region2) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
        tracksY_region2(:,numTracks_region2) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
    elseif ( dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i) < ...
            xEnd_3regions.region3(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) ) ...
            && ( xStart_3regions.region3(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) <= ...
            dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i)  )
        numTracks_region3 = numTracks_region3 +1;
        lengthTracks_region3(numTracks_region3) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
        indexXEnd_region3(numTracks_region3) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
        indexXStart_region3(numTracks_region3) = dataTracks_rotated.entireEmbryo.indexXStart(i);
        tracksX_region3(:,numTracks_region3) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
        tracksY_region3(:,numTracks_region3) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
    else
        numTracks_outside = numTracks_outside +1;
        
    end
end


%save data from tracking in structure file
dataTracks_rotated.region1.tracksX = tracksX_region1;
dataTracks_rotated.region1.tracksY = tracksY_region1;
dataTracks_rotated.region1.lengthTracks = lengthTracks_region1;
dataTracks_rotated.region1.indexXStart = indexXStart_region1;
dataTracks_rotated.region1.indexXEnd = indexXEnd_region1;
dataTracks_rotated.region1.numTracks = numTracks_region1;

dataTracks_rotated.region2.tracksX = tracksX_region2;
dataTracks_rotated.region2.tracksY = tracksY_region2;
dataTracks_rotated.region2.lengthTracks = lengthTracks_region2;
dataTracks_rotated.region2.indexXStart = indexXStart_region2;
dataTracks_rotated.region2.indexXEnd = indexXEnd_region2;
dataTracks_rotated.region2.numTracks = numTracks_region2;

dataTracks_rotated.region3.tracksX = tracksX_region3;
dataTracks_rotated.region3.tracksY = tracksY_region3;
dataTracks_rotated.region3.lengthTracks = lengthTracks_region3;
dataTracks_rotated.region3.indexXStart = indexXStart_region3;
dataTracks_rotated.region3.indexXEnd = indexXEnd_region3;
dataTracks_rotated.region3.numTracks = numTracks_region3;

dataTracks_rotated.numTracks_outside = numTracks_outside;


name = strcat('dataTracks_rotated_3regions-', short_name, '.mat');
save(fullfile(pathMainDirectory,name), '-struct', 'dataTracks_rotated');

clear tracksX_region1
clear tracksY_region1
clear lengthTracks_region1
clear indexXStart_region1
clear indexXEnd_region1

clear tracksX_region2
clear tracksY_region2
clear lengthTracks_region2
clear indexXStart_region2
clear indexXEnd_region2

clear tracksX_region3
clear tracksY_region3
clear lengthTracks_region3
clear indexXStart_region3
clear indexXEnd_region3


%% before/after ingression furrow onset

if time_reference_choice == 0
    
    % after apparition furrow at cortex
    
    numTracks_region1a = 0;
    numTracks_region2a= 0;
    numTracks_region3a = 0;
    numTracks_outsideA = 0;
    numTracks_region1b = 0;
    numTracks_region2b= 0;
    numTracks_region3b = 0;
    numTracks_outsideB = 0;
    
    
    dataTracks_rotated_beforeFurrow.region1.numTracks = NaN;
    dataTracks_rotated_beforeFurrow.region2.numTracks = NaN;
    dataTracks_rotated_beforeFurrow.region3.numTracks = NaN;
    dataTracks_rotated_afterFurrow.region1.numTracks = NaN;
    dataTracks_rotated_afterFurrow.region2.numTracks = NaN;
    dataTracks_rotated_afterFurrow.region3.numTracks = NaN;
    
    if ~isnan(furrow_position.image_start_detection) && furrow_position.image_start_detection < param.sp3
        
        for i = 1 : dataTracks_rotated.entireEmbryo.numTracks
            if dataTracks_rotated.entireEmbryo.indexXStart(1,i) > furrow_position.image_start_detection
                
                if ( dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i) < ...
                        xEnd_3regions.region1(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) ) ...
                        && ( xStart_3regions.region1(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) <= ...
                        dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i)  )
                    numTracks_region1a = numTracks_region1a +1;
                    lengthTracks_region1a(numTracks_region1a) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
                    indexXEnd_region1a(numTracks_region1a) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
                    indexXStart_region1a(numTracks_region1a) = dataTracks_rotated.entireEmbryo.indexXStart(i);
                    tracksX_region1a(:,numTracks_region1a) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
                    tracksY_region1a(:,numTracks_region1a) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
                elseif ( dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i) < ...
                        xEnd_3regions.region2(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) ) ...
                        && ( xStart_3regions.region2(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) <= ...
                        dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i)  )
                    numTracks_region2a = numTracks_region2a +1;
                    lengthTracks_region2a(numTracks_region2a) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
                    indexXEnd_region2a(numTracks_region2a) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
                    indexXStart_region2a(numTracks_region2a) = dataTracks_rotated.entireEmbryo.indexXStart(i);
                    tracksX_region2a(:,numTracks_region2a) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
                    tracksY_region2a(:,numTracks_region2a) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
                elseif ( dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i) < ...
                        xEnd_3regions.region3(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) ) ...
                        && ( xStart_3regions.region3(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) <= ...
                        dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i)  )
                    numTracks_region3a = numTracks_region3a +1;
                    lengthTracks_region3a(numTracks_region3a) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
                    indexXEnd_region3a(numTracks_region3a) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
                    indexXStart_region3a(numTracks_region3a) = dataTracks_rotated.entireEmbryo.indexXStart(i);
                    tracksX_region3a(:,numTracks_region3a) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
                    tracksY_region3a(:,numTracks_region3a) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
                else
                    numTracks_outsideA = numTracks_outsideA +1;
                    
                end
                
            else
                
                if ( dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i) < ...
                        xEnd_3regions.region1(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) ) ...
                        && ( xStart_3regions.region1(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) <= ...
                        dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i)  )
                    numTracks_region1b = numTracks_region1b +1;
                    lengthTracks_region1b(numTracks_region1b) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
                    indexXEnd_region1b(numTracks_region1b) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
                    indexXStart_region1b(numTracks_region1b) = dataTracks_rotated.entireEmbryo.indexXStart(i);
                    tracksX_region1b(:,numTracks_region1b) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
                    tracksY_region1b(:,numTracks_region1b) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
                elseif ( dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i) < ...
                        xEnd_3regions.region2(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) ) ...
                        && ( xStart_3regions.region2(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) <= ...
                        dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i)  )
                    numTracks_region2b = numTracks_region2b +1;
                    lengthTracks_region2b(numTracks_region2b) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
                    indexXEnd_region2b(numTracks_region2b) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
                    indexXStart_region2b(numTracks_region2b) = dataTracks_rotated.entireEmbryo.indexXStart(i);
                    tracksX_region2b(:,numTracks_region2b) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
                    tracksY_region2b(:,numTracks_region2b) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
                elseif ( dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i) < ...
                        xEnd_3regions.region3(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) ) ...
                        && ( xStart_3regions.region3(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) <= ...
                        dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i)  )
                    numTracks_region3b = numTracks_region3b +1;
                    lengthTracks_region3b(numTracks_region3b) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
                    indexXEnd_region3b(numTracks_region3b) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
                    indexXStart_region3b(numTracks_region3b) = dataTracks_rotated.entireEmbryo.indexXStart(i);
                    tracksX_region3b(:,numTracks_region3b) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
                    tracksY_region3b(:,numTracks_region3b) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
                else
                    numTracks_outsideB = numTracks_outsideB +1;
                    
                end
                
                
            end
        end
        
        %save data from tracking in structure file
        if numTracks_region1a ~= 0
            dataTracks_rotated_afterFurrow.region1.tracksX = tracksX_region1a;
            dataTracks_rotated_afterFurrow.region1.tracksY = tracksY_region1a;
            dataTracks_rotated_afterFurrow.region1.lengthTracks = lengthTracks_region1a;
            dataTracks_rotated_afterFurrow.region1.indexXStart = indexXStart_region1a;
            dataTracks_rotated_afterFurrow.region1.indexXEnd = indexXEnd_region1a;
            dataTracks_rotated_afterFurrow.region1.numTracks = numTracks_region1a;
            
            clear tracksX_region1a
            clear tracksY_region1a
            clear lengthTracks_region1a
            clear indexXStart_region1a
            clear indexXEnd_region1a
        else
            dataTracks_rotated_afterFurrow.region1.numTracks = 0;
        end
        
        if numTracks_region2a ~= 0
            dataTracks_rotated_afterFurrow.region2.tracksX = tracksX_region2a;
            dataTracks_rotated_afterFurrow.region2.tracksY = tracksY_region2a;
            dataTracks_rotated_afterFurrow.region2.lengthTracks = lengthTracks_region2a;
            dataTracks_rotated_afterFurrow.region2.indexXStart = indexXStart_region2a;
            dataTracks_rotated_afterFurrow.region2.indexXEnd = indexXEnd_region2a;
            dataTracks_rotated_afterFurrow.region2.numTracks = numTracks_region2a;
            
            clear tracksX_region2a
            clear tracksY_region2a
            clear lengthTracks_region2a
            clear indexXStart_region2a
            clear indexXEnd_region2a
        else
            dataTracks_rotated_afterFurrow.region2.numTracks = 0;
        end
        
        if numTracks_region2a ~= 0
            dataTracks_rotated_afterFurrow.region3.tracksX = tracksX_region3a;
            dataTracks_rotated_afterFurrow.region3.tracksY = tracksY_region3a;
            dataTracks_rotated_afterFurrow.region3.lengthTracks = lengthTracks_region3a;
            dataTracks_rotated_afterFurrow.region3.indexXStart = indexXStart_region3a;
            dataTracks_rotated_afterFurrow.region3.indexXEnd = indexXEnd_region3a;
            dataTracks_rotated_afterFurrow.region3.numTracks = numTracks_region3a;
            
            clear tracksX_region3a
            clear tracksY_region3a
            clear lengthTracks_region3a
            clear indexXStart_region3a
            clear indexXEnd_region3a
        else
            dataTracks_rotated_afterFurrow.region3.numTracks = 0;
        end
        
        dataTracks_rotated_afterFurrow.numTracks_outside = numTracks_outsideA;
        
        name2 = strcat('dataTracks_rotated_afterFurrowDetection_3regions-', short_name, '.mat');
        save(fullfile(pathMainDirectory,name2), '-struct', 'dataTracks_rotated_afterFurrow');
        
        if numTracks_region1b ~= 0
            dataTracks_rotated_beforeFurrow.region1.tracksX = tracksX_region1b;
            dataTracks_rotated_beforeFurrow.region1.tracksY = tracksY_region1b;
            dataTracks_rotated_beforeFurrow.region1.lengthTracks = lengthTracks_region1b;
            dataTracks_rotated_beforeFurrow.region1.indexXStart = indexXStart_region1b;
            dataTracks_rotated_beforeFurrow.region1.indexXEnd = indexXEnd_region1b;
            dataTracks_rotated_beforeFurrow.region1.numTracks = numTracks_region1b;
            
            clear tracksX_region1b
            clear tracksY_region1b
            clear lengthTracks_region1b
            clear indexXStart_region1b
            clear indexXEnd_region1b
        else
            dataTracks_rotated_beforeFurrow.region1.numTracks = 0;
        end
        
        if numTracks_region2b ~= 0
            dataTracks_rotated_beforeFurrow.region2.tracksX = tracksX_region2b;
            dataTracks_rotated_beforeFurrow.region2.tracksY = tracksY_region2b;
            dataTracks_rotated_beforeFurrow.region2.lengthTracks = lengthTracks_region2b;
            dataTracks_rotated_beforeFurrow.region2.indexXStart = indexXStart_region2b;
            dataTracks_rotated_beforeFurrow.region2.indexXEnd = indexXEnd_region2b;
            dataTracks_rotated_beforeFurrow.region2.numTracks = numTracks_region2b;
            
            clear tracksX_region2b
            clear tracksY_region2b
            clear lengthTracks_region2b
            clear indexXStart_region2b
            clear indexXEnd_region2b
        else
            dataTracks_rotated_beforeFurrow.region2.numTracks = 0;
        end
        
        if numTracks_region3b ~= 0
            dataTracks_rotated_beforeFurrow.region3.tracksX = tracksX_region3b;
            dataTracks_rotated_beforeFurrow.region3.tracksY = tracksY_region3b;
            dataTracks_rotated_beforeFurrow.region3.lengthTracks = lengthTracks_region3b;
            dataTracks_rotated_beforeFurrow.region3.indexXStart = indexXStart_region3b;
            dataTracks_rotated_beforeFurrow.region3.indexXEnd = indexXEnd_region3b;
            dataTracks_rotated_beforeFurrow.region3.numTracks = numTracks_region3b;
            
            clear tracksX_region3b
            clear tracksY_region3b
            clear lengthTracks_region3b
            clear indexXStart_region3b
            clear indexXEnd_region3b
        else
            dataTracks_rotated_beforeFurrow.region3.numTracks = 0;
        end
        
        dataTracks_rotated_beforeFurrow.numTracks_outside = numTracks_outsideB;
        
        name2 = strcat('dataTracks_rotated_beforeFurrowDetection_3regions-', short_name, '.mat');
        save(fullfile(pathMainDirectory,name2), '-struct', 'dataTracks_rotated_beforeFurrow');
        
    end
    
else
    
    dataTracks_rotated_beforeFurrow.region1.numTracks = NaN;
    dataTracks_rotated_beforeFurrow.region2.numTracks = NaN;
    dataTracks_rotated_beforeFurrow.region3.numTracks = NaN;
    dataTracks_rotated_afterFurrow.region1.numTracks = NaN;
    dataTracks_rotated_afterFurrow.region2.numTracks = NaN;
    dataTracks_rotated_afterFurrow.region3.numTracks = NaN;
    
end


%% prophase / metaphase / anaphase / late

% in metaphase and anaphase and late anaphase

numTracks_region1a = 0; % anaphase region 1
numTracks_region2a= 0; % anaphase region 2
numTracks_region3a = 0; % anaphase & region 3
numTracks_outsideA = 0;
numTracks_region1b = 0; % metaphase
numTracks_region2b= 0;
numTracks_region3b = 0;
numTracks_outsideB = 0;
numTracks_region1c = 0; % late
numTracks_region2c= 0;
numTracks_region3c = 0;
numTracks_outsideC = 0;
numTracks_region1d = 0; % prometaphase
numTracks_region2d= 0;
numTracks_region3d = 0;
numTracks_outsideD = 0;
numTracks_region1e = 0; % prophase
numTracks_region2e= 0;
numTracks_region3e = 0;
numTracks_outsideE = 0;

numTracks_prometaphase = 0;
numTracks_prophase = 0;
numTracks_metaphase = 0;
numTracks_anaphase = 0;
numTracks_late = 0;

dataTracks_rotated_prometaphase.entireEmbryo.numTracks = NaN;
dataTracks_rotated_prometaphase.region1.numTracks = NaN;
dataTracks_rotated_prometaphase.region2.numTracks = NaN;
dataTracks_rotated_prometaphase.region3.numTracks = NaN;
dataTracks_rotated_prophase.entireEmbryo.numTracks = NaN;
dataTracks_rotated_prophase.region1.numTracks = NaN;
dataTracks_rotated_prophase.region2.numTracks = NaN;
dataTracks_rotated_prophase.region3.numTracks = NaN;
dataTracks_rotated_metaphase.entireEmbryo.numTracks = NaN;
dataTracks_rotated_metaphase.region1.numTracks = NaN;
dataTracks_rotated_metaphase.region2.numTracks = NaN;
dataTracks_rotated_metaphase.region3.numTracks = NaN;
dataTracks_rotated_anaphase.entireEmbryo.numTracks = NaN;
dataTracks_rotated_anaphase.region1.numTracks = NaN;
dataTracks_rotated_anaphase.region2.numTracks = NaN;
dataTracks_rotated_anaphase.region3.numTracks = NaN;
dataTracks_rotated_late.entireEmbryo.numTracks = NaN;
dataTracks_rotated_late.region1.numTracks = NaN;
dataTracks_rotated_late.region2.numTracks = NaN;
dataTracks_rotated_late.region3.numTracks = NaN;


%-------------------------------------------------------------------------------
% time reference = furrow onset

if time_reference_choice == 0
    
    if ~isnan(furrow_position.image_start_detection) && ...
            ( furrow_position.image_start_detection - param.delta_furrowDetection_anaphaseOnset *param.sp6 > 0 ) && ...
            ( furrow_position.image_start_detection < param.sp3 )
        
        for i = 1 : dataTracks_rotated.entireEmbryo.numTracks
            if dataTracks_rotated.entireEmbryo.indexXStart(1,i) <= ( furrow_position.image_start_detection - param.delta_furrowDetection_anaphaseOnset *param.sp6 )
                
                if early_separation == 1  % to look at prophase / metaphase
                    
                    if dataTracks_rotated.entireEmbryo.indexXStart(1,i) < ( furrow_position.image_start_detection - param.delta_furrowDetection_anaphaseOnset *param.sp6 ) ...
                            && dataTracks_rotated.entireEmbryo.indexXStart(1,i) >= ( furrow_position.image_start_detection - param.delta_furrowDetection_anaphaseOnset *param.sp6 ...
                            - param.delta_early_metaphase *param.sp6 ) % metaphase
                        
                        numTracks_metaphase = numTracks_metaphase +1;
                        lengthTracks_metaphase(numTracks_metaphase) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
                        indexXEnd_metaphase(numTracks_metaphase) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
                        indexXStart_metaphase(numTracks_metaphase) = dataTracks_rotated.entireEmbryo.indexXStart(i);
                        tracksX_metaphase(:,numTracks_metaphase) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
                        tracksY_metaphase(:,numTracks_metaphase) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
                        
                        if ( dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i) < ...
                                xEnd_3regions.region1(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) ) ...
                                && ( xStart_3regions.region1(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) <= ...
                                dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i)  )
                            numTracks_region1b = numTracks_region1b +1;
                            lengthTracks_region1b(numTracks_region1b) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
                            indexXEnd_region1b(numTracks_region1b) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
                            indexXStart_region1b(numTracks_region1b) = dataTracks_rotated.entireEmbryo.indexXStart(i);
                            tracksX_region1b(:,numTracks_region1b) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
                            tracksY_region1b(:,numTracks_region1b) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
                        elseif ( dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i) < ...
                                xEnd_3regions.region2(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) ) ...
                                && ( xStart_3regions.region2(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) <= ...
                                dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i)  )
                            numTracks_region2b = numTracks_region2b +1;
                            lengthTracks_region2b(numTracks_region2b) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
                            indexXEnd_region2b(numTracks_region2b) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
                            indexXStart_region2b(numTracks_region2b) = dataTracks_rotated.entireEmbryo.indexXStart(i);
                            tracksX_region2b(:,numTracks_region2b) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
                            tracksY_region2b(:,numTracks_region2b) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
                        elseif ( dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i) < ...
                                xEnd_3regions.region3(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) ) ...
                                && ( xStart_3regions.region3(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) <= ...
                                dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i)  )
                            numTracks_region3b = numTracks_region3b +1;
                            lengthTracks_region3b(numTracks_region3b) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
                            indexXEnd_region3b(numTracks_region3b) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
                            indexXStart_region3b(numTracks_region3b) = dataTracks_rotated.entireEmbryo.indexXStart(i);
                            tracksX_region3b(:,numTracks_region3b) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
                            tracksY_region3b(:,numTracks_region3b) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
                        else
                            numTracks_outsideB = numTracks_outsideB +1;
                            
                        end
                        
                        
                    elseif dataTracks_rotated.entireEmbryo.indexXStart(1,i) < ( furrow_position.image_start_detection - param.delta_furrowDetection_anaphaseOnset *param.sp6 ...
                            - param.delta_early_metaphase *param.sp6 ) % prophase
                        
                        numTracks_prophase = numTracks_prophase +1;
                        lengthTracks_prophase(numTracks_prophase) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
                        indexXEnd_prophase(numTracks_prophase) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
                        indexXStart_prophase(numTracks_prophase) = dataTracks_rotated.entireEmbryo.indexXStart(i);
                        tracksX_prophase(:,numTracks_prophase) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
                        tracksY_prophase(:,numTracks_prophase) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
                        
                        if ( dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i) < ...
                                xEnd_3regions.region1(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) ) ...
                                && ( xStart_3regions.region1(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) <= ...
                                dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i)  )
                            numTracks_region1e = numTracks_region1e +1;
                            lengthTracks_region1e(numTracks_region1e) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
                            indexXEnd_region1e(numTracks_region1e) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
                            indexXStart_region1e(numTracks_region1e) = dataTracks_rotated.entireEmbryo.indexXStart(i);
                            tracksX_region1e(:,numTracks_region1e) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
                            tracksY_region1e(:,numTracks_region1e) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
                        elseif ( dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i) < ...
                                xEnd_3regions.region2(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) ) ...
                                && ( xStart_3regions.region2(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) <= ...
                                dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i)  )
                            numTracks_region2e = numTracks_region2e +1;
                            lengthTracks_region2e(numTracks_region2e) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
                            indexXEnd_region2e(numTracks_region2e) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
                            indexXStart_region2e(numTracks_region2e) = dataTracks_rotated.entireEmbryo.indexXStart(i);
                            tracksX_region2e(:,numTracks_region2e) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
                            tracksY_region2e(:,numTracks_region2e) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
                        elseif ( dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i) < ...
                                xEnd_3regions.region3(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) ) ...
                                && ( xStart_3regions.region3(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) <= ...
                                dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i)  )
                            numTracks_region3e = numTracks_region3e +1;
                            lengthTracks_region3e(numTracks_region3e) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
                            indexXEnd_region3e(numTracks_region3e) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
                            indexXStart_region3e(numTracks_region3e) = dataTracks_rotated.entireEmbryo.indexXStart(i);
                            tracksX_region3e(:,numTracks_region3e) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
                            tracksY_region3e(:,numTracks_region3e) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
                        else
                            numTracks_outsideE = numTracks_outsideE +1;
                            
                        end
                        
                    end
                else %prometaphase
                    
                    numTracks_prometaphase = numTracks_prometaphase +1;
                    lengthTracks_prometaphase(numTracks_prometaphase) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
                    indexXEnd_prometaphase(numTracks_prometaphase) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
                    indexXStart_prometaphase(numTracks_prometaphase) = dataTracks_rotated.entireEmbryo.indexXStart(i);
                    tracksX_prometaphase(:,numTracks_prometaphase) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
                    tracksY_prometaphase(:,numTracks_prometaphase) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
                    
                    if ( dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i) < ...
                            xEnd_3regions.region1(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) ) ...
                            && ( xStart_3regions.region1(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) <= ...
                            dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i)  )
                        numTracks_region1d = numTracks_region1d +1;
                        lengthTracks_region1d(numTracks_region1d) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
                        indexXEnd_region1d(numTracks_region1d) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
                        indexXStart_region1d(numTracks_region1d) = dataTracks_rotated.entireEmbryo.indexXStart(i);
                        tracksX_region1d(:,numTracks_region1d) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
                        tracksY_region1d(:,numTracks_region1d) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
                    elseif ( dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i) < ...
                            xEnd_3regions.region2(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) ) ...
                            && ( xStart_3regions.region2(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) <= ...
                            dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i)  )
                        numTracks_region2d = numTracks_region2d +1;
                        lengthTracks_region2d(numTracks_region2d) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
                        indexXEnd_region2d(numTracks_region2d) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
                        indexXStart_region2d(numTracks_region2d) = dataTracks_rotated.entireEmbryo.indexXStart(i);
                        tracksX_region2d(:,numTracks_region2d) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
                        tracksY_region2d(:,numTracks_region2d) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
                    elseif ( dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i) < ...
                            xEnd_3regions.region3(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) ) ...
                            && ( xStart_3regions.region3(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) <= ...
                            dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i)  )
                        numTracks_region3d = numTracks_region3d +1;
                        lengthTracks_region3d(numTracks_region3d) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
                        indexXEnd_region3d(numTracks_region3d) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
                        indexXStart_region3d(numTracks_region3d) = dataTracks_rotated.entireEmbryo.indexXStart(i);
                        tracksX_region3d(:,numTracks_region3d) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
                        tracksY_region3d(:,numTracks_region3d) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
                    else
                        numTracks_outsideD = numTracks_outsideD +1;
                        
                    end
                    
                end
            else
                
                if late_separation == 1 % to look at anaphase / late mitosis
                    
                    if dataTracks_rotated.entireEmbryo.indexXStart(1,i) > ( furrow_position.image_start_detection - param.delta_furrowDetection_anaphaseOnset *param.sp6 ) ...
                            && dataTracks_rotated.entireEmbryo.indexXStart(1,i) <= ( furrow_position.image_start_detection - param.delta_furrowDetection_anaphaseOnset *param.sp6 ...
                            + param.delta_late_anaphase *param.sp6 )
                        
                        numTracks_anaphase = numTracks_anaphase +1;
                        lengthTracks_anaphase(numTracks_anaphase) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
                        indexXEnd_anaphase(numTracks_anaphase) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
                        indexXStart_anaphase(numTracks_anaphase) = dataTracks_rotated.entireEmbryo.indexXStart(i);
                        tracksX_anaphase(:,numTracks_anaphase) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
                        tracksY_anaphase(:,numTracks_anaphase) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
                        
                        if ( dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i) < ...
                                xEnd_3regions.region1(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) ) ...
                                && ( xStart_3regions.region1(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) <= ...
                                dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i)  )
                            numTracks_region1a = numTracks_region1a +1;
                            lengthTracks_region1a(numTracks_region1a) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
                            indexXEnd_region1a(numTracks_region1a) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
                            indexXStart_region1a(numTracks_region1a) = dataTracks_rotated.entireEmbryo.indexXStart(i);
                            tracksX_region1a(:,numTracks_region1a) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
                            tracksY_region1a(:,numTracks_region1a) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
                        elseif ( dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i) < ...
                                xEnd_3regions.region2(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) ) ...
                                && ( xStart_3regions.region2(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) <= ...
                                dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i)  )
                            numTracks_region2a = numTracks_region2a +1;
                            lengthTracks_region2a(numTracks_region2a) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
                            indexXEnd_region2a(numTracks_region2a) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
                            indexXStart_region2a(numTracks_region2a) = dataTracks_rotated.entireEmbryo.indexXStart(i);
                            tracksX_region2a(:,numTracks_region2a) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
                            tracksY_region2a(:,numTracks_region2a) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
                        elseif ( dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i) < ...
                                xEnd_3regions.region3(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) ) ...
                                && ( xStart_3regions.region3(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) <= ...
                                dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i)  )
                            numTracks_region3a = numTracks_region3a +1;
                            lengthTracks_region3a(numTracks_region3a) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
                            indexXEnd_region3a(numTracks_region3a) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
                            indexXStart_region3a(numTracks_region3a) = dataTracks_rotated.entireEmbryo.indexXStart(i);
                            tracksX_region3a(:,numTracks_region3a) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
                            tracksY_region3a(:,numTracks_region3a) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
                        else
                            numTracks_outsideA = numTracks_outsideA +1;
                            
                        end
                        
                    elseif dataTracks_rotated.entireEmbryo.indexXStart(1,i) > ( furrow_position.image_start_detection - param.delta_furrowDetection_anaphaseOnset *param.sp6 ...
                            + param.delta_late_anaphase *param.sp6 ) % late_anaphase
                        
                        numTracks_late = numTracks_late +1;
                        lengthTracks_late(numTracks_late) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
                        indexXEnd_late(numTracks_late) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
                        indexXStart_late(numTracks_late) = dataTracks_rotated.entireEmbryo.indexXStart(i);
                        tracksX_late(:,numTracks_late) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
                        tracksY_late(:,numTracks_late) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
                        
                        if ( dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i) < ...
                                xEnd_3regions.region1(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) ) ...
                                && ( xStart_3regions.region1(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) <= ...
                                dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i)  )
                            numTracks_region1c = numTracks_region1c +1;
                            lengthTracks_region1c(numTracks_region1c) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
                            indexXEnd_region1c(numTracks_region1c) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
                            indexXStart_region1c(numTracks_region1c) = dataTracks_rotated.entireEmbryo.indexXStart(i);
                            tracksX_region1c(:,numTracks_region1c) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
                            tracksY_region1c(:,numTracks_region1c) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
                        elseif ( dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i) < ...
                                xEnd_3regions.region2(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) ) ...
                                && ( xStart_3regions.region2(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) <= ...
                                dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i)  )
                            numTracks_region2c = numTracks_region2c +1;
                            lengthTracks_region2c(numTracks_region2c) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
                            indexXEnd_region2c(numTracks_region2c) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
                            indexXStart_region2c(numTracks_region2c) = dataTracks_rotated.entireEmbryo.indexXStart(i);
                            tracksX_region2c(:,numTracks_region2c) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
                            tracksY_region2c(:,numTracks_region2c) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
                        elseif ( dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i) < ...
                                xEnd_3regions.region3(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) ) ...
                                && ( xStart_3regions.region3(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) <= ...
                                dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i)  )
                            numTracks_region3c = numTracks_region3c +1;
                            lengthTracks_region3c(numTracks_region3c) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
                            indexXEnd_region3c(numTracks_region3c) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
                            indexXStart_region3c(numTracks_region3c) = dataTracks_rotated.entireEmbryo.indexXStart(i);
                            tracksX_region3c(:,numTracks_region3c) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
                            tracksY_region3c(:,numTracks_region3c) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
                        else
                            numTracks_outsideC = numTracks_outsideC +1;
                            
                        end
                        
                    end
                else
                    if dataTracks_rotated.entireEmbryo.indexXStart(1,i) > ( furrow_position.image_start_detection - param.delta_furrowDetection_anaphaseOnset *param.sp6 )
                        
                        numTracks_anaphase = numTracks_anaphase +1;
                        lengthTracks_anaphase(numTracks_anaphase) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
                        indexXEnd_anaphase(numTracks_anaphase) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
                        indexXStart_anaphase(numTracks_anaphase) = dataTracks_rotated.entireEmbryo.indexXStart(i);
                        tracksX_anaphase(:,numTracks_anaphase) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
                        tracksY_anaphase(:,numTracks_anaphase) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
                        
                        if ( dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i) < ...
                                xEnd_3regions.region1(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) ) ...
                                && ( xStart_3regions.region1(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) <= ...
                                dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i)  )
                            numTracks_region1a = numTracks_region1a +1;
                            lengthTracks_region1a(numTracks_region1a) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
                            indexXEnd_region1a(numTracks_region1a) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
                            indexXStart_region1a(numTracks_region1a) = dataTracks_rotated.entireEmbryo.indexXStart(i);
                            tracksX_region1a(:,numTracks_region1a) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
                            tracksY_region1a(:,numTracks_region1a) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
                        elseif ( dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i) < ...
                                xEnd_3regions.region2(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) ) ...
                                && ( xStart_3regions.region2(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) <= ...
                                dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i)  )
                            numTracks_region2a = numTracks_region2a +1;
                            lengthTracks_region2a(numTracks_region2a) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
                            indexXEnd_region2a(numTracks_region2a) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
                            indexXStart_region2a(numTracks_region2a) = dataTracks_rotated.entireEmbryo.indexXStart(i);
                            tracksX_region2a(:,numTracks_region2a) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
                            tracksY_region2a(:,numTracks_region2a) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
                        elseif ( dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i) < ...
                                xEnd_3regions.region3(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) ) ...
                                && ( xStart_3regions.region3(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) <= ...
                                dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i)  )
                            numTracks_region3a = numTracks_region3a +1;
                            lengthTracks_region3a(numTracks_region3a) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
                            indexXEnd_region3a(numTracks_region3a) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
                            indexXStart_region3a(numTracks_region3a) = dataTracks_rotated.entireEmbryo.indexXStart(i);
                            tracksX_region3a(:,numTracks_region3a) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
                            tracksY_region3a(:,numTracks_region3a) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
                        else
                            numTracks_outsideA = numTracks_outsideA +1;
                            
                        end
                        
                        
                    end
                end
            end
        end
        
    end
end



%-----------------------------------------------------------------------------
% refrence time = pseudocleavage end

if time_reference_choice == 1
    
    for i = 1 : dataTracks_rotated.entireEmbryo.numTracks
        
        if dataTracks_rotated.entireEmbryo.indexXStart(1,i) <= ( (param.pseudoCleavage_endingTime + param.delta_NEBD_pseudoCleavage + param.delta_early_metaphase )*param.sp6 )
            
            if early_separation == 1  % to look at prophase / metaphase
                
                if dataTracks_rotated.entireEmbryo.indexXStart(1,i) <= ( (param.pseudoCleavage_endingTime + param.delta_NEBD_pseudoCleavage + param.delta_early_metaphase ) *param.sp6 ) ...
                        && dataTracks_rotated.entireEmbryo.indexXStart(1,i) > ( (param.pseudoCleavage_endingTime + param.delta_NEBD_pseudoCleavage )*param.sp6 ) % metaphase
                    
                    numTracks_metaphase = numTracks_metaphase +1;
                    lengthTracks_metaphase(numTracks_metaphase) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
                    indexXEnd_metaphase(numTracks_metaphase) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
                    indexXStart_metaphase(numTracks_metaphase) = dataTracks_rotated.entireEmbryo.indexXStart(i);
                    tracksX_metaphase(:,numTracks_metaphase) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
                    tracksY_metaphase(:,numTracks_metaphase) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
                    
                    if ( dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i) < ...
                            xEnd_3regions.region1(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) ) ...
                            && ( xStart_3regions.region1(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) <= ...
                            dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i)  )
                        numTracks_region1b = numTracks_region1b +1;
                        lengthTracks_region1b(numTracks_region1b) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
                        indexXEnd_region1b(numTracks_region1b) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
                        indexXStart_region1b(numTracks_region1b) = dataTracks_rotated.entireEmbryo.indexXStart(i);
                        tracksX_region1b(:,numTracks_region1b) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
                        tracksY_region1b(:,numTracks_region1b) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
                    elseif ( dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i) < ...
                            xEnd_3regions.region2(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) ) ...
                            && ( xStart_3regions.region2(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) <= ...
                            dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i)  )
                        numTracks_region2b = numTracks_region2b +1;
                        lengthTracks_region2b(numTracks_region2b) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
                        indexXEnd_region2b(numTracks_region2b) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
                        indexXStart_region2b(numTracks_region2b) = dataTracks_rotated.entireEmbryo.indexXStart(i);
                        tracksX_region2b(:,numTracks_region2b) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
                        tracksY_region2b(:,numTracks_region2b) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
                    elseif ( dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i) < ...
                            xEnd_3regions.region3(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) ) ...
                            && ( xStart_3regions.region3(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) <= ...
                            dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i)  )
                        numTracks_region3b = numTracks_region3b +1;
                        lengthTracks_region3b(numTracks_region3b) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
                        indexXEnd_region3b(numTracks_region3b) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
                        indexXStart_region3b(numTracks_region3b) = dataTracks_rotated.entireEmbryo.indexXStart(i);
                        tracksX_region3b(:,numTracks_region3b) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
                        tracksY_region3b(:,numTracks_region3b) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
                    else
                        numTracks_outsideB = numTracks_outsideB +1;
                        
                    end
                    
                    
                elseif dataTracks_rotated.entireEmbryo.indexXStart(1,i) <= ( (param.pseudoCleavage_endingTime + param.delta_NEBD_pseudoCleavage )*param.sp6 ) % prophase
                    
                    numTracks_prophase = numTracks_prophase +1;
                    lengthTracks_prophase(numTracks_prophase) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
                    indexXEnd_prophase(numTracks_prophase) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
                    indexXStart_prophase(numTracks_prophase) = dataTracks_rotated.entireEmbryo.indexXStart(i);
                    tracksX_prophase(:,numTracks_prophase) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
                    tracksY_prophase(:,numTracks_prophase) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
                    
                    if ( dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i) < ...
                            xEnd_3regions.region1(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) ) ...
                            && ( xStart_3regions.region1(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) <= ...
                            dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i)  )
                        numTracks_region1e = numTracks_region1e +1;
                        lengthTracks_region1e(numTracks_region1e) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
                        indexXEnd_region1e(numTracks_region1e) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
                        indexXStart_region1e(numTracks_region1e) = dataTracks_rotated.entireEmbryo.indexXStart(i);
                        tracksX_region1e(:,numTracks_region1e) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
                        tracksY_region1e(:,numTracks_region1e) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
                    elseif ( dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i) < ...
                            xEnd_3regions.region2(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) ) ...
                            && ( xStart_3regions.region2(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) <= ...
                            dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i)  )
                        numTracks_region2e = numTracks_region2e +1;
                        lengthTracks_region2e(numTracks_region2e) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
                        indexXEnd_region2e(numTracks_region2e) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
                        indexXStart_region2e(numTracks_region2e) = dataTracks_rotated.entireEmbryo.indexXStart(i);
                        tracksX_region2e(:,numTracks_region2e) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
                        tracksY_region2e(:,numTracks_region2e) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
                    elseif ( dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i) < ...
                            xEnd_3regions.region3(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) ) ...
                            && ( xStart_3regions.region3(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) <= ...
                            dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i)  )
                        numTracks_region3e = numTracks_region3e +1;
                        lengthTracks_region3e(numTracks_region3e) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
                        indexXEnd_region3e(numTracks_region3e) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
                        indexXStart_region3e(numTracks_region3e) = dataTracks_rotated.entireEmbryo.indexXStart(i);
                        tracksX_region3e(:,numTracks_region3e) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
                        tracksY_region3e(:,numTracks_region3e) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
                    else
                        numTracks_outsideE = numTracks_outsideE +1;
                        
                    end
                    
                end
                
            else %prometaphase
                
                numTracks_prometaphase = numTracks_prometaphase +1;
                lengthTracks_prometaphase(numTracks_prometaphase) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
                indexXEnd_prometaphase(numTracks_prometaphase) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
                indexXStart_prometaphase(numTracks_prometaphase) = dataTracks_rotated.entireEmbryo.indexXStart(i);
                tracksX_prometaphase(:,numTracks_prometaphase) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
                tracksY_prometaphase(:,numTracks_prometaphase) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
                
                if ( dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i) < ...
                        xEnd_3regions.region1(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) ) ...
                        && ( xStart_3regions.region1(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) <= ...
                        dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i)  )
                    numTracks_region1d = numTracks_region1d +1;
                    lengthTracks_region1d(numTracks_region1d) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
                    indexXEnd_region1d(numTracks_region1d) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
                    indexXStart_region1d(numTracks_region1d) = dataTracks_rotated.entireEmbryo.indexXStart(i);
                    tracksX_region1d(:,numTracks_region1d) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
                    tracksY_region1d(:,numTracks_region1d) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
                elseif ( dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i) < ...
                        xEnd_3regions.region2(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) ) ...
                        && ( xStart_3regions.region2(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) <= ...
                        dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i)  )
                    numTracks_region2d = numTracks_region2d +1;
                    lengthTracks_region2d(numTracks_region2d) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
                    indexXEnd_region2d(numTracks_region2d) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
                    indexXStart_region2d(numTracks_region2d) = dataTracks_rotated.entireEmbryo.indexXStart(i);
                    tracksX_region2d(:,numTracks_region2d) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
                    tracksY_region2d(:,numTracks_region2d) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
                elseif ( dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i) < ...
                        xEnd_3regions.region3(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) ) ...
                        && ( xStart_3regions.region3(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) <= ...
                        dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i)  )
                    numTracks_region3d = numTracks_region3d +1;
                    lengthTracks_region3d(numTracks_region3d) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
                    indexXEnd_region3d(numTracks_region3d) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
                    indexXStart_region3d(numTracks_region3d) = dataTracks_rotated.entireEmbryo.indexXStart(i);
                    tracksX_region3d(:,numTracks_region3d) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
                    tracksY_region3d(:,numTracks_region3d) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
                else
                    numTracks_outsideD = numTracks_outsideD +1;
                    
                end
                
            end
        else
            
            if late_separation == 1 % to look at anaphase / late mitosis
                
                if dataTracks_rotated.entireEmbryo.indexXStart(1,i) > ( (param.pseudoCleavage_endingTime + param.delta_NEBD_pseudoCleavage + param.delta_early_metaphase ) *param.sp6 ) ...
                        && dataTracks_rotated.entireEmbryo.indexXStart(1,i) <= ( (param.pseudoCleavage_endingTime + param.delta_NEBD_pseudoCleavage + param.delta_early_metaphase + ...
                        param.delta_late_anaphase ) *param.sp6 ) % anaphase
                    
                    numTracks_anaphase = numTracks_anaphase +1;
                    lengthTracks_anaphase(numTracks_anaphase) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
                    indexXEnd_anaphase(numTracks_anaphase) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
                    indexXStart_anaphase(numTracks_anaphase) = dataTracks_rotated.entireEmbryo.indexXStart(i);
                    tracksX_anaphase(:,numTracks_anaphase) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
                    tracksY_anaphase(:,numTracks_anaphase) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
                    
                    if ( dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i) < ...
                            xEnd_3regions.region1(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) ) ...
                            && ( xStart_3regions.region1(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) <= ...
                            dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i)  )
                        numTracks_region1a = numTracks_region1a +1;
                        lengthTracks_region1a(numTracks_region1a) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
                        indexXEnd_region1a(numTracks_region1a) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
                        indexXStart_region1a(numTracks_region1a) = dataTracks_rotated.entireEmbryo.indexXStart(i);
                        tracksX_region1a(:,numTracks_region1a) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
                        tracksY_region1a(:,numTracks_region1a) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
                    elseif ( dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i) < ...
                            xEnd_3regions.region2(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) ) ...
                            && ( xStart_3regions.region2(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) <= ...
                            dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i)  )
                        numTracks_region2a = numTracks_region2a +1;
                        lengthTracks_region2a(numTracks_region2a) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
                        indexXEnd_region2a(numTracks_region2a) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
                        indexXStart_region2a(numTracks_region2a) = dataTracks_rotated.entireEmbryo.indexXStart(i);
                        tracksX_region2a(:,numTracks_region2a) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
                        tracksY_region2a(:,numTracks_region2a) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
                    elseif ( dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i) < ...
                            xEnd_3regions.region3(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) ) ...
                            && ( xStart_3regions.region3(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) <= ...
                            dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i)  )
                        numTracks_region3a = numTracks_region3a +1;
                        lengthTracks_region3a(numTracks_region3a) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
                        indexXEnd_region3a(numTracks_region3a) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
                        indexXStart_region3a(numTracks_region3a) = dataTracks_rotated.entireEmbryo.indexXStart(i);
                        tracksX_region3a(:,numTracks_region3a) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
                        tracksY_region3a(:,numTracks_region3a) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
                    else
                        numTracks_outsideA = numTracks_outsideA +1;
                        
                    end
                    
                elseif dataTracks_rotated.entireEmbryo.indexXStart(1,i) > ( param.pseudoCleavage_endingTime + param.delta_NEBD_pseudoCleavage + param.delta_early_metaphase + ...
                        param.delta_late_anaphase ) *param.sp6 % late_anaphase
                    
                    numTracks_late = numTracks_late +1;
                    lengthTracks_late(numTracks_late) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
                    indexXEnd_late(numTracks_late) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
                    indexXStart_late(numTracks_late) = dataTracks_rotated.entireEmbryo.indexXStart(i);
                    tracksX_late(:,numTracks_late) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
                    tracksY_late(:,numTracks_late) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
                    
                    if ( dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i) < ...
                            xEnd_3regions.region1(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) ) ...
                            && ( xStart_3regions.region1(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) <= ...
                            dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i)  )
                        numTracks_region1c = numTracks_region1c +1;
                        lengthTracks_region1c(numTracks_region1c) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
                        indexXEnd_region1c(numTracks_region1c) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
                        indexXStart_region1c(numTracks_region1c) = dataTracks_rotated.entireEmbryo.indexXStart(i);
                        tracksX_region1c(:,numTracks_region1c) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
                        tracksY_region1c(:,numTracks_region1c) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
                    elseif ( dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i) < ...
                            xEnd_3regions.region2(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) ) ...
                            && ( xStart_3regions.region2(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) <= ...
                            dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i)  )
                        numTracks_region2c = numTracks_region2c +1;
                        lengthTracks_region2c(numTracks_region2c) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
                        indexXEnd_region2c(numTracks_region2c) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
                        indexXStart_region2c(numTracks_region2c) = dataTracks_rotated.entireEmbryo.indexXStart(i);
                        tracksX_region2c(:,numTracks_region2c) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
                        tracksY_region2c(:,numTracks_region2c) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
                    elseif ( dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i) < ...
                            xEnd_3regions.region3(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) ) ...
                            && ( xStart_3regions.region3(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) <= ...
                            dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i)  )
                        numTracks_region3c = numTracks_region3c +1;
                        lengthTracks_region3c(numTracks_region3c) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
                        indexXEnd_region3c(numTracks_region3c) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
                        indexXStart_region3c(numTracks_region3c) = dataTracks_rotated.entireEmbryo.indexXStart(i);
                        tracksX_region3c(:,numTracks_region3c) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
                        tracksY_region3c(:,numTracks_region3c) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
                    else
                        numTracks_outsideC = numTracks_outsideC +1;
                        
                    end
                    
                end
            else
                if dataTracks_rotated.entireEmbryo.indexXStart(1,i) > ( param.pseudoCleavage_endingTime + param.delta_NEBD_pseudoCleavage + param.delta_early_metaphase ) *param.sp6
                    
                    numTracks_anaphase = numTracks_anaphase +1;
                    lengthTracks_anaphase(numTracks_anaphase) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
                    indexXEnd_anaphase(numTracks_anaphase) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
                    indexXStart_anaphase(numTracks_anaphase) = dataTracks_rotated.entireEmbryo.indexXStart(i);
                    tracksX_anaphase(:,numTracks_anaphase) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
                    tracksY_anaphase(:,numTracks_anaphase) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
                    
                    if ( dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i) < ...
                            xEnd_3regions.region1(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) ) ...
                            && ( xStart_3regions.region1(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) <= ...
                            dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i)  )
                        numTracks_region1a = numTracks_region1a +1;
                        lengthTracks_region1a(numTracks_region1a) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
                        indexXEnd_region1a(numTracks_region1a) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
                        indexXStart_region1a(numTracks_region1a) = dataTracks_rotated.entireEmbryo.indexXStart(i);
                        tracksX_region1a(:,numTracks_region1a) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
                        tracksY_region1a(:,numTracks_region1a) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
                    elseif ( dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i) < ...
                            xEnd_3regions.region2(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) ) ...
                            && ( xStart_3regions.region2(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) <= ...
                            dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i)  )
                        numTracks_region2a = numTracks_region2a +1;
                        lengthTracks_region2a(numTracks_region2a) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
                        indexXEnd_region2a(numTracks_region2a) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
                        indexXStart_region2a(numTracks_region2a) = dataTracks_rotated.entireEmbryo.indexXStart(i);
                        tracksX_region2a(:,numTracks_region2a) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
                        tracksY_region2a(:,numTracks_region2a) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
                    elseif ( dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i) < ...
                            xEnd_3regions.region3(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) ) ...
                            && ( xStart_3regions.region3(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2) <= ...
                            dataTracks_rotated.entireEmbryo.tracksX(dataTracks_rotated.entireEmbryo.indexXStart(1,i)-param.sp2,i)  )
                        numTracks_region3a = numTracks_region3a +1;
                        lengthTracks_region3a(numTracks_region3a) = dataTracks_rotated.entireEmbryo.lengthTracks(i);
                        indexXEnd_region3a(numTracks_region3a) = dataTracks_rotated.entireEmbryo.indexXEnd(i);
                        indexXStart_region3a(numTracks_region3a) = dataTracks_rotated.entireEmbryo.indexXStart(i);
                        tracksX_region3a(:,numTracks_region3a) = dataTracks_rotated.entireEmbryo.tracksX(:,i);
                        tracksY_region3a(:,numTracks_region3a) = dataTracks_rotated.entireEmbryo.tracksY(:,i);
                    else
                        numTracks_outsideA = numTracks_outsideA +1;
                        
                    end
                                      
                end
            end
        end
    end
    
end
    
 
%------------------------------------------------------------------------- 
%save data from tracking in structure file

if numTracks_anaphase ~= 0
    dataTracks_rotated_anaphase.entireEmbryo.tracksX = tracksX_anaphase;
    dataTracks_rotated_anaphase.entireEmbryo.tracksY = tracksY_anaphase;
    dataTracks_rotated_anaphase.entireEmbryo.lengthTracks = lengthTracks_anaphase;
    dataTracks_rotated_anaphase.entireEmbryo.indexXStart = indexXStart_anaphase;
    dataTracks_rotated_anaphase.entireEmbryo.indexXEnd = indexXEnd_anaphase;
    dataTracks_rotated_anaphase.entireEmbryo.numTracks = numTracks_anaphase;
    
    clear tracksX_anaphase
    clear tracksY_anaphase
    clear lengthTracks_anaphase
    clear indexXStart_anaphase
    clear indexXEnd_anaphase
else
    dataTracks_rotated_anaphase.entireEmbryo.numTracks = 0;
end

if numTracks_metaphase ~= 0
    dataTracks_rotated_metaphase.entireEmbryo.tracksX = tracksX_metaphase;
    dataTracks_rotated_metaphase.entireEmbryo.tracksY = tracksY_metaphase;
    dataTracks_rotated_metaphase.entireEmbryo.lengthTracks = lengthTracks_metaphase;
    dataTracks_rotated_metaphase.entireEmbryo.indexXStart = indexXStart_metaphase;
    dataTracks_rotated_metaphase.entireEmbryo.indexXEnd = indexXEnd_metaphase;
    dataTracks_rotated_metaphase.entireEmbryo.numTracks = numTracks_metaphase;
    
    clear tracksX_metaphase
    clear tracksY_metaphase
    clear lengthTracks_metaphase
    clear indexXStart_metaphase
    clear indexXEnd_metaphase
else
    dataTracks_rotated_metaphase.entireEmbryo.numTracks = 0;
end

if numTracks_prometaphase ~= 0
    dataTracks_rotated_prometaphase.entireEmbryo.tracksX = tracksX_prometaphase;
    dataTracks_rotated_prometaphase.entireEmbryo.tracksY = tracksY_prometaphase;
    dataTracks_rotated_prometaphase.entireEmbryo.lengthTracks = lengthTracks_prometaphase;
    dataTracks_rotated_prometaphase.entireEmbryo.indexXStart = indexXStart_prometaphase;
    dataTracks_rotated_prometaphase.entireEmbryo.indexXEnd = indexXEnd_prometaphase;
    dataTracks_rotated_prometaphase.entireEmbryo.numTracks = numTracks_prometaphase;
    
    clear tracksX_prometaphase
    clear tracksY_prometaphase
    clear lengthTracks_prometaphase
    clear indexXStart_prometaphase
    clear indexXEnd_prometaphase
else
    dataTracks_rotated_prometaphase.entireEmbryo.numTracks = 0;
end

if numTracks_prophase ~= 0
    dataTracks_rotated_prophase.entireEmbryo.tracksX = tracksX_prophase;
    dataTracks_rotated_prophase.entireEmbryo.tracksY = tracksY_prophase;
    dataTracks_rotated_prophase.entireEmbryo.lengthTracks = lengthTracks_prophase;
    dataTracks_rotated_prophase.entireEmbryo.indexXStart = indexXStart_prophase;
    dataTracks_rotated_prophase.entireEmbryo.indexXEnd = indexXEnd_prophase;
    dataTracks_rotated_prophase.entireEmbryo.numTracks = numTracks_prophase;
    
    clear tracksX_prophase
    clear tracksY_prophase
    clear lengthTracks_prophase
    clear indexXStart_prophase
    clear indexXEnd_prophase
else
    dataTracks_rotated_prophase.entireEmbryo.numTracks = 0;
end

if numTracks_late ~= 0
    dataTracks_rotated_late.entireEmbryo.tracksX = tracksX_late;
    dataTracks_rotated_late.entireEmbryo.tracksY = tracksY_late;
    dataTracks_rotated_late.entireEmbryo.lengthTracks = lengthTracks_late;
    dataTracks_rotated_late.entireEmbryo.indexXStart = indexXStart_late;
    dataTracks_rotated_late.entireEmbryo.indexXEnd = indexXEnd_late;
    dataTracks_rotated_late.entireEmbryo.numTracks = numTracks_late;
    
    clear tracksX_late
    clear tracksY_late
    clear lengthTracks_late
    clear indexXStart_late
    clear indexXEnd_late
else
    dataTracks_rotated_late.entireEmbryo.numTracks = 0;
end


%------------------------------------

if numTracks_region1a ~= 0
    dataTracks_rotated_anaphase.region1.tracksX = tracksX_region1a;
    dataTracks_rotated_anaphase.region1.tracksY = tracksY_region1a;
    dataTracks_rotated_anaphase.region1.lengthTracks = lengthTracks_region1a;
    dataTracks_rotated_anaphase.region1.indexXStart = indexXStart_region1a;
    dataTracks_rotated_anaphase.region1.indexXEnd = indexXEnd_region1a;
    dataTracks_rotated_anaphase.region1.numTracks = numTracks_region1a;
    
    clear tracksX_region1a
    clear tracksY_region1a
    clear lengthTracks_region1a
    clear indexXStart_region1a
    clear indexXEnd_region1a
else
    dataTracks_rotated_anaphase.region1.numTracks = 0;
end


if numTracks_region2a ~= 0
    dataTracks_rotated_anaphase.region2.tracksX = tracksX_region2a;
    dataTracks_rotated_anaphase.region2.tracksY = tracksY_region2a;
    dataTracks_rotated_anaphase.region2.lengthTracks = lengthTracks_region2a;
    dataTracks_rotated_anaphase.region2.indexXStart = indexXStart_region2a;
    dataTracks_rotated_anaphase.region2.indexXEnd = indexXEnd_region2a;
    dataTracks_rotated_anaphase.region2.numTracks = numTracks_region2a;
    
    clear tracksX_region2a
    clear tracksY_region2a
    clear lengthTracks_region2a
    clear indexXStart_region2a
    clear indexXEnd_region2a
else
    dataTracks_rotated_anaphase.region2.numTracks = 0;
end

if numTracks_region3a ~= 0
    dataTracks_rotated_anaphase.region3.tracksX = tracksX_region3a;
    dataTracks_rotated_anaphase.region3.tracksY = tracksY_region3a;
    dataTracks_rotated_anaphase.region3.lengthTracks = lengthTracks_region3a;
    dataTracks_rotated_anaphase.region3.indexXStart = indexXStart_region3a;
    dataTracks_rotated_anaphase.region3.indexXEnd = indexXEnd_region3a;
    dataTracks_rotated_anaphase.region3.numTracks = numTracks_region3a;
    
    clear tracksX_region3a
    clear tracksY_region3a
    clear lengthTracks_region3a
    clear indexXStart_region3a
    clear indexXEnd_region3a
else
    dataTracks_rotated_anaphase.region3.numTracks = 0;
end

dataTracks_rotated_anaphase.numTracks_outside = numTracks_outsideA;

name2 = strcat('dataTracks_rotated_anaphase_3regions-', short_name, '.mat');
save(fullfile(pathMainDirectory,name2), '-struct', 'dataTracks_rotated_anaphase');

if numTracks_region1b ~= 0
    dataTracks_rotated_metaphase.region1.tracksX = tracksX_region1b;
    dataTracks_rotated_metaphase.region1.tracksY = tracksY_region1b;
    dataTracks_rotated_metaphase.region1.lengthTracks = lengthTracks_region1b;
    dataTracks_rotated_metaphase.region1.indexXStart = indexXStart_region1b;
    dataTracks_rotated_metaphase.region1.indexXEnd = indexXEnd_region1b;
    dataTracks_rotated_metaphase.region1.numTracks = numTracks_region1b;
    
    clear tracksX_region1b
    clear tracksY_region1b
    clear lengthTracks_region1b
    clear indexXStart_region1b
    clear indexXEnd_region1b
else
    dataTracks_rotated_metaphase.region1.numTracks = 0;
end

if numTracks_region2b ~= 0
    dataTracks_rotated_metaphase.region2.tracksX = tracksX_region2b;
    dataTracks_rotated_metaphase.region2.tracksY = tracksY_region2b;
    dataTracks_rotated_metaphase.region2.lengthTracks = lengthTracks_region2b;
    dataTracks_rotated_metaphase.region2.indexXStart = indexXStart_region2b;
    dataTracks_rotated_metaphase.region2.indexXEnd = indexXEnd_region2b;
    dataTracks_rotated_metaphase.region2.numTracks = numTracks_region2b;
    
    clear tracksX_region2b
    clear tracksY_region2b
    clear lengthTracks_region2b
    clear indexXStart_region2b
    clear indexXEnd_region2b
else
    dataTracks_rotated_metaphase.region2.numTracks = 0;
end

if numTracks_region3b ~= 0
    dataTracks_rotated_metaphase.region3.tracksX = tracksX_region3b;
    dataTracks_rotated_metaphase.region3.tracksY = tracksY_region3b;
    dataTracks_rotated_metaphase.region3.lengthTracks = lengthTracks_region3b;
    dataTracks_rotated_metaphase.region3.indexXStart = indexXStart_region3b;
    dataTracks_rotated_metaphase.region3.indexXEnd = indexXEnd_region3b;
    dataTracks_rotated_metaphase.region3.numTracks = numTracks_region3b;
    
    clear tracksX_region3b
    clear tracksY_region3b
    clear lengthTracks_region3b
    clear indexXStart_region3b
    clear indexXEnd_region3b
else
    dataTracks_rotated_metaphase.region3.numTracks = 0;
end

dataTracks_rotated_metaphase.numTracks_outside = numTracks_outsideB;

name2 = strcat('dataTracks_rotated_metaphase_3regions-', short_name, '.mat');
save(fullfile(pathMainDirectory,name2), '-struct', 'dataTracks_rotated_metaphase');

if numTracks_region1c ~= 0
    dataTracks_rotated_late.region1.tracksX = tracksX_region1c;
    dataTracks_rotated_late.region1.tracksY = tracksY_region1c;
    dataTracks_rotated_late.region1.lengthTracks = lengthTracks_region1c;
    dataTracks_rotated_late.region1.indexXStart = indexXStart_region1c;
    dataTracks_rotated_late.region1.indexXEnd = indexXEnd_region1c;
    dataTracks_rotated_late.region1.numTracks = numTracks_region1c;
    
    clear tracksX_region1c
    clear tracksY_region1c
    clear lengthTracks_region1c
    clear indexXStart_region1c
    clear indexXEnd_region1c
else
    dataTracks_rotated_late.region1.numTracks = 0;
end


if numTracks_region2c ~= 0
    dataTracks_rotated_late.region2.tracksX = tracksX_region2c;
    dataTracks_rotated_late.region2.tracksY = tracksY_region2c;
    dataTracks_rotated_late.region2.lengthTracks = lengthTracks_region2c;
    dataTracks_rotated_late.region2.indexXStart = indexXStart_region2c;
    dataTracks_rotated_late.region2.indexXEnd = indexXEnd_region2c;
    dataTracks_rotated_late.region2.numTracks = numTracks_region2c;
    
    clear tracksX_region2c
    clear tracksY_region2c
    clear lengthTracks_region2c
    clear indexXStart_region2c
    clear indexXEnd_region2c
else
    dataTracks_rotated_late.region2.numTracks = 0;
end

if numTracks_region3c ~= 0
    dataTracks_rotated_late.region3.tracksX = tracksX_region3c;
    dataTracks_rotated_late.region3.tracksY = tracksY_region3c;
    dataTracks_rotated_late.region3.lengthTracks = lengthTracks_region3c;
    dataTracks_rotated_late.region3.indexXStart = indexXStart_region3c;
    dataTracks_rotated_late.region3.indexXEnd = indexXEnd_region3c;
    dataTracks_rotated_late.region3.numTracks = numTracks_region3c;
    
    clear tracksX_region3c
    clear tracksY_region3c
    clear lengthTracks_region3c
    clear indexXStart_region3c
    clear indexXEnd_region3c
else
    dataTracks_rotated_late.region3.numTracks = 0;
end

dataTracks_rotated_late.numTracks_outside = numTracks_outsideA;

name2 = strcat('dataTracks_rotated_late_3regions-', short_name, '.mat');
save(fullfile(pathMainDirectory,name2), '-struct', 'dataTracks_rotated_late');



if numTracks_region1d ~= 0
    dataTracks_rotated_prometaphase.region1.tracksX = tracksX_region1d;
    dataTracks_rotated_prometaphase.region1.tracksY = tracksY_region1d;
    dataTracks_rotated_prometaphase.region1.lengthTracks = lengthTracks_region1d;
    dataTracks_rotated_prometaphase.region1.indexXStart = indexXStart_region1d;
    dataTracks_rotated_prometaphase.region1.indexXEnd = indexXEnd_region1d;
    dataTracks_rotated_prometaphase.region1.numTracks = numTracks_region1d;
    
    clear tracksX_region1d
    clear tracksY_region1d
    clear lengthTracks_region1d
    clear indexXStart_region1d
    clear indexXEnd_region1d
else
    dataTracks_rotated_prometaphase.region1.numTracks = 0;
end

if numTracks_region2d ~= 0
    dataTracks_rotated_prometaphase.region2.tracksX = tracksX_region2d;
    dataTracks_rotated_prometaphase.region2.tracksY = tracksY_region2d;
    dataTracks_rotated_prometaphase.region2.lengthTracks = lengthTracks_region2d;
    dataTracks_rotated_prometaphase.region2.indexXStart = indexXStart_region2d;
    dataTracks_rotated_prometaphase.region2.indexXEnd = indexXEnd_region2d;
    dataTracks_rotated_prometaphase.region2.numTracks = numTracks_region2d;
    
    clear tracksX_region2d
    clear tracksY_region2d
    clear lengthTracks_region2d
    clear indexXStart_region2d
    clear indexXEnd_region2d
else
    dataTracks_rotated_prometaphase.region2.numTracks = 0;
end

if numTracks_region3d ~= 0
    dataTracks_rotated_prometaphase.region3.tracksX = tracksX_region3d;
    dataTracks_rotated_prometaphase.region3.tracksY = tracksY_region3d;
    dataTracks_rotated_prometaphase.region3.lengthTracks = lengthTracks_region3d;
    dataTracks_rotated_prometaphase.region3.indexXStart = indexXStart_region3d;
    dataTracks_rotated_prometaphase.region3.indexXEnd = indexXEnd_region3d;
    dataTracks_rotated_prometaphase.region3.numTracks = numTracks_region3d;
    
    clear tracksX_region3d
    clear tracksY_region3d
    clear lengthTracks_region3d
    clear indexXStart_region3d
    clear indexXEnd_region3d
else
    dataTracks_rotated_prometaphase.region3.numTracks = 0;
end

dataTracks_rotated_prometaphase.numTracks_outside = numTracks_outsideB;

name2 = strcat('dataTracks_rotated_prometaphase_3regions-', short_name, '.mat');
save(fullfile(pathMainDirectory,name2), '-struct', 'dataTracks_rotated_prometaphase');


if numTracks_region1e ~= 0
    dataTracks_rotated_prophase.region1.tracksX = tracksX_region1e;
    dataTracks_rotated_prophase.region1.tracksY = tracksY_region1e;
    dataTracks_rotated_prophase.region1.lengthTracks = lengthTracks_region1e;
    dataTracks_rotated_prophase.region1.indexXStart = indexXStart_region1e;
    dataTracks_rotated_prophase.region1.indexXEnd = indexXEnd_region1e;
    dataTracks_rotated_prophase.region1.numTracks = numTracks_region1e;
    
    clear tracksX_region1e
    clear tracksY_region1e
    clear lengthTracks_region1e
    clear indexXStart_region1e
    clear indexXEnd_region1e
else
    dataTracks_rotated_prophase.region1.numTracks = 0;
end

if numTracks_region2e ~= 0
    dataTracks_rotated_prophase.region2.tracksX = tracksX_region2e;
    dataTracks_rotated_prophase.region2.tracksY = tracksY_region2e;
    dataTracks_rotated_prophase.region2.lengthTracks = lengthTracks_region2e;
    dataTracks_rotated_prophase.region2.indexXStart = indexXStart_region2e;
    dataTracks_rotated_prophase.region2.indexXEnd = indexXEnd_region2e;
    dataTracks_rotated_prophase.region2.numTracks = numTracks_region2e;
    
    clear tracksX_region2e
    clear tracksY_region2e
    clear lengthTracks_region2e
    clear indexXStart_region2e
    clear indexXEnd_region2e
else
    dataTracks_rotated_prophase.region2.numTracks = 0;
end

if numTracks_region3e ~= 0
    dataTracks_rotated_prophase.region3.tracksX = tracksX_region3e;
    dataTracks_rotated_prophase.region3.tracksY = tracksY_region3e;
    dataTracks_rotated_prophase.region3.lengthTracks = lengthTracks_region3e;
    dataTracks_rotated_prophase.region3.indexXStart = indexXStart_region3e;
    dataTracks_rotated_prophase.region3.indexXEnd = indexXEnd_region3e;
    dataTracks_rotated_prophase.region3.numTracks = numTracks_region3e;
    
    clear tracksX_region3e
    clear tracksY_region3e
    clear lengthTracks_region3e
    clear indexXStart_region3e
    clear indexXEnd_region3e
else
    dataTracks_rotated_prophase.region3.numTracks = 0;
end

dataTracks_rotated_prophase.numTracks_outside = numTracks_outsideE;

name2 = strcat('dataTracks_rotated_prophase_3regions-', short_name, '.mat');
save(fullfile(pathMainDirectory,name2), '-struct', 'dataTracks_rotated_prophase');



end

