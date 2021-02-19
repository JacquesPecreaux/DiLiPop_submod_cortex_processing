function [ dataTracks, dataTracks_afterFurrow, dataTracks_beforeFurrow, dataTracks_metaphase, dataTracks_anaphase, dataTracks_late, dataTracks_early ] = ...
    get_data_tracks( tracksFinal, furrow_position, embryo_name )

%load data obtained from tracking process and found interesting parameters
%for next tasks that will be performed

%all the trajectories data are saved in a structure in a mat file, whose name is
% dataTracks-nameRecording.mat


global pathMainDirectory;
global param;
global general_param

if nargin < 3
    embryo_name = short_name;
end
% if nargin < 2
%     furrow_position = '';
% end


%%
%get number of tracks and number of time points
if isstruct(tracksFinal) %if tracks are in structure format
    numTracks = length(tracksFinal);
    tmp = vertcat(tracksFinal.seqOfEvents);
    numTimePoints = max(tmp(:,1));
    clear tmp
else %if tracks are in matrix format
    [numTracks,numTimePoints] = size(tracksFinal);
    numTimePoints = numTimePoints/8;
end


%%
%get different infos from structure
if isstruct(tracksFinal) %if tracks are input in structure format
    
    %store the input structure as a variable with a different name
    inputStructure = tracksFinal;
    clear Track.tracksFinal;
    
    %get number of segments making each track
    numSegments = zeros(numTracks,1);
    for i = 1 : numTracks
        numSegments(i) = size(inputStructure(i).tracksCoordAmpCG,1);
    end
    
    %if all tracks have only one segment ...
    if max(numSegments) == 1
        
        %indicate that there are no compound tracks with merging and splitting branches
        mergeSplit = 0;
        
        %locate the row of the first track of each compound track in the
        %big matrix of all tracks (to be constructed in the next step)
        %in this case of course every compound track is simply one track
        %without branches
        trackStartRow = (1:numTracks)';
        %   rowcol=contorno(image);
        
        % L = imerode(L,seD);
        %store tracks in a matrix
        tracksFinal = NaN*ones(numTracks,8*numTimePoints);
        for i = 1 : numTracks
            startTime = inputStructure(i).seqOfEvents(1,1);
            endTime   = inputStructure(i).seqOfEvents(end,1);
            tracksFinal(i,8*(startTime-1)+1:8*endTime) = inputStructure(i).tracksCoordAmpCG;
        end
        
    else %if some tracks have merging/splitting branches
        
        %indicate that in the variable mergeSplit
        mergeSplit = 1;
        
        %locate the row of the first track of each compound track in the
        %big matrix of all tracks (to be constructed in the next step)
        trackStartRow = ones(numTracks,1);
        for iTrack = 2 : numTracks
            trackStartRow(iTrack) = trackStartRow(iTrack-1) + numSegments(iTrack-1);
        end
        
        %put all tracks together in a matrix
        tracksFinal = NaN*ones(trackStartRow(end)+numSegments(end)-1,8*numTimePoints);
        for i = 1 : numTracks
            startTime = inputStructure(i).seqOfEvents(1,1);
            endTime   = inputStructure(i).seqOfEvents(end,1);
            tracksFinal(trackStartRow(i):trackStartRow(i)+...
                numSegments(i)-1,8*(startTime-1)+1:8*endTime) = ...
                inputStructure(i).tracksCoordAmpCG;
        end
        
    end
    
else %if tracks are not input in structure format
    
    %indicate that there are no compound tracks with merging and splitting branches
    mergeSplit = 0;
    
    %indicate that each track consists of one segment
    numSegments = ones(numTracks,1);
    
    %locate the row of the first track of each compound track in the
    %big matrix of all tracks
    %in this case of course every compound track is simply one track
    %without branches
    trackStartRow = (1:numTracks)';
    
end


%%
%get the x,y-coordinates of features in all tracks

tracksX = tracksFinal(:,1:8:end)' ;
tracksY = tracksFinal(:,2:8:end)' ;


%%
%length of each track, index of start and index of end
for i = 1 : trackStartRow(end) + numSegments(end) - 1
    obsAvail = find(~isnan(tracksX(:,i)));
    lengthTracks(i) = length((obsAvail(1):obsAvail(end)));
    indexXStart(i) = obsAvail(1);
    indexXEnd(i) = obsAvail(end);
end

if size(lengthTracks,2) >100
    residencyTime_mean = mean(lengthTracks)*param.decimate/param.sp6;
else
    residencyTime_mean = NaN;
end


%%
%save data from tracking in structure file
dataTracks.entireEmbryo.tracksX = tracksX;
dataTracks.entireEmbryo.tracksY = tracksY;
dataTracks.entireEmbryo.lengthTracks = lengthTracks;
dataTracks.entireEmbryo.indexXStart = indexXStart+param.sp2;
dataTracks.entireEmbryo.indexXEnd = indexXEnd+param.sp2;
dataTracks.entireEmbryo.numTracks = numTracks;
dataTracks.entireEmbryo.residencyTime = residencyTime_mean;
dataTracks.numTimePoints = numTimePoints;

name = strcat('dataTracks-', embryo_name, param.extra, '.mat');
save(fullfile(pathMainDirectory,name), '-struct', 'dataTracks');

clear tracksX
clear tracksY
clear lengthTracks
clear indexXStart
clear indexXEnd
clear numTracks

%% set number of tracks in different regions to NaN
% necessary since check on numTracks to perform or not analysis on tracks
% ex: fitting duration distribution, study motion, ...

dataTracks_beforeFurrow.entireEmbryo.numTracks = NaN;
dataTracks_afterFurrow.entireEmbryo.numTracks = NaN;
dataTracks_metaphase.entireEmbryo.numTracks = NaN;
dataTracks_anaphase.entireEmbryo.numTracks = NaN;
dataTracks_late.entireEmbryo.numTracks = NaN;
dataTracks_early.entireEmbryo.numTracks = NaN;


%% separate tracks before/after furrow detection


if exist('furrow_position','var') == 1 && general_param.cortex_analysis.minimal_analysis == 0
    % after apparition furrow at cortex
    
    numTracks_after = 0;
    numTracks_before = 0;
    
    if ~isnan(furrow_position.image_start_detection) && furrow_position.image_start_detection < param.sp3
        
        for i = 1 : dataTracks.entireEmbryo.numTracks
            if dataTracks.entireEmbryo.indexXStart(1,i) > furrow_position.image_start_detection
                numTracks_after = numTracks_after +1;
                lengthTracks_after(numTracks_after) = dataTracks.entireEmbryo.lengthTracks(i);
                indexXEnd_after(numTracks_after) = dataTracks.entireEmbryo.indexXEnd(i);
                indexXStart_after(numTracks_after) = dataTracks.entireEmbryo.indexXStart(i);
                tracksX_after(:,numTracks_after) = dataTracks.entireEmbryo.tracksX(:,i);
                tracksY_after(:,numTracks_after) = dataTracks.entireEmbryo.tracksY(:,i);
            else
                numTracks_before = numTracks_before +1;
                lengthTracks_before(numTracks_before) = dataTracks.entireEmbryo.lengthTracks(i);
                indexXEnd_before(numTracks_before) = dataTracks.entireEmbryo.indexXEnd(i);
                indexXStart_before(numTracks_before) = dataTracks.entireEmbryo.indexXStart(i);
                tracksX_before(:,numTracks_before) = dataTracks.entireEmbryo.tracksX(:,i);
                tracksY_before(:,numTracks_before) = dataTracks.entireEmbryo.tracksY(:,i);
            end
        end
        
        %save data from tracking in structure file
        if numTracks_after > 0
            dataTracks_afterFurrow.entireEmbryo.tracksX = tracksX_after;
            dataTracks_afterFurrow.entireEmbryo.tracksY = tracksY_after;
            dataTracks_afterFurrow.entireEmbryo.lengthTracks = lengthTracks_after;
            dataTracks_afterFurrow.entireEmbryo.indexXStart = indexXStart_after;
            dataTracks_afterFurrow.entireEmbryo.indexXEnd = indexXEnd_after;
            dataTracks_afterFurrow.entireEmbryo.numTracks = numTracks_after;
            dataTracks_afterFurrow.numTimePoints = numTimePoints;
        end
        
        
        name2 = strcat('dataTracks_afterFurrowDetection-', embryo_name, param.extra, '.mat');
        save(fullfile(pathMainDirectory,name2), '-struct', 'dataTracks_afterFurrow');
        
        if numTracks_before > 0
            dataTracks_beforeFurrow.entireEmbryo.tracksX = tracksX_before;
            dataTracks_beforeFurrow.entireEmbryo.tracksY = tracksY_before;
            dataTracks_beforeFurrow.entireEmbryo.lengthTracks = lengthTracks_before;
            dataTracks_beforeFurrow.entireEmbryo.indexXStart = indexXStart_before;
            dataTracks_beforeFurrow.entireEmbryo.indexXEnd = indexXEnd_before;
            dataTracks_beforeFurrow.entireEmbryo.numTracks = numTracks_before;
            dataTracks_beforeFurrow.numTimePoints = numTimePoints;
        end
        
        name2 = strcat('dataTracks_beforeFurrowDetection-', embryo_name, param.extra, '.mat');
        save(fullfile(pathMainDirectory,name2), '-struct', 'dataTracks_beforeFurrow');
        
        clear tracksX_before
        clear tracksY_before
        clear lengthTracks_before
        clear indexXStart_before
        clear indexXEnd_before
        
        clear tracksX_after
        clear tracksY_after
        clear lengthTracks_after
        clear indexXStart_after
        clear indexXEnd_after
        
        numTracks_after = 0;
        numTracks_before = 0;
        numTracks_late = 0;
        numTracks_early = 0;
        
        if ( furrow_position.image_start_detection - param.delta_furrowDetection_anaphaseOnset *param.sp6 ) > 100 && ...
                ( furrow_position.image_start_detection - param.delta_furrowDetection_anaphaseOnset *param.sp6 > param.sp2 ) && ...
                furrow_position.image_start_detection < param.sp3
            
            %    numTracks_after = 0;
            %    numTracks_before = 0;
            
            for i = 1 : dataTracks.entireEmbryo.numTracks
                if dataTracks.entireEmbryo.indexXStart(1,i) > ( furrow_position.image_start_detection - param.delta_furrowDetection_anaphaseOnset *param.sp6 ) ...
                        && dataTracks.entireEmbryo.indexXStart(1,i) <= ( furrow_position.image_start_detection - param.delta_furrowDetection_anaphaseOnset *param.sp6 ...
                        + param.delta_late_anaphase *param.sp6 ) % anaphse
                    numTracks_after = numTracks_after +1;
                    lengthTracks_after(numTracks_after) = dataTracks.entireEmbryo.lengthTracks(i);
                    indexXEnd_after(numTracks_after) = dataTracks.entireEmbryo.indexXEnd(i);
                    indexXStart_after(numTracks_after) = dataTracks.entireEmbryo.indexXStart(i);
                    tracksX_after(:,numTracks_after) = dataTracks.entireEmbryo.tracksX(:,i);
                    tracksY_after(:,numTracks_after) = dataTracks.entireEmbryo.tracksY(:,i);
                    
                elseif dataTracks.entireEmbryo.indexXStart(1,i) <= ( furrow_position.image_start_detection - param.delta_furrowDetection_anaphaseOnset *param.sp6 ) 
                    if dataTracks.entireEmbryo.indexXStart(1,i) >= ( furrow_position.image_start_detection - param.delta_furrowDetection_anaphaseOnset *param.sp6 ...
                            - param.delta_early_metaphase *param.sp6 ) % metaphase
                        numTracks_before = numTracks_before +1;
                        lengthTracks_before(numTracks_before) = dataTracks.entireEmbryo.lengthTracks(i);
                        indexXEnd_before(numTracks_before) = dataTracks.entireEmbryo.indexXEnd(i);
                        indexXStart_before(numTracks_before) = dataTracks.entireEmbryo.indexXStart(i);
                        tracksX_before(:,numTracks_before) = dataTracks.entireEmbryo.tracksX(:,i);
                        tracksY_before(:,numTracks_before) = dataTracks.entireEmbryo.tracksY(:,i);
                    elseif dataTracks.entireEmbryo.indexXStart(1,i) < ( furrow_position.image_start_detection - param.delta_furrowDetection_anaphaseOnset *param.sp6 ...
                            - param.delta_early_metaphase *param.sp6 ) % prophase
                        numTracks_early = numTracks_early +1;
                        lengthTracks_early(numTracks_early) = dataTracks.entireEmbryo.lengthTracks(i);
                        indexXEnd_early(numTracks_early) = dataTracks.entireEmbryo.indexXEnd(i);
                        indexXStart_early(numTracks_early) = dataTracks.entireEmbryo.indexXStart(i);
                        tracksX_early(:,numTracks_early) = dataTracks.entireEmbryo.tracksX(:,i);
                        tracksY_early(:,numTracks_early) = dataTracks.entireEmbryo.tracksY(:,i);
                    end
                    
                elseif dataTracks.entireEmbryo.indexXStart(1,i) > ( furrow_position.image_start_detection - param.delta_furrowDetection_anaphaseOnset *param.sp6 ...
                        + param.delta_late_anaphase *param.sp6 ) % late anaphase
                    numTracks_late = numTracks_late +1;
                    lengthTracks_late(numTracks_late) = dataTracks.entireEmbryo.lengthTracks(i);
                    indexXEnd_late(numTracks_late) = dataTracks.entireEmbryo.indexXEnd(i);
                    indexXStart_late(numTracks_late) = dataTracks.entireEmbryo.indexXStart(i);
                    tracksX_late(:,numTracks_late) = dataTracks.entireEmbryo.tracksX(:,i);
                    tracksY_late(:,numTracks_late) = dataTracks.entireEmbryo.tracksY(:,i);                    
                end
            end
            
            %save data from tracking in structure file
            if numTracks_after > 0
                dataTracks_anaphase.entireEmbryo.tracksX = tracksX_after;
                dataTracks_anaphase.entireEmbryo.tracksY = tracksY_after;
                dataTracks_anaphase.entireEmbryo.lengthTracks = lengthTracks_after;
                dataTracks_anaphase.entireEmbryo.indexXStart = indexXStart_after;
                dataTracks_anaphase.entireEmbryo.indexXEnd = indexXEnd_after;
                dataTracks_anaphase.entireEmbryo.numTracks = numTracks_after;
                dataTracks_anaphase.numTimePoints = numTimePoints;
            end
            
            
            name2 = strcat('dataTracks_anaphase-', embryo_name, param.extra, '.mat');
            save(fullfile(pathMainDirectory,name2), '-struct', 'dataTracks_anaphase');
            
            if numTracks_before > 0
                dataTracks_metaphase.entireEmbryo.tracksX = tracksX_before;
                dataTracks_metaphase.entireEmbryo.tracksY = tracksY_before;
                dataTracks_metaphase.entireEmbryo.lengthTracks = lengthTracks_before;
                dataTracks_metaphase.entireEmbryo.indexXStart = indexXStart_before;
                dataTracks_metaphase.entireEmbryo.indexXEnd = indexXEnd_before;
                dataTracks_metaphase.entireEmbryo.numTracks = numTracks_before;
                dataTracks_metaphase.numTimePoints = numTimePoints;
            end
            
            name2 = strcat('dataTracks_metaphase-', embryo_name, param.extra, '.mat');
            save(fullfile(pathMainDirectory,name2), '-struct', 'dataTracks_metaphase');

            if numTracks_early > 0
                dataTracks_early.entireEmbryo.tracksX = tracksX_early;
                dataTracks_early.entireEmbryo.tracksY = tracksY_early;
                dataTracks_early.entireEmbryo.lengthTracks = lengthTracks_early;
                dataTracks_early.entireEmbryo.indexXStart = indexXStart_early;
                dataTracks_early.entireEmbryo.indexXEnd = indexXEnd_early;
                dataTracks_early.entireEmbryo.numTracks = numTracks_early;
                dataTracks_early.numTimePoints = numTimePoints;
            end
            
            name2 = strcat('dataTracks_early-', embryo_name, param.extra, '.mat');
            save(fullfile(pathMainDirectory,name2), '-struct', 'dataTracks_early');            
            
            if numTracks_late > 0
                dataTracks_late.entireEmbryo.tracksX = tracksX_late;
                dataTracks_late.entireEmbryo.tracksY = tracksY_late;
                dataTracks_late.entireEmbryo.lengthTracks = lengthTracks_late;
                dataTracks_late.entireEmbryo.indexXStart = indexXStart_late;
                dataTracks_late.entireEmbryo.indexXEnd = indexXEnd_late;
                dataTracks_late.entireEmbryo.numTracks = numTracks_late;
                dataTracks_late.numTimePoints = numTimePoints;
            end
            
            name2 = strcat('dataTracks_late-', embryo_name, param.extra, '.mat');
            save(fullfile(pathMainDirectory,name2), '-struct', 'dataTracks_late');            
            
            clear tracksX_before
            clear tracksY_before
            clear lengthTracks_before
            clear indexXStart_before
            clear indexXEnd_before
            
            clear tracksX_early
            clear tracksY_early
            clear lengthTracks_early
            clear indexXStart_early
            clear indexXEnd_early            
            
            clear tracksX_after
            clear tracksY_after
            clear lengthTracks_after
            clear indexXStart_after
            clear indexXEnd_after
            
            clear tracksX_late
            clear tracksY_late
            clear lengthTracks_late
            clear indexXStart_late
            clear indexXEnd_late            
            
        end
        
    end
    
end

clear numTimePoints

end

