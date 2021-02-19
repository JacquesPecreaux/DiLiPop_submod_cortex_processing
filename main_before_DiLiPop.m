function main_before_DiLiPop(xmlfile,tag)

% This main script has to be done before DiLiPop analysis to generate files
% that will be used to get DiLiPop parameters and maps.

% 1/ Download of the raw stack
% 2/ Tracking of embryo contour along time: this will generate the following mat files:
%   - regionArea
%   - regionXlength
%   - regionXlimit
% and if a furrow ingression is present and investigated (general_param.furrow_detection.status_performed = 1  )
%   - area_blastomere
%   - furrow_position_convexity

global param
global general_param
global pathMainDirectory

%% input args
% 1st arg: path of a xmlfile (mat file) for a given condition that contains the parameters associlated
% 2nd arg: possibility to add a tag in the name of the file generated to disciminate them from existing files, if
% you wish to perform analysis of the same dataset, with different setting.

if nargin<1
    [xmlfile, p] = uigetfile('*.mat','Please choose a job file to process');
    [~,xmlfile_bkp,~] = fileparts(xmlfile);
    xmlfile = fullfile(p, xmlfile);
end

if nargin<2
    tag = input_perso([' Please, give a tag to embryo (empty = no tag): '], ''  );
end

load(xmlfile); % load structure containing general_param and param
general_param = saveVarsMat_new.general_param;% read general_param from structure


%% for each embryo

for k = 1:length(saveVarsMat_new.params) % all embryo
    
    param = saveVarsMat_new.params{k};
    
    if param.status >= 0
        
        disp( param.sp1);
        
        param.extra = tag;
        %________________________________________________________________________________________________
        % DOWNLOAD RAW STACK
        
        % format = 'mtif' (tif stack) ;padding ( = param.sp7) = -1; param_set = '';
        path_tiff_stack = strcat(param.basepath , '/' , param.sp1 , '/');
        name_tiff_stack = param.stem_name;
        name = fullfile(path_tiff_stack,name_tiff_stack);
        
        [end_reading,siz,fitsnom] = read_init(name,param.format_image,param.sp7,param.sp2,param.sp3,'');
        
        % possibility to skeep some image from the stack with decimate, as
        % contour takes time
        number_images = (param.sp3 - param.sp2) +1;
        
        if param.landing_analysis == 0 % will download only images for contour detection
            if param.cortex_pass2.decimate > 0
                loop_array=0:param.cortex_pass2.decimate:(number_images-1);
            else
                loop_array=(number_images-1):param.cortex_pass2.decimate:0;
            end          
            % loop
            for index_ = 1 : 1 : length(loop_array)
                [Image,~,error_reading]=read_with_trial(param.sp2+loop_array(index_),param.sp7,param.format_image,siz,fitsnom,'none',param.sp3,...
                    param.cortex_pass2.channel_interest_AC,'');
                if error_reading || isempty(Image)
                    error('JACQ:FAILREAD','fail to read the file');
                end
                if index_ == 1
                    Imagee_ref = Image;
                end
                Image2 = mat2gray(Image,stretchlim(Imagee_ref,0));
                imageStack_raw (:,:,index_) = single(Image2);
                clear Image2 Image
            end
            
        elseif param.landing_analysis == 1 % will download all images
            % loop
            for index_ = 0 : 1 : number_images-1
                [Image,~,error_reading]=read_with_trial(param.sp2+index_,param.sp7,param.format_image,siz,fitsnom,'none',param.sp3,...
                    param.cortex_pass2.channel_interest_AC,'');
                if error_reading || isempty(Image)
                    error('JACQ:FAILREAD','fail to read the file');
                end
                if index_ == 0
                    Imagee_ref = Image;
                end
                Image2 = mat2gray(Image,stretchlim(Imagee_ref,0));
                imageStack_raw (:,:,index_+1) = single(Image2);
                clear Image2 Image
            end                    
        end
        disp('Raw stack downloaded');
        
        %_____________________________________________________________________________________________________________
        % load the binary mask generated using the function to_perform_DiLiPop_ROI.m

        mask_BW_AC = logical(imread(param.cortex_pass1.mask_image));

               
        %______________________________________________________________________________________________________
        % perform the contour detection along time and generate mask needed for coming quantifications
        
        mainDirectory = 'contour detection';
        mkdir(path_tiff_stack,mainDirectory);
        pathMainDirectory = strcat(path_tiff_stack , '/' , mainDirectory, '/');
        
        %--------------------------
        % save image to check correct images
        fig_firstImage = figure;
        imshow(imageStack_raw(:,:,1));
        figureName = strcat('firstImage-', short_name,  param.extra, '.tif');
        saveas(gcf,fullfile(pathMainDirectory,figureName));
        close(fig_firstImage)
        
        fig_lastImage = figure;
        imshow(imageStack_raw(:,:,end));
        figureName = strcat('lastImage-', short_name,  param.extra, '.tif');
        saveas(gcf,fullfile(pathMainDirectory,figureName));
        close(fig_lastImage)
        
        fig_mask_BW = figure;
        imshow(mask_BW_AC);
        figureName = strcat('mask_BW_AC-', short_name,  param.extra, '.tif');
        saveas(gcf,fullfile(pathMainDirectory,figureName));
        close(fig_mask_BW)
        
        %-------------------------
        % to get active contours of the embryo
        [~,segmentation] = perform_active_contours(imageStack_raw);
        disp('active contours obtained');
        
        %---------------------------
        %convert AC into mask
        rectangular_mask = get_rectangular_mask(mask_BW_AC,Imagee_ref);
        [ maskedStack] = from_contour_to_mask(segmentation,imageStack_raw,rectangular_mask);
        disp('AC converted into mask');
        clear segmentation
        clear Imagee_ref
        
        %--------------------------
        % save image to check correct images
        fig_firstMask = figure;
        imshow(maskedStack(:,:,1));
        figureName = strcat('firstMask-', short_name,  param.extra, '.tif');
        saveas(gcf,fullfile(pathMainDirectory,figureName));
        close(fig_firstMask)
        
        fig_lastMask = figure;
        imshow(maskedStack(:,:,end));
        figureName = strcat('lastMask-', short_name,  param.extra, '.tif');
        saveas(gcf,fullfile(pathMainDirectory,figureName));
        close(fig_lastMask)
        
        %----------------------
        % generate stacks needed to quantify embryo area/limits/lengths
        
        % crop stack according to ROI polygonal selection
        for i = 1 : size(imageStack_raw,3)
            imageStack(:,:,i) = imcrop(imageStack_raw(:,:,i),rectangular_mask);
        end
        clear imageStack_raw
        
        % generate and save image masked stack that could be used by user
        % to track fluorescent spots, masking preventing to consider signal
        % outisde from the embryo
        maskedImageStack = imageStack.*maskedStack;
        name_stack = strcat(param.stem_name,'_imageMaskedStack.tif');
        name_path_imageMaskedStack = fullfile(path_tiff_stack,name_stack);
        write_tif_stack(name_path_imageMaskedStack,double(maskedImageStack));
        clear maskedImageStack
   
        % generate stack use
        maskedStack_NaN = single(maskedStack);
        maskedStack_NaN(maskedStack_NaN==0)=NaN;
        imageMaskedStack_NaN = single(imageStack.* maskedStack_NaN);
        clear maskedStack_NaN
        clear imageStack
        
        %--------------------------
        % save image to check correct images
        fig_firstMaskedImage = figure;
        imshow(imageMaskedStack_NaN(:,:,1));
        figureName = strcat('firstMaskedImage-', short_name,  param.extra, '.tif');
        saveas(gcf,fullfile(pathMainDirectory,figureName));
        close(fig_firstMaskedImage)
        
        fig_lastMaskedImage = figure;
        imshow(imageMaskedStack_NaN(:,:,end));
        figureName = strcat('lastMaskedImage-', short_name,  param.extra, '.tif');
        saveas(gcf,fullfile(pathMainDirectory,figureName));
        close(fig_lastMaskedImage)
        
        
        %_________________________________________________________________________________________________
        % caclulate the area of the embryo, the limits over different region width
        % NB: area in pixel squared
        
        [ regionXlimit,regionXlength,regionArea] = calculate_area_and_limits( imageMaskedStack_NaN,maskedStack );
        clear imageMaskedStack_NaN
        disp('area and limits of embryo calculated');
        
        
        %____________________________________________________________________________________________________
        % study the furrow ingression: position and timing
        
        if general_param.furrow_detection.status_perform == 1
            
            mainDirectory = 'furrow_characterization';
            mkdir_perso(path_tiff_stack,mainDirectory);
            pathMainDirectory = strcat(path_tiff_stack , '/' , mainDirectory, '/');
            
            %------------------
            % characterization of furrow
            [ furrow_onset_sec,furrow_position ] = characterize_furrow(name,regionXlimit,regionXlength,maskedStack,pathMainDirectory); % get furrow onset time in sec
            disp('characterization of furrow done');
            param.furrow_detection_time = furrow_onset_sec;
            %-------------------------------
            % determine area of each blastomere
            if ~isnan( sum(furrow_position.xCoordinate.timeDependence(:)) )
                [maskedStack_anterior, ~,~] ...
                    = mask_anterior_posterior_side(maskedStack,furrow_position,regionXlimit,regionXlength);
                dbMaskedStack_anterior = logical(maskedStack_anterior.* maskedStack);
                clear maskedStack_anterior
                to_get_area_blastomere( dbMaskedStack_anterior,regionArea,pathMainDirectory );
            end
            disp('blastomere areas obtained');
            clear dbMaskedStack_anterior 
        end
        clear maskedStack regionXlimit regionXlength regionArea furrow_position
        
        
        %___________________________________________________________________________________________________
        % create a txt file with data about contour detection
        
        mainDirectory = 'contour detection';
        mkdir(path_tiff_stack,mainDirectory);
        pathMainDirectory = strcat(path_tiff_stack , '/' , mainDirectory, '/');
        
        id = param.stem_name;
        contour_detection_name=[sprintf('%s%s_%s',pathMainDirectory,id,tag) '_contourDetection.txt'];
        
        fp1=fopen(contour_detection_name,'w');
        fprintf(fp1,'Active contour parameter \n');
        fprintf(fp1,'Cortex pass1 \n');
        fprintf(fp1,'kappa: \t%e\t%s\t%f\t%s\t%f\n',param.cortex_pass1.kap, 'Lin', param.cortex_pass1.Lin, ...
            'Lout', param.cortex_pass1.Lout);
        fprintf(fp1,'clahe: \t%f\t%s\t%f\n',param.cortex_pass1.clahe_image,'cliplimit: ',param.cortex_pass1.clahe_cliplimit);
        fprintf(fp1,'filter: \t%f\t%s\t%s\n',param.cortex_pass1.filter_image,'kind: ',param.cortex_pass1.kind_filter);
        fprintf(fp1,'kalman: \t%f\n',param.cortex_pass1.kalman_image);
        fprintf(fp1,'recompute_c1_c2_mode: \t%f\t%s\t%f\n', param.cortex_pass1.recompute_c1_c2_mode, ...
            'reuse_initial_contour: ',param.cortex_pass1.reuse_initial_contour);
        fprintf(fp1,'mask: initialisation: \t%f\t%s\t%f\n',param.cortex_pass1.initialization_mask,...
            'maskFromInit: ' , param.cortex_pass1.maskFromInitialization_image);
        fprintf(fp1,'Cortex pass2 \n');
        fprintf(fp1,'kappa: \t%e\t%s\t%f\t%s\t%f\n',param.cortex_pass2.kap, 'Lin', param.cortex_pass2.Lin, ...
            'Lout', param.cortex_pass1.Lout);
        fprintf(fp1,'clahe: \t%f\t%s\t%f\n',param.cortex_pass2.clahe_image,'cliplimit: ',param.cortex_pass2.clahe_cliplimit);
        fprintf(fp1,'filter: \t%f\t%s\t%s\n',param.cortex_pass2.filter_image,'kind: ',param.cortex_pass2.kind_filter);
        fprintf(fp1,'kalman: \t%f\n',param.cortex_pass2.kalman_image);
        fprintf(fp1,'recompute_c1_c2_mode: \t%f\t%s\t%f\n', param.cortex_pass2.recompute_c1_c2_mode, ...
            'reuse_initial_contour: ',param.cortex_pass2.reuse_initial_contour);
        fprintf(fp1,'mask: initialisation: \t%f\t%s\t%f\n',param.cortex_pass2.initialization_mask,...
            'maskFromInit: ', param.cortex_pass2.maskFromInitialization_image);
        fprintf(fp1,'\n');
        fclose(fp1);
        
        clear fp1
        clear contour_detection_name
        clear rectangular_mask
        clear name_
        clear id
        
        clear_global_variables

        saveVarsMat_new.params{k} = param ;
        
        %______________________________________________________________________________________________________
        % possibility to edit the script below to perform Kalman denoising
        % or utrack analysis

        
        
        
    end
end

new_name = [xmlfile_bkp, '_updated' ];
save(fullfile(p, new_name), 'saveVarsMat_new');
close all
clear all

end

