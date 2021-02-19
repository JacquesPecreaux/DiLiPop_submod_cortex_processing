function to_perform_DiLiPop_ROI(xml_matfile)

% function that will enable to draw a Rregion of interest around the
% studied emrbyo to exclude bright structure that will impair the contiur
% detection.
% For each embryo of the dataset, the user will be ask to draw the roi.

% 1st arg: path of a xmlfile (mat file) for a given condition that contains the parameters associlated

global param
global general_param

%% download xmlfile and load param and general_param

if nargin < 1
    [xmlfile, p] = uigetfile('*.mat','Please choose a job file to process');
    [~,xmlfile_bkp,~] = fileparts(xmlfile);
    xml_matfile = fullfile(p, xmlfile);
end

load(xml_matfile); % load structure containing general_param and param
general_param = saveVarsMat_new.general_param;% read general_param from structure

%% for each embryo

for k = 1:length(saveVarsMat_new.params) % all embryo
    
    param = saveVarsMat_new.params{k};
    
    if param.status >= 0
                
        disp(param.sp1)
        
        % give the possibility to user to generate a ROI selection using polygonal selection tool: the mask will be saved
        path_tiff_stack = strcat(param.basepath , '/' , param.sp1 , '/');
        name_tiff_stack = param.stem_name;
        name = fullfile(path_tiff_stack,name_tiff_stack);
        
        [end_reading,siz,fitsnom] = read_init(name,param.format_image,param.sp7,param.sp2,param.sp3,'');
        [Image,~,error_reading] = read_with_trial(param.sp2,param.sp7,param.format_image,siz,fitsnom,'none',param.sp3,param.cortex_pass2.channel_interest_AC);
        
        % user will do a selection using roipoly tool
        mask_BW_AC = roipoly(imadjust_perso(Image));
        name_mask = strcat(param.stem_name,'_mask.tif');
        name_path_mask = fullfile(path_tiff_stack,name_mask);
        imwrite(mask_BW_AC,name_path_mask);
        param.cortex_pass1.mask_image = name_path_mask;
        param.cortex_pass2.mask_image = name_path_mask;
        clear Image
        
        saveVarsMat_new.params{k} = param ;
    end
    
end

save(xml_matfile, 'saveVarsMat_new');
close all
clear all

end

