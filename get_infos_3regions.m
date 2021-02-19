function [ area_3regions,xStart_3regions,xEnd_3regions ] = get_infos_3regions( regionXlimit,regionXlength,regionArea,limit1,limit2 )

% to get limits and area of each 3 region

    area = regionArea.entireEmbryo.nbR10;
    xStart = regionXlimit.entireEmbryo.nbR1;
    deltaX = regionXlength.entireEmbryo.nbR10;
    
    for iTime = 1 : size(area,1)
        xStart_region1(iTime) = xStart(iTime);
        xStart_region2(iTime) = xStart(iTime) +limit1/10*deltaX(iTime);
        xStart_region3(iTime) = xStart(iTime) + limit2/10*deltaX(iTime);
        
        xEnd_region1(iTime) = xStart(iTime)+limit1/10*deltaX(iTime);
        xEnd_region2(iTime) = xStart(iTime) + limit2/10*deltaX(iTime);
        xEnd_region3(iTime) = xStart(iTime) + 10*deltaX(iTime);
        
        area_region1(iTime)=0;
        for i = 1 : floor(limit1/10)
            area_region1(iTime) = area_region1(iTime) + area(iTime,i);
        end
        area_region1(iTime) = area_region1(iTime) + mod(limit1/10,1)*area(iTime,floor(limit1/10)+1);
        
        area_region2(iTime)=0;
        for i = floor(limit1/10)+2 : floor(limit2/10)
            area_region2(iTime) = area_region2(iTime) + area(iTime,i);
        end
        area_region2(iTime) = area_region2(iTime) + (1-mod(limit1/10,1))*area(iTime,floor(limit1/10)+1);
        area_region2(iTime) = area_region2(iTime) + mod(limit2/10,1)*area(iTime,floor(limit2/10)+1);
        
        area_region3(iTime)=0;
        for i = floor(limit2/10)+2 : 10
            area_region3(iTime) = area_region3(iTime) + area(iTime,i);
        end
        area_region3(iTime) = area_region3(iTime) + (1-mod(limit2/10,1))*area(iTime,floor(limit2/10)+1);

    end

    area_3regions.region1 = area_region1; % area in squared pixels
    area_3regions.region2 = area_region2;
    area_3regions.region3 = area_region3;
    
    xStart_3regions.region1 = xStart_region1;
    xStart_3regions.region2 = xStart_region2;
    xStart_3regions.region3 = xStart_region3;
        
    xEnd_3regions.region1 = xEnd_region1;
    xEnd_3regions.region2 = xEnd_region2;
    xEnd_3regions.region3 = xEnd_region3;    
    
    
end

