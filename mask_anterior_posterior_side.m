function [maskedStack_rotated_anterior, maskedStack_rotated_posterior,maskedStack_rotated_anterior_reduced] = ...
    mask_anterior_posterior_side(maskedStack_rotated,furrow_position,regionXlimit,regionXlength)


% to mask the posterior or anterior blastomere

XSize = size(maskedStack_rotated,1);
YSize = size(maskedStack_rotated,2);

% below commented line that are necessary when trying to redo tracking on
% old movies because of delta time frame in sp2
% end_ = length(regionXlimit.entireEmbryo.nbR1);
% for i = 1 : 25
%     regionXlimit.entireEmbryo.nbR1(end_+i) = regionXlimit.entireEmbryo.nbR1(end_);
% end
% for i = 1 : 25
%     regionXlength.entireEmbryo.nbR1(end_+i) = regionXlength.entireEmbryo.nbR1(end_);
% end

%for k = 1 : size(maskedStack_rotated,3)
for k = 1 : length(regionXlimit.entireEmbryo.nbR1) % c briggsae
    furrow_position.xCoordinate.timeDependence_inverse(k,1) = regionXlimit.entireEmbryo.nbR1(k) + ...
        ( 100 - furrow_position.percent_length.timeDependence(k) )./100 .* regionXlength.entireEmbryo.nbR1(k) ;
end

maskedStack_rotated_anterior = zeros(XSize,YSize,size(maskedStack_rotated,3));
for i = 1 : size(maskedStack_rotated,3)
    for j = 1 : round( furrow_position.xCoordinate.timeDependence(i) )
        maskedStack_rotated_anterior (:,j,i) = ones(XSize,1); % posterior blastomere is masked
    end
end

maskedStack_rotated_anterior_reduced = zeros(XSize,YSize,size(maskedStack_rotated,3));
% for i = 1 : size(maskedStack_rotated,3)
%     for j = 1 : round( furrow_position.xCoordinate.timeDependence_inverse(i) )
%         maskedStack_rotated_anterior_reduced(:,j,i) = ones(XSize,1);
%     end
% end

maskedStack_rotated_posterior = NaN(XSize,YSize,size(maskedStack_rotated,3));
maskedStack_rotated_posterior(maskedStack_rotated_anterior == 1) = 0;
maskedStack_rotated_posterior(maskedStack_rotated_anterior == 0) = 1;

maskedStack_rotated_anterior = logical(maskedStack_rotated_anterior);
%maskedStack_rotated_anterior_reduced = logical(maskedStack_rotated_anterior_reduced);
maskedStack_rotated_posterior = logical(maskedStack_rotated_posterior);


end

