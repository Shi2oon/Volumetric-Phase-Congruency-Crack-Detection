function [factor]=factorUnits(input_unit,pixel_size)
%function return pixeled data to unit'ed to meter
switch input_unit
        case 'um' 
        % convert um --> m
            factor=pixel_size*10^-6;
        case 'mm'
        % convert mm --> m
            factor=pixel_size*10^-3;  
        case 'm'% convert pixels --> m
            factor =pixel_size; 
        otherwise
            error('Invalid input_unit')
end