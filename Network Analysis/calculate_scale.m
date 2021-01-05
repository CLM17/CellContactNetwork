function scale = calculate_scale(magnification, fieldSize)
    %----------------------------------------------------------------------
    % This function calculates the scale in um/pixel.
    
    % INPUT PARAMETERS
    %-----------------
    
    % magnification: string
    % Magnification with which the image is made, either '10x' or '20x'.
    
    % fieldSize: double
    % Size of 1 field in CX7. Examples: 1104 or 2208.
    
    % OUTPUT PARAMETER
    %-----------------
    
    % scale: double
    % size of 1 pixel in um.
    %----------------------------------------------------------------------

    if strcmp(magnification, 'M10')
        scale = 873.96 / fieldSize; % 873.96 is field size for 10x in CX7
    elseif strcmp(magnification, 'M20')
        scale = 441.41 / fieldSize; % 441.41 is field size for 10x in CX7
    else
        error('Invalid choice for magnification')
    end
    
end