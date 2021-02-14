% rotation_minimum.m
% Automatic contact point determination

% Basically introduce a virtual deflection then pick the "minimum point
% forcescale/lengthscale is the virtual deflection magnitude, this can be
% optimized somehow (for future work?)

% Input parameters -- SI units [Pa,m,radians]:
%  incurve: the force-curve structure
%  degree: amount of rotation
% Output parameters -- SI units [m,Pa]:
%  contact: the contact point

function contact = rotation_minimum (incurve, degree)
    if nargin < 2
        degree = 2;
    end
    
    forcescale = (incurve.ext(end,2) - incurve.ext(1,2)) / degree;
    lengthscale = incurve.ext(end,1) - incurve.ext(1,1);
    
    % Basically fabricate a virtual deflection for the enitre curve
    scale = forcescale / lengthscale;
    rotcurve = incurve;
    rotcurve.ext(:,2) = rotcurve.ext(:,2) - rotcurve.ext(:,1) * scale;
    
    % Find the minimum in the rotated curve, which is our contact        
    % point
    [~, index] = min(rotcurve.ext(:,2));        
    contact = rotcurve.ext(index,1);
    
    % Optional: plot
%     plot(incurve.ext(:,1),incurve.ext(:,2), rotcurve.ext(:,1),rotcurve.ext(:,2));
%     hold on
%     plot(contact, incurve.ext(index,2), 'kx', 'LineWidth', 2);
%     
end
