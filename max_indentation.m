function [indentation] = max_indentation(curve, contact)
    
    i = get_index(curve.ext(:,1), contact);
    i = i(1);
    
    indentation = curve.ext(end,1) - curve.ext(i,1);
%     indentation = curve.ext(end,2) - min(curve.ext(:,2));
end
