function curve = zero_contact(incurve,contact)
    index = get_index(incurve.ext(:,1), contact);
%     ext = incurve.ext(index:end, :);
%     ext(:,1) = ext(:,1) - min(ext(:,1));
%     ext(:,2) = ext(:,2) - min(ext(:,2));
%     curve = incurve;
%     curve.ext = ext;
    curve = incurve;
    curve.ext(:,1) = curve.ext(:,1) - curve.ext(index,1);
    curve.ext(:,2) = curve.ext(:,2) - curve.ext(index,2);
%     x_diff = curve.ext(end,1) - curve.ret(end,1);
%     y_diff = curve.ext(end,2) - curve.ret(end,2);
%     curve.ret(:,1) = curve.ret(:,1) + x_diff;
%     curve.ret(:,2) = curve.ret(:,2) + y_diff;
end

