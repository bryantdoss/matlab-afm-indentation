function z = flatten(curves, contacts, E, E_cutoff)

x = nan(1,length(curves));
y = nan(1,length(curves));
z = nan(1,length(curves));

for c=1:length(curves)
    if E(c) > E_cutoff
        x(c) = curves(c).tip_x;
        y(c) = curves(c).tip_y;
        z(c) = contacts(c);
    end
end
    
x = excise(x);
y = excise(y);
z = excise(z);

A = [sum(x.*x) sum(x.*y) sum(x);
sum(x.*y) sum(y.*y) sum(y);
sum(x) sum(y) length(z)];

B = [sum(x.*z) sum(y.*z) sum(z)];

R = A\B';

for c=1:length(curves)
    x(c) = curves(c).tip_x;
    y(c) = curves(c).tip_y;
    z(c) = contacts(c);
end

% close all

% figure(1)
% scatter(x,y,480,z,'filled')
% colorbar

z = z-x*R(1)-y*R(2)-R(3);

% figure(2)
% scatter(x,y,480,z,'filled')
% colorbar

% figure(3)
% xlin = linspace(min(x),max(x),14);
% ylin = linspace(min(y),max(y),14);
% [X,Y] = meshgrid(xlin,ylin);
% f = scatteredInterpolant(x',y',z');
% Z = f(X,Y);
% mesh(X,Y,-Z) %interpolated
% axis tight; hold on
% plot3(x,y,-z,'.','MarkerSize',15) %nonuniform

end