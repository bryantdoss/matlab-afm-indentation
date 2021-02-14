function i = get_index_map (curves, xpos, ypos)
x = [curves.tip_x]-(xpos);
y = [curves.tip_y]-(ypos);
[~,i]=min(abs(x)+abs(y));
end