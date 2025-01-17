function [tform_x,tform_y,tform_z] = projection_maps_scatter(p_dat,indx);


lim = 2.5;
ax = gca;
tform_x = hgtransform(ax); 
tform_y = hgtransform(ax); 
tform_z = hgtransform(ax);


tform_x.Matrix =  makehgtform('translate', ...
    [-lim,0,0],...
    'xrotate',pi/2,'yrotate',pi/2);

tform_y.Matrix =  makehgtform('translate', ...
    [0,-lim,0],...
    'zrotate',-pi/2,'yrotate',-pi/2);

tform_z.Matrix =  makehgtform('translate', ...
    [0,0,-lim]);

cmap = jet;
c_list = cmap(round(linspace(1,length(cmap),7)),:);

for k = indx
    x = p_dat{k}(:,1);
    y = p_dat{k}(:,2);
    z = p_dat{k}(:,3);
    plot(tform_x,y,z,'.',"Color",[c_list(k,:) 1]); hold on;
    plot(tform_y,z,x,'.',"Color",[c_list(k,:) 1]); hold on;
    plot(tform_z,x,y,'.',"Color",[c_list(k,:) 1]); hold on;
    text(tform_z,mean(x),mean(y),num2str(k))
    text(tform_x,mean(y),mean(z),num2str(k))
    text(tform_y,mean(z),mean(x),num2str(k))

end
xlim([-lim,lim]); ylim([-lim,lim]); zlim([-lim,lim]);
axis square

view([1 1 1])
xlabel('Px'); ylabel('Py'); zlabel('Pz');

end