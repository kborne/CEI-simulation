function [tform_x,tform_y,tform_z] = projection_maps(x,y,z,...
                        edges,scatter_q);


lim = 2;
% fig = figure(1);
ax = gca;

if scatter_q

% 
rand_ind = randperm(length(x));
plot_ind = rand_ind(1:20000);
    scatter3(x(plot_ind),y(plot_ind),z(plot_ind),...
        2.5,...
        'filled','MarkerFaceAlpha',0.15,"MarkerEdgeAlpha",0.5,"MarkerEdgeColor","none")
end
% grid off
hold(ax, 'on')

tform_x = hgtransform(ax);
tform_y = hgtransform(ax);
tform_z = hgtransform(ax);


yz = hist3([y,z], 'Edges',{edges,edges});
xy = hist3([x,y], 'Edges',{edges,edges});
zx = hist3([z,x], 'Edges',{edges,edges});

image(tform_x,edges,edges, yz)
tform_x.Matrix =  makehgtform('translate', ...
    [-lim,0,0],...
    'yrotate',-pi/2);

image(tform_y,edges,edges, zx)
tform_y.Matrix =  makehgtform('translate', ...
    [0,-lim,0],...
    'xrotate',pi/2);

image(tform_z,edges,edges, xy)
tform_z.Matrix =  makehgtform('translate', ...
    [0,0,-lim],...
    'zrotate',-pi/2,'yrotate',pi);

xlim([-lim,lim]); ylim([-lim,lim]); zlim([-lim,lim])
axis square

% view(2)
view([1 1 1])
xlabel('Px'); ylabel('Py'); zlabel('Pz');


end