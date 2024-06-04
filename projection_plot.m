figure

% Define the mean and covariance matrix for the normal distribution
mu = [0 0 0]; % mean
sigma = [0.2 0 0; 0 2 0; 0 0 4]; % covariance matrix

% Generate random samples from the normal distribution
nSamples = 10000; % number of samples
data = mvnrnd(mu, sigma, nSamples); % generates random samples from a multivariate normal distribution

x = data(:,1);
y = data(:,2);
z = data(:,3);
edges = -5:0.1:5;


fig = figure(1);
ax = gca(fig);
scatter3(x,y,z,...
    'filled','MarkerFaceAlpha',0.1,"MarkerEdgeAlpha",0.1,"MarkerEdgeColor","none")
grid off
hold(ax, 'on')

tform_x = hgtransform(ax);
tform_y = hgtransform(ax);
tform_z = hgtransform(ax);


yz = hist3([y,z], 'Edges',{edges,edges});
xy = hist3([x,y], 'Edges',{edges,edges});
zx = hist3([z,x], 'Edges',{edges,edges});

image(tform_x,edges,edges, yz)
tform_x.Matrix =  makehgtform('translate', ...
    [-6,0,0],...
    'yrotate',pi/2);

image(tform_y,edges,edges, zx)
tform_y.Matrix =  makehgtform('translate', ...
    [0,-6,0],...
    'xrotate',pi/2);

image(tform_z,edges,edges, xy)
tform_z.Matrix =  makehgtform('translate', ...
    [0,0,-6],...
    'zrotate',pi/2);

xlim([-6,6])
ylim([-6,6])
zlim([-6,6])
axis square

view([1 1 1])
xlabel('Px');
ylabel('Py');
zlabel('Pz');
