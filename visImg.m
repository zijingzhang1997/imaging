function [imgBrightness,g] = visImg(imgComplex,roomSize,voxelSize)

xVoxel = roomSize(1,1):voxelSize(1): roomSize(1,2); 
yVoxel = roomSize(2,1):voxelSize(2): roomSize(2,2); 
zVoxel = roomSize(3,1):voxelSize(3): roomSize(3,2); 
% Combination of all of these to get coordinates for all the voxels
xyzVoxelCoord = combvec(xVoxel,yVoxel,zVoxel)';

nx = length(xVoxel);
ny = length(yVoxel);
nz = length(zVoxel);

imgBrightness = (abs(imgComplex).^2);
maxBrightness = max(imgBrightness(:)); 
minBrightness = min(imgBrightness(:));
imgBrightness = (imgBrightness-minBrightness)/(maxBrightness-minBrightness);

% Visualizing reconstructed image
x = reshape(xyzVoxelCoord(:,1),nx,ny,nz);
y = reshape(xyzVoxelCoord(:,2),nx,ny,nz);
z = reshape(xyzVoxelCoord(:,3),nx,ny,nz);

xslice = roomSize(1,1):voxelSize(1):roomSize(1,2);
yslice = roomSize(2,1):voxelSize(2):roomSize(2,2);
zslice = roomSize(3,1):voxelSize(3):roomSize(3,2);

p = [2 1 3];
x = permute(x, p);
y = permute(y, p);
z = permute(z, p);
imgBrightness = reshape(imgBrightness,[nx ny nz]);
imgBrightness = permute(imgBrightness, p);

%gradient detection instead of threshold


[gx,gy,gz] = gradient(imgBrightness);
g=sqrt(gx.^2 +gy.^2 +gz.^2);
% gthresh=0.2;
% g(g<gthresh) = 0;
% g=g./imgBrightness;
gthresh=0.1;
g(g<gthresh) = 0;
maxg=max(g(:));

% [~, clusters,a] = i4block_components(g, roomSize, voxelSize);
% 
% fprintf('Initial cluster number = %d\n',size(clusters.centroid,1));
% for i = 1:size(clusters.centroid,1)
%     fprintf('Clusters centroid: [%3.2f, %3.2f, %3.2f] with element number %d\n',clusters.centroid(i,:),clusters.elemNum(i));
% end
% fprintf('\n');
% opts.distTh = 0.4; % distance threshold, clusters with centers closer than this will be combined
% opts.XYdistTh = 0.3;
% opts.elemNumTh = 0.8; % clusters with element number less than 60% of the maximum will be rejected
% opts.minHeightRatio = 0.6; % Minimum height ratio compared to largest object, exact ht depends on voxel size etc.
% clusterOut = clusterProcess(clusters,opts);
% 
% centroid = clusterOut.centroid;
% elemNum = clusterOut.elemNum;
% for i = 1:size(centroid,1)
%     fprintf('Clusters centroid: [%3.2f, %3.2f, %3.2f] with element number %d\n',centroid(i,:),elemNum(i));
% end




figure
h = slice(x,y,z,g,xslice,yslice,zslice);
%h = slice(x,y,z,imgBrightness,xslice,yslice,zslice);
xlabel('x (m)','FontSize',14)
ylabel('y (m)','FontSize',14)
zlabel('z (m)','FontSize',14)
xlim([roomSize(1,1),  roomSize(1,2)])
ylim([roomSize(2,1),  roomSize(2,2)])
zlim([roomSize(3,1), roomSize(3,2)])

set(h, 'EdgeColor','none',...
    'FaceColor','interp',...
    'FaceAlpha','interp');
alpha('color')
end