function [imgBrightness] = visImg(imgComplex,roomSize,voxelSize)

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
imgthresh=0.6;
imgBrightness(imgBrightness<imgthresh) = 0;

imgBrightnessf = csaps({xVoxel,yVoxel,zVoxel},imgBrightness,0.9995);
imgSpline = fnval(imgBrightnessf,{xVoxel,yVoxel,zVoxel}) ;
figure
h = slice(x,y,z,imgSpline,xslice,yslice,zslice);
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
a = alphamap('rampup',256);
imgThresh = 150;
a(1:imgThresh)=0;
alphamap(a); 






figure
%h = slice(x,y,z,imgBrightness,[],[],[0.1]);
h = slice(x,y,z,imgBrightness,xslice,yslice,zslice);
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
a = alphamap('rampup',256);
imgThresh = 100;
a(1:imgThresh)=0;
alphamap(a); 
%title('11a6 R20    Zmin-50 Zmax 150   Xcenter600 Ycenter600 ')

end