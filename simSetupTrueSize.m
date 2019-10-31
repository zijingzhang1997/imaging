% RF Imaging with tags and receivers. 
% Pragya Sharma (ps847@cornell.edu)
% April 10, 2018

%% Loading S-parameters and Receiver position

% function [imgMF,imgFISTA,roomSize,voxelSize] = simSetup16tag(dataName)
dataPath = ['E:\RFID imaging\pragya code\']; % Enter path to data
dataName = '5b_openaddspace1030';
fileName1 = ['RFimaging_3D_Planar_TrueSized_',dataName,'.s44p']; % With object
fileName2 = 'RFimaging_3D_Planar_TrueSized_5_openaddspace_noObj.s44p'; % Without object for calibration
% rxXYZname = [dataName,'.mat'];
posRxTxNum = '5';
opts.isNoisy = 0;
opts.usePhi = 0; opts.genPhi = 0; opts.loadPhi = 0; 
opts.genLib = 0; 
opts.sc = 1; % Stopping k of dictionary based OMP
opts.kVal = 1; % Number of objects + 1
opts1.fileName = 'lib9';
opts.RRfig = 0;

% Specify start, stop frequencies and number of frequencies used for
% reconstruction. Frequency range is [.6, 1.1) GHz.
freq = struct;
freq.Start = 0.9e9; %0.86 %0.3175
freq.Stop = 0.93e9; %0.9 %0.3265
freq.Num = 6;

[sParObjStruct,freqIdx] = readTouchstone(fileName1,dataPath,freq);
[sParNoObjStruct,    ~] = readTouchstone(fileName2,dataPath,freq);

%% Tag and Rx port numbers and corresponding locations
% Tags and receivers are arranged anti-clockwise. Unit is in meters
[rxPortNum,tagPortNum,rxPosition,tagPosition] = posRxTxTrueSizeSim(posRxTxNum);   

% Perturbing antenna locations
 maxLocErr = 0*0.001; % 20 mm * 0.001 = 0.02 m error
 % Introducing random error in recording of tag and receiver position,
 % bounded by maxLocErr. As error can be positive or negative, choosing that
 % randomly as well.
tagPosition = tagPosition + maxLocErr.*sign(randn(size(tagPosition,1),3)).*rand(size(tagPosition,1),3); 
%rxPosition = rxPosition + maxLocErr.*sign(randn(size(rxPosition,1),3)).*rand(size(rxPosition,1),3); 


%% Inverse object reflectivity reconstruction for each voxel
% codePath1 = 'E:\ArpaE2018\3DImaging_Simulation\CST_SimulationDataAnalysis\Algorithms\LeastSquare';
% codePath2 = 'E:\ArpaE2018\3DImaging_Simulation\CST_SimulationDataAnalysis\Algorithms\Fourier';
% codePath3 = 'E:\ArpaE2018\3DImaging_Simulation\CST_SimulationDataAnalysis\Algorithms\MP';
% codePath4 = 'E:\ArpaE2018\3DImaging_Simulation\CST_SimulationDataAnalysis\Algorithms\Cluster';
% 
% addpath(codePath1);addpath(codePath2);addpath(codePath3);addpath(codePath4);
% savePathLS = ['E:\ArpaE2018\3DImaging_Simulation\CST_SimulationDataAnalysis\',...
%     'Algorithms\LeastSquare\ProcessedData'];
% Room size in meters, each row indicates x, y and z limits
roomSize = [0.2, 3.4; 0.2, 3.4; -0.6, 0.6]; 
% voxel size in meters, each row indicates x, y and z direction values
voxelSize = [0.04;0.04;0.04]; 
xVoxel = roomSize(1,1):voxelSize(1): roomSize(1,2); 
yVoxel = roomSize(2,1):voxelSize(2): roomSize(2,2); 
zVoxel = roomSize(3,1):voxelSize(3): roomSize(3,2); 
nVoxel = [length(xVoxel), length(yVoxel), length(zVoxel)];

% [imgFou,~,~,~,~] = ...
%     fourierImagSim(sParObjStruct,sParNoObjStruct,tagPortNum,rxPortNum,...
%     tagPosition,rxPosition,freqIdx,roomSize,voxelSize);

% [imgBrightness,xyzVoxelCoord,sParamObj,sParamNoObj,freq] = ...
%     leastSquare1(sParObjStruct,sParNoObjStruct,tagPortNum,rxPortNum,...
%     tagPosition,rxPosition,freqIdx,dataName);

% A = expConst, b = sParamCalib
opts.phNoise = 1;  opts.phErrStdDev = 0;
%%
rng(0)
opts.snr = 10;
opts.calibType = 2;
[A,b,sParamObj,sParamNoObj,freq,opts] = ...
    genSimAb(sParObjStruct,sParNoObjStruct,tagPortNum,rxPortNum,...
    tagPosition,rxPosition,freqIdx,roomSize,voxelSize,opts);
if opts.isNoisy
    fprintf('Calculated SNR for obj %f and no obj cases %f. \n',...
        opts.calcSNRObj,opts.calcSNRNoObj);
end

% clear sParamObj sParamNoObj sParObjStruct sParNoObjStruct 
% save([dataPath,'\data_',dataName,'_pos',num2str(posRxTxNum)],'freq','tagPosition','rxPosition','sParamObj','sParamNoObj');
% close all


% K space
%  posScat = [1.8,1.8,0.3]; % Considering middle of the capture volume
%  opts3 = [];
%  [K,infoTagRxFreq] = kspaceXYZ(tagPosition, rxPosition, freq, posScat,opts3);
%  
%  opts3.weightType = 'lp';
%  opts3.wtFactor = 0.4;
%  opts3.threshRel = 0.7;
%  [bWeighted,kWeight,K,infoTagRxFreq] = kSpaceWeighting(b, ...
%                                 tagPosition,rxPosition, freq, posScat, opts3);
% close all;
%% Algorithms for image reconstruction or location detection

[imgMF,imgMFcomp,g] = matchFilt(A,b,roomSize,voxelSize);view(0,90); 

% xticks(0.2:0.1:0.7)
% savefig([dataPath,'\data_16tx8rx_',dataName,'_mf_',num2str(opts.snr)]);
[xyzVoxelCoord,~,~,~] = genXYZ(roomSize,voxelSize);

[maxImg,maxIdx] = (max(abs(imgMFcomp)));
pkVoxel = xyzVoxelCoord(maxIdx,:);

thresh = maxImg*0.8;
imgMFabs = reshape(abs(imgMFcomp),nVoxel);
imgMFabs(imgMFabs<thresh) = 0;
[~, clusters,a] = i4block_components(g, roomSize, voxelSize);

fprintf('Initial cluster number = %d\n',size(clusters.centroid,1));
for i = 1:size(clusters.centroid,1)
    fprintf('Clusters centroid: [%3.2f, %3.2f, %3.2f] with element number %d\n',clusters.centroid(i,:),clusters.elemNum(i));
end
fprintf('\n');
opts.distTh = 0.4; % distance threshold, clusters with centers closer than this will be combined
opts.XYdistTh = 0.3;
opts.elemNumTh = 0.8; % clusters with element number less than 60% of the maximum will be rejected
opts.minHeightRatio = 0.6; % Minimum height ratio compared to largest object, exact ht depends on voxel size etc.
clusterOut = clusterProcess(clusters,opts);

centroid = clusterOut.centroid;
elemNum = clusterOut.elemNum;
for i = 1:size(centroid,1)
    fprintf('Clusters centroid: [%3.2f, %3.2f, %3.2f] with element number %d\n',centroid(i,:),elemNum(i));
end
% -------------------------------------------------------------------------
% imgCG = ...
%     lsCG(A,b,roomSize,voxelSize,savePathLS,dataName);
% savefig([dataPath,'\data_16tx8rx_',dataName,'_cg']);

% ------------------------------------------------------------------------- 
% imgFISTA = lsFISTA_l1(A,b,roomSize, voxelSize,savePathLS,dataName);view(0,0);
% xticks(0.2:0.1:0.7)
% savefig([dataPath,'\data_16tx8rx_',dataName,'_fista_',num2str(opts.snr)]);

% -------------------------------------------------------------------------
% savePathBP = ['E:\ArpaE2018\3DImaging_Simulation\CST_Simulation',...
%         'DataAnalysis\Algorithms\MP\LibATrueSize\'];
% 
% % Voxel based OMP
% opts.RRfig = 0;
% if opts.usePhi == 1
%     if opts.genPhi == 1
%         opts2.imSz = roomSize;
%         opts2.vxSz = voxelSize;
%         opts2.A = A;
%         opts2.reconstruction = 0;
%         opts2.wName = 'db20';
%         phi = wavFilt(opts2);
%         save(['phi_',opts2.wName,'_pos',num2str(posRxTxNum)],'phi');
%     elseif opts.loadPhi == 1
%            load('phi_db20_pos1.mat');
%     end
%     opts2.imgWav = voxelOMP(phi,b,opts,roomSize,voxelSize,savePathBP,dataName);view(0,90);
%     % Take inverse wavelet Tx
%     opts2.reconstruct = 1;
%     imgVoxOMPphiC = wavFilt(opts2);
%     imgVoxOMPphi = visImg(imgVoxOMPphiC,roomSize,voxelSize);
%     title('OMP after inverse wavelet transform')
%     a = alphamap('rampup',256);
%     imgThresh = 0;
%     a(1:imgThresh)=0;
%     alphamap(a); 
%     clear phi opts.imgWav opts2
% 
% else
%     imgVoxOMP = voxelOMP(A,b,opts,roomSize,voxelSize,savePathBP,dataName);view(0,90);
% end
% imgVoxOMPabs = reshape(abs(imgVoxOMP),nVoxel);
% imgVoxOMPabs(imgVoxOMPabs>0) = 1;
% [c, clusters] = i4block_components(imgVoxOMPabs, roomSize, voxelSize );
% 
% % -------------------------------------------------------------------------
% % Dictionary based OMP
% if opts.genLib == 1
%     % This is optional. If library is already created, load it.    
%     opts1.imSize = roomSize;
%     opts1.freq = freq(:);
%     opts1.posRxTx = posRxTxNum;
%     opts1.voxSz = [0.04;0.04;0.04]; % meters
%     opts1.objSz = [0.4,0.4,0.5];
%     opts1.seeImg = 0;
%     opts1.savepath = savePathBP;
%     [objCenterGrid] = genLibraryTrueSizeMP(opts1); 
% end
% opts.RRfig = 1;
% [objDictOMP,xDictOMP,normResDictOMP] = dictOMP(b,opts,savePathBP,opts1.fileName);
% % clear A;
% % rmpath(codePath1); rmpath(codePath2);
% 
