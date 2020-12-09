% dataPath = ['C:\RFID imaging\smallsize sim\ML dataset\TOUCHSTONE files\'];
% 
% 
N=80*4*51;
% dataNum=121;
% class=6;
% feat=zeros(dataNum,4,N);
% label=zeros(dataNum,2);
% a=linspace(50,150,11);
% b=linspace(50,150,11);
% [A,B] = meshgrid(a,b);
% A=A(:);
% B=B(:);
% iter=0;
% for i=1:dataNum
%     dataname = sprintf('%04d',i);
%     feat(i,:,:) = feature(dataname,dataPath);
%     iter=iter+1;
% end
% %  square =a*b  linear categorize 
% % for j=1:dataNum   
% %     volume(j)=A(j)*B(j);
% %     label(j)=(volume(j)-A(1)*B(1))/(A(dataNum)*B(dataNum)-A(1)*B(1));
% %     label(j)=(label(j)*(class-1))+1;
% %     label(j)=round(label(j));
% % end
% for j=1:dataNum   
%     
%     label(j,1)=A(j);
%     label(j,2)=B(j);
%     
% end
% data.featVec=feat;
% data.labelVec=label;
% save('data_train.mat','data');

%% test data
dataPath_test = ['C:\RFID imaging\smallsize sim\ML dataset\test\TOUCHSTONE files\'];

dataNum_test=25;
class=6;
feat_test=zeros(dataNum_test,4,N);
label_test=zeros(dataNum_test,2);
a_test=linspace(65,145,5);
b_test=linspace(65,145,5);
[A_test,B_test] = meshgrid(a_test,b_test);
A_test=A_test(:);
B_test=B_test(:);
iter_test=0;
for i=1:dataNum_test
    dataname = sprintf('%04d',i);
    feat_test(i,:,:) = feature(dataname,dataPath_test);
    iter_test=iter_test+1;
end
%  square =a*b  linear categorize 
for j=1:dataNum_test   
    
    label_test(j,1)=A_test(j);
    label_test(j,2)=B_test(j);
    
end
clear data
data.featVec=feat_test;
data.labelVec=label_test;
save('data_test.mat','data');


function feat = feature(dataName,dataPath)
%dataPath = ['C:\RFID imaging\smallsize sim\ML dataset\TOUCHSTONE files\']; % Enter path to data

fileName1 = ['MWS-sweep-01-',dataName,'.s84p']; % With object
fileName2 = 'noobj.s84p'; % Without object for calibration
% rxXYZname = [dataName,'.mat'];
posRxTxNum = '11';
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
freq.Start = 0.86e9; %0.86 %0.3175
freq.Stop = 0.91e9; %0.9 %0.3265
freq.Num = 51;

[sParObjStruct,freqIdx] = readTouchstone(fileName1,dataPath,freq);
[sParNoObjStruct,    ~] = readTouchstone(fileName2,dataPath,freq);

%% Tag and Rx port numbers and corresponding locations
% Tags and receivers are arranged anti-clockwise. Unit is in meters
[rxPortNum,tagPortNum,rxPosition,tagPosition] = posRxTxTrueSizeSim(posRxTxNum);   

% Perturbing antenna locations
 maxLocErr = 0; % 20 mm * 0.001 = 0.02 m error
 % Introducing random error in recording of tag and receiver position,
 % bounded by maxLocErr. As error can be positive or negative, choosing that
 % randomly as well.
%tagPosition = tagPosition + maxLocErr.*sign(randn(size(tagPosition,1),3)).*rand(size(tagPosition,1),3); 
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
roomSize = [0.2, 1.0; 0.2, 1.0; -0.6, 0.6]; 
% voxel size in meters, each row indicates x, y and z direction values
voxelSize = [0.02;0.02;0.02]; 
pixel=voxelSize(1);
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
opts.phNoise = 0;  opts.phErrStdDev = 0;
%%
rng(0)
opts.snr = 10;
opts.calibType = 2;
[sParamCalib,sParamObj,sParamNoObj,freq,opts] = ...
    genSimAb(sParObjStruct,sParNoObjStruct,tagPortNum,rxPortNum,...
    tagPosition,rxPosition,freqIdx,roomSize,voxelSize,opts);

a=1e5;
sParamObj_ph=angle(sParamObj);
sParamObj_abs=abs(sParamObj);
sParamCalib_ph=angle(sParamCalib);
sParamCalib_abs=abs(sParamCalib);
sParamObj_abs = sParamObj_abs(:);
sParamObj_ph = (sParamObj_ph(:));
sParamCalib_abs = (sParamCalib_abs(:));
sParamCalib_ph = (sParamCalib_ph(:));
sParamObj_abs=sParamObj_abs*a;
sParamCalib_abs=sParamCalib_abs*a;

feat=[sParamObj_abs sParamObj_ph sParamCalib_abs sParamCalib_ph];
feat=permute(feat,[2 1]);

end
