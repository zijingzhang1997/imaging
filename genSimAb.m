function [sParamCalib,sParamObj,sParamNoObj,freq,opts] = ...
    genSimAb(sParObjStruct,sParNoObjStruct,tagPortNum,rxPortNum,...
    tagPosition,rxPosition,freqIdx,roomSize,voxelSize,opts)
% Input -------------------------------------------------------------------
% sParamObj: S parameter structure with objects (people)
% sParamNoObj: S parameter structure with no object - for calibration
% tagPortNum: CST port numbers corresponding to the tags 
% rxPortNum: CST port numbers corresponding to the receivers
% posTag: Position of tags corresponding to tagPortNum
% posRecv: Position of receivers corresponding to rxPortNum
% freqIdx: Indices of frequencies to be extracted from sParamObj
% Output ------------------------------------------------------------------
% expConst == A in b = Ax, x is reflectivity at each voxel. Complex single.
% sParamCalib, b in above equation. Complex single float.
% -------------------------------------------------------------------------
%%

if ~isfield(opts,'saveData')
    opts.saveData = 0;
end
if ~isfield(opts,'isNoisy')
    opts.isNoisy = 0;
    if ~isfield(opts,'snr')
        opts.snr = 0;
    end
end
if ~isfield(opts,'phNoise')
    opts.phNoise = 0;
end
% if ~isfield(opts,'calibType')
%     opts.calibType = 1; % Different calibration types
%     fprintf('Using default calibration type 1.\n');
% end


%%
c = 3e8;                                 % Speed of light
nFreq = length(freqIdx);                 % Number of frequencies that we use for reconstruction using IFT
[nTag, ~] = size(tagPosition);                % Number of tags
[nRecv, ~] = size(rxPosition);              % Number of receivers
freq = sParObjStruct.Frequencies(freqIdx); % Extracting frequencies from S-parameter object in Hz
sParObjAllPort = sParObjStruct.Parameters; 
sParNoObjAllPort = sParNoObjStruct.Parameters;

% This is to the results with both S12 and S21 in simulation
% Only because Rx number is small - less simulation time.
% sParObj = permute(sParObj,[2 1 3]); 
% sParNoObj = permute(sParNoObj,[2 1 3]);
% % 
% figure;
% scatter3(tagPosition(:,1),tagPosition(:,2),tagPosition(:,3),'ro');
% hold on
% scatter3(rxPosition(:,1),rxPosition(:,2),rxPosition(:,3),'k*');
% axis 'equal'; axis 'tight';
% legend('Tags','Receivers');

%% Extract signal at each receiver corresponding to each tag and frequency
% and perform calibration

% sParameter arranged as (tag,Rx,freq), truncated from entire object
sParamObj = zeros(nTag,nRecv,nFreq); 
sParamNoObj = zeros(nTag,nRecv,nFreq); 

if opts.phNoise == 1
    if ~isfield(opts,'phErrStdDev')
        opts.phErrStdDev = 15;
        fprintf('Introducind phase error with standard dev of %2.1f degrees\n',opts.phErrStdDev);
    end
    %% Perturbing phase 
    fprintf('Introducing phase perturbation.\n');
    %phaseNoise = normrnd(0,opts.phErrStdDev,size(sParObjAllPort))*pi/180;
    phaseNoise = opts.phErrStdDev.*sign(randn(size(sParObjAllPort))).*rand(size(sParObjAllPort))*pi/180;
    sParamObjR = abs(sParObjAllPort);
    sParamObjTheta = angle(sParObjAllPort) + phaseNoise;
    sParObjAllPort = sParamObjR.*(cos(sParamObjTheta) + 1j * sin(sParamObjTheta));

    %phaseNoise = normrnd(0,opts.phErrStdDev,size(sParNoObjAllPort))*pi/180;
    %phaseNoise = opts.phErrStdDev.*sign(randn(size(sParObjAllPort))).*rand(size(sParObjAllPort))*pi/180;
    sParamNoObjR = abs(sParNoObjAllPort);
    sParamNoObjTheta = angle(sParNoObjAllPort) + phaseNoise;
    sParNoObjAllPort = sParamNoObjR.*(cos(sParamNoObjTheta) + 1j * sin(sParamNoObjTheta));
end

for freqNum = 1:nFreq
    for recvNum = 1:nRecv
        for tagNum = 1:nTag
            % sParamObj has S(2,1,:) == S21. Output at 2, input at 1. 
            % Hence, the extraction is sParamObj(Recv,Tx,frequencies).
            sParamObj(tagNum,recvNum,freqNum) = ...
                sParObjAllPort(rxPortNum(recvNum),tagPortNum(tagNum),freqIdx(freqNum));
            sParamNoObj(tagNum,recvNum,freqNum) = ...
                sParNoObjAllPort(rxPortNum(recvNum),tagPortNum(tagNum),freqIdx(freqNum));
        end
    end
end

% Need to clear unused variables, as this is huge data
clearvars sParObj sParNoObj sParObjStruct sParNoObjStruct sParNoObjAllPort sParObjAllPort

if opts.isNoisy == 1
    [sParamObj,sParamNoObj,opts.calcSNRObj,opts.calcSNRNoObj] = addNoise(sParamObj,sParamNoObj,opts.snr);
end
%%
% Considering S-parameters as not just the gain, but the signal received at
% receiver from a particular tag at a particular frequuency. Performing
% calibration, that includes background signal removal and initial phase
% cancellation, based on Yunfei's paper.
sParamCalib = zeros(nTag,nRecv,nFreq);

switch(opts.calibType)
    case 1
        % This is from Yunfei's paper and updated for 0 tag reading case
        fprintf('Using calibration type 1\n');
        for freqNum = 1:length(freq)
            lamda = c/freq(freqNum);
            for recvNum = 1:nRecv
                for tagNum = 1:nTag
                    distTagRx = norm(tagPosition(tagNum,:)-rxPosition(recvNum,:));
                    for recvPairNum = 1:nRecv
                        if recvPairNum ~= recvNum
                            sParamObjRatio = sParamObj(tagNum,recvNum,freqNum)...
                                /sParamObj(tagNum,recvPairNum,freqNum);
                            sParamNoObjRatio = sParamNoObj(tagNum,recvNum,freqNum)...
                                /sParamNoObj(tagNum,recvPairNum,freqNum);
                            sParamCalib(tagNum,recvNum,freqNum) = ...
                                sParamCalib(tagNum,recvNum,freqNum) + ...
                                (sParamObjRatio/sParamNoObjRatio-1)*...
                                exp(-1j*(2*pi/lamda)*distTagRx);
                        end
                    end
                end
            end
        end
        
            % ---------------------------------------------------------------------    
    case 2
        % This is based on Prof. Hysell's method for 0 tag reading case
        % Calibrate receivers, for each freq, rx
        scale = 1e-3;
        sParamObj  = scale.*sParamObj;
        sParamNoObj = scale.*sParamNoObj;
        fprintf('Using calibration type 2\n');
        rcal = zeros(length(freq),nRecv);
        for freqNum = 1:length(freq)
            dummy1 = zeros(nTag,nRecv);

            wn = 2*pi*freq(freqNum)/c;
            rcal(freqNum,1) = 1; % Using Rx1 as standard reference
            for recvNum = 1:nRecv
                sumComp = 0;
                for tagNum = 1:nTag % Average over tags
                    distTagRx1 = norm(tagPosition(tagNum,:) - rxPosition(1,:));
                    distTagRx = norm(tagPosition(tagNum,:) - rxPosition(recvNum,:));
                    % Expected relative phase
                    pse = -wn*(distTagRx - distTagRx1);
                    fac = cos(pse) + 1j*sin(pse);
                    % Measured relative phase - they are weighted by the
                    % power for that particular tag?
                    poff = sParamNoObj(tagNum,recvNum,freqNum)*sParamNoObj(tagNum,1,freqNum)';
                    
                    dummy1(tagNum,recvNum) = fac*poff';
                    
                    sumComp = sumComp+fac*poff';
                end
                % For final calibration, normalize it
                rcal(freqNum,recvNum) = sumComp/abs(sumComp);
            end
        end
        % Calibrate tags making use of rx calibration; for each freq, tag
        tcal = zeros(length(freq),nTag);
        for freqNum = 1:length(freq)
            wn = 2*pi*freq(freqNum)/c;
            for tagNum = 1:nTag
                sumComp = 0;
                for recvNum = 1:nRecv
                    distTagRx = norm(tagPosition(tagNum,:) - rxPosition(recvNum,:));
                    % Anticipated accumulated phase
                    pse = -wn*distTagRx;
                    fac = cos(pse) + 1j*sin(pse);
                    % Measured accumulated phase with rx correction
                    poff = sParamNoObj(tagNum,recvNum,freqNum)*rcal(freqNum,recvNum);
                    sumComp = sumComp + fac*poff';% cancel starting phase  of each tag
                end
                if abs(sumComp)>0
                    tcal(freqNum,tagNum) = sumComp/abs(sumComp);
                else
                    tcal(freqNum,tagNum) = 1;
                end
            end
        end
        for tagNum = 1:nTag
            for freqNum = 1:nFreq
                for recvNum = 1:nRecv
                    % This simply means removing LOS information. What if
                    % there is no LOS?
                    sParamCalib(tagNum,recvNum,freqNum) = sParamObj(tagNum,recvNum,freqNum) - sParamNoObj(tagNum,recvNum,freqNum);
                    fac = rcal(freqNum,recvNum) * tcal(freqNum,tagNum);
                    sParamCalib(tagNum,recvNum,freqNum) = sParamCalib(tagNum,recvNum,freqNum)*fac;
                end
            end
        end
        % --------------------------------------------------------------
    otherwise
        fprintf('Wrong calibration option selection, select a valid method. \n');
end




end