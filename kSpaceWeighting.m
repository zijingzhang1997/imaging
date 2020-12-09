% K-space weightings: Adding weights to K-space points
% Pragya Sharma, ps847@cornell.edu
% 10 June 2019

function [dataWeighted,kWeight,K,infoTagRxFreq] = kSpaceWeighting(data, ...
                               tagPosition,rxPosition, freq, posScat, opts)
% -------------------------------------------------------------------------
% Input: 
% data: Data input vector, in the same order as infoTagRxFreq.
% tagPosition: Position of tags in a column format with [x,y,z] coordinates
% rxPosition: Position of Rx in a column format with [x,y,z] coordinates
% freq: Frequency in column vector
% posScat: Position of scatterers. Their centroid is taken as origin.
% opts: Options for weighting
% -------------------------------------------------------------------------
% Output:
% dataWeighted: Data weighted according to the rules as defined in options.
% kWeight: Weights for each tag, rx, freq pair. Information for tag, rx,
% freq is given in infoTagRxFreq. In column format.
% K: [Kx,Ky,Kz] coordinates of K space in column format
% infoTagRxFreq: This is same length as K, with each column corresponding
% to respective [tag, rx, freq] row number from the input tagPosition,
% rxPosition and freq vectors.
% -------------------------------------------------------------------------
c = physconst('LightSpeed');

kMax = 2*max(freq)/c;

if ~isfield(opts,'weightType') % Wighting type of K-space
    opts.weightType = 'lp'; % Low pass: less weights for higher freq
end
if ~isfield(opts,'threshRel') % Knee point for lp/ hp relative to kMax
    opts.threshRel = 0.6;
end
if ~isfield(opts,'wtFactor') % How much weight factor after threshold. 
    opts.wtFactor = 0.6; % 0 = no weight. Maximum weightage is 1. 
end
if ~isfield(opts,'viewFigK')
    opts.viewFigK = 0; % By default do not view K-space figures.
end
if ~isfield(opts,'viewFig')
    opts.viewFig = 1; % By default, viewing figures of this function.
end
fprintf('K-space weighting: %s\n',opts.weightType);
fprintf('K-space weighting relative threshold: %3.2f\n',opts.threshRel);
fprintf('K-space attenuation: %3.2f\n',opts.wtFactor);

% Finding [Kx, Ky, Kz]

[K, infoTagRxFreq] = kspaceXYZ(tagPosition, rxPosition, freq, posScat, opts); 

kThresh = opts.threshRel*kMax;
kDist = vecnorm(K,2,2);  % 2????? 2????????
kWeight = ones(size(K,1),1); % Equal weights to all.

switch (opts.weightType)
    case 'lp'
        attnIdx = kDist >= kThresh; % Finding indices more than threshold
    case 'hp'
        attnIdx = kDist <= kThresh;
    otherwise
        warning('kSpaceWeighting: Enter correct opts.weightType.');        
end

% Currently, doing a strict filter, without any attenuation slope. Either
% weight is 1 or weight factor.
kWeight(attnIdx) = opts.wtFactor;
dataWeighted = kWeight .* data;
   
phi = 0:0.01:2*pi;
x = 2*max(freq)/c*cos(phi);
y = 2*max(freq)/c*sin(phi);

xLim = kThresh * cos(phi);
yLim = kThresh * sin(phi);

nonAttnIdx = not(attnIdx);

if opts.viewFig
    figure('Position',[400,300,300,300]);
    %ax(1) = subplot(3,1,1);
    plot(K(nonAttnIdx,1),K(nonAttnIdx,2),'bs'); hold on
    plot(K(attnIdx,1),K(attnIdx,2),'r+');
    xlabel('Kx'); ylabel('Ky');hold on;
    plot(x, y, 'k'); plot(xLim, yLim, 'r');
    axis square; axis equal; 

    figure('Position',[800,300,300,300]);
    plot(K(nonAttnIdx,1),K(nonAttnIdx,3),'bs'); hold on
    plot(K(attnIdx,1),K(attnIdx,3),'r+');
    xlabel('Kx'); ylabel('Kz');hold on;
    plot(x, y, 'k'); plot(xLim, yLim, 'r');
    axis square; axis equal;

    figure('Position',[1200,300,300,300]);
    plot(K(nonAttnIdx,2),K(nonAttnIdx,3),'bs'); hold on
    plot(K(attnIdx,2),K(attnIdx,3),'r+');
    xlabel('Ky'); ylabel('Kz');hold on;
    plot(x, y, 'k');plot(xLim, yLim, 'r');
    axis square; axis equal;

    figure
    plot3(K(nonAttnIdx,1),K(nonAttnIdx,2),K(nonAttnIdx,3),'bs');hold on
    plot3(K(attnIdx,1),K(attnIdx,2),K(attnIdx,3),'r+');
    plot(x,y,'k');
    plot3(x,zeros(length(x),1),y,'k');
    plot3(zeros(length(x),1),x,y,'k');
    plot(xLim,yLim,'r');
    plot3(xLim,zeros(length(xLim),1),yLim,'r');
    plot3(zeros(length(xLim),1),xLim,yLim,'r');
    xlabel('Kx'); ylabel('Ky'); zlabel('Kz');
end

end
