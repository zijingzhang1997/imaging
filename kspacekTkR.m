% K-space 
% Pragya Sharma, ps847@cornell.edu
% 20 May 2019, Updated on 9 May 2019.

function [K,infoTagRxFreq] = kspacekTkR(tagPosition, rxPosition, freq, ...
    posScat, opts)
% kspaceView is to visualize K-space.
% **This depends on scatterer position.**
% -------------------------------------------------------------------------
% Input: 
% tagPosition: Position of tags in a column format with [x,y,z] coordinates
% rxPosition: Position of Rx in a column format with [x,y,z] coordinates
% freq: Frequency in column vector
% posScat: Position of scatterers. Their centroid is taken as origin.
% -------------------------------------------------------------------------
% Output:
% K: [Kx,Ky,Kz] coordinates of K space in column format
% infoTagRxFreq: This is same length as K, with each column corresponding
% to respective [tag, rx, freq] row number from the input tagPosition,
% rxPosition and freq vectors.
% -------------------------------------------------------------------------

% Options
if ~isfield(opts,'viewFigK')
    opts.viewFigK = 1;
end

c = physconst('LightSpeed');

nTag = size(tagPosition,1);
nRecv = size(rxPosition,1);
nFreq = size(freq,1);


xScat = mean(posScat(:,1));
yScat = mean(posScat(:,2));
zScat = mean(posScat(:,3)); 

tagNum = (1:nTag)';
rxNum = (1:nRecv)';
freqNum = (1:nFreq)';

kT = tagPosition - [xScat, yScat, zScat];
kT = kT./vecnorm(kT,2,2);
kR = rxPosition - [xScat, yScat, zScat];
kR = kR./vecnorm(kR,2,2);

kT = repmat(kT,nRecv,1); 
tagNum = repmat(tagNum(:),nRecv,1);
kR = repelem(kR,nTag,1);
rxNum = repelem(rxNum(:),nTag,1);

kTkR = repmat((kT+kR),nFreq,1);
tagNum = repmat(tagNum,nFreq,1);
rxNum = repmat(rxNum,nFreq,1);
fScale = freq./c;
fScale = repelem(fScale,nTag*nRecv,1);
freqNum = repelem(freqNum(:),nTag*nRecv,1);

infoTagRxFreq = [tagNum,rxNum,freqNum];

K = fScale.*kTkR;
clearvars kT kR fScale tagNum rxNum freqNum

phi = 0:0.01:2*pi;
x = 2*max(freq)/c*cos(phi);
y = 2*max(freq)/c*sin(phi);

if opts.viewFigK == 1
    figure('Position',[400,300,300,300]);
    %ax(1) = subplot(3,1,1);
    plot(K(:,1),K(:,2),'r+');
    xlabel('Kx'); ylabel('Ky');hold on;
    plot(x, y, 'k');
    axis square; axis equal; 

    figure('Position',[800,300,300,300]);
    plot(K(:,1),K(:,3),'go');
    xlabel('Kx'); ylabel('Kz');hold on;
    plot(x, y, 'k');
    axis square; axis equal;

    figure('Position',[1200,300,300,300]);
    plot(K(:,2),K(:,3),'bs');
    xlabel('Ky'); ylabel('Kz');hold on;
    plot(x, y, 'k');
    axis square; axis equal;

    figure
    plot3(K(:,1),K(:,2),K(:,3),'r+');hold on
    plot(x,y,'k');
    plot3(x,zeros(length(x),1),y,'k');
    plot3(zeros(length(x),1),x,y,'k');
    xlabel('Kx'); ylabel('Ky'); zlabel('Kz');
end

