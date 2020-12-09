% K-space 
% Pragya Sharma, ps847@cornell.edu
% 20 May 2019, Updated on 9 May 2019.

function [K,infoTagRxFreq] = kspaceXYZ(tagPosition, rxPosition, freq, ...
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

thetaTag = zeros(nTag,1);  % Angle from XY plane between XY plane and tag.
phiTag = zeros(nTag,1);  % Angle from X-axis in X-Y plane.
for i=1:nTag
    x_d = tagPosition(i,1)-xScat;
    y_d = tagPosition(i,2)-yScat;
    z_d = tagPosition(i,3)-zScat;
%     thetaTag(i) = acos(sqrt(x_d^2 + y_d^2)/sqrt(x_d^2 + y_d^2 + z_d^2));
    thetaTag(i) = atan(z_d/sqrt(x_d^2 + y_d^2));
    phiTag(i) = angle(x_d + 1j*y_d);
end

thetaRecv = zeros(nRecv,1);
phiRecv=zeros(nRecv,1);
for i=1:nRecv
    x_d = rxPosition(i,1)-xScat;
    y_d = rxPosition(i,2)-yScat;
    z_d = rxPosition(i,3)-zScat;
%     thetaRecv(i) = acos(sqrt(x_d^2 + y_d^2)/sqrt(x_d^2 + y_d^2 + z_d^2));
    thetaRecv(i) = atan(z_d/sqrt(x_d^2 + y_d^2));
    phiRecv(i)= angle(x_d + 1j*y_d);
end

phi = 0:0.01:2*pi;
x = 2*max(freq)/c*cos(phi);
y = 2*max(freq)/c*sin(phi);

count=1;
K = zeros(nTag*nRecv*nFreq,3);
infoTagRxFreq = zeros(nTag*nRecv*nFreq,3);

for freqNum=1:length(freq)
   for rxNum=1:length(phiRecv)
       for tagNum=1:length(phiTag)
            K(count,1) = (freq(freqNum)/c)*(cos(thetaTag(tagNum))*cos(phiTag(tagNum)) + ...
                cos(thetaRecv(rxNum))*cos(phiRecv(rxNum)));
            K(count,2) = (freq(freqNum)/c)*(cos(thetaTag(tagNum))*sin(phiTag(tagNum)) + ...
                cos(thetaRecv(rxNum))*sin(phiRecv(rxNum)));
            K(count,3) = (freq(freqNum)/c)*(sin(thetaTag(tagNum)) + sin(thetaRecv(rxNum)));
            infoTagRxFreq(count,:) = [tagNum,rxNum,freqNum];
            K_rad = norm(K(count,:)); 
            count=count+1;
        end
    end
end

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

