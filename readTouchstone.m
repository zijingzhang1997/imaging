% Reach touchstone file v1 and v2 formats containing S parameters.
% Version 1: By convention, it uses the extension of '.snp', where "n" is
% the number of network ports of the device. With 36 tags and 12 receivers,
% we have file format of '.s48p'. The data is saved in RI (real - imag)
% format.
% Pragya Sharma (ps847@cornell.edu)
% April 10, 2018

function [sobj,freqIdx] = readTouchstone(fileName1,dataPath,freq)
% Input: 
% fileName = 'RFimaging_36Tx_12Rx_3D_RxXYZ1.s48p';
% dataPath = ['D:\Research\ARPAE_RFImaging\ARPAE_SharedFiles\DataShare',...
% '\CodeData2\Data']; % Enter path to data

%%  
addpath(dataPath);

%% Get S-parameters
sobj = sparameters(fileName1); % Read the touchstone file
freqSampled = sobj.Frequencies; % Extract the frequency vector from data structure
freqReqd = linspace(freq.Start, freq.Stop, freq.Num)'; % Required frequencies

freqIdx = zeros(freq.Num,1); % Initializing frequency indices
for iter = 1:freq.Num
    % Finding first index greater than each required frequency. Range might
    % be sligtly more than the stop frequency.
    freqIdx(iter) = find(freqSampled == freqReqd(iter),1); 
end

%%
rmpath(dataPath);

end
%% 