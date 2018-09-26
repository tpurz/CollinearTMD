function [timeDelay,data] = PhaseInterp(data)

%% READ ME %%

%Takes the new data format (9 lined: 0-4 as usual, 5&6 give phase for t, 7&8 for tau)
% Returns Lines 1-4 (X&Y for Linear and FWM) and new t spacing (in ps)
%Example use: [pos.t,dataInterp(1:4,:)] = PhaseInterp(data);

%{
%Get data
parameters.Folder = 'C:\Users\Torben Purz\Documents\Universitaet\USA\Fall18\Physics515\Hanna_Matlab\'; 
parameters.parameters.Run = '010'; parameters.Run = '010';
data = csvread(strcat(parameters.Folder,'XCorrPP',parameters.Run,'.csv'),1,0);

s3 = false;
shift = false;
refshift = 1586.3;
underS = 0;
hcn = 1.00029*1239842; %hc in air 
cn = 1.00029*299792458;
sampling = 916; %Hz from DAQ, should save this value in LV
plotextra = true;
stgvel = data(1,7); % stage velocity in mm/s

[rows, columns] = size(data); paramRows = 1:5:rows; 
dataParams= data(paramRows,1:10); 
stgvel = data(1,7); % stage velocity in mm/s
numt =  length(data(1,:));
pos.t = linspace(0,numt*stgvel*0.002*10^12/(cn*sampling),numt); % t steps starting from zero, in units of ps?
numtau = length(data(:,1))/9; 
pos.tauSteps = unique(dataParams(:,5)); 
pos.tau = 2*(pos.tauSteps(1:numtau)-pos.tauSteps(1))/0.1499; % what is this factor? should be 0.3?
pos.T = data(1,4);
procDataT= zeros(numt,6,numtau);
procData = zeros(numt,6, numtau);
counter = 0;
for n = 1:8:(rows-rows/9);
    counter = counter + 1;
    dataRows = (n+counter):(n+counter+5);
    targetRows = (n:(n+3));
    procDataT(:,:,counter) = data(dataRows,:)';
end

%}



%%


%%

%Get Phase
%phase=data(:,7,:)+data(:,6,:)/(2^30); 
phase=squeeze(data(:,6,:)+data(:,5,:)/(2^30)); 




%Translate into stage position using wavelength from phase-stage pos.
phase_tot = phase(size(phase,1),:)-phase(1,:); %Array with value for each tau
distance_tot = 2*size(phase,1)/1000*data(1,7,1)*1e-3; %Factor 2 because of path passed twice, data(1,7,:) always the same
lambda = distance_tot./(phase_tot);
%Average over all wavelengths to get a good estimate (Doesn't work for
%chirped pulses?)
%lambdaMean=mean(lambda)
translation=phase.*lambda; %Is different for each tau position, therefore dotwise product


%Declare translation "new x value" and interpolate Lines 1-4 (2-5) with it
step=(translation(size(data,2),1)-translation(1,1))/(size(data,2)-1); %Determined by first tau position
translation_new=translation(1,1):step:translation(size(data,2),1); %Determined by first tau position




for j=1:size(data,3) %along tau 
    for i=1:4
        data(:,i,j) = interp1(translation(:,j),data(:,i,j),translation_new);
    end
end

%Convert distance to time delay
c = 299792458;
timeDelay = (translation_new/c)*1e12;

%Return interpolated data and new t.pos in ps
end