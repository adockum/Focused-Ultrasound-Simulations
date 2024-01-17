% Slicer2Kwave - Free Field Multielement Steering
% Written by Allison Dockum 
% Last Updated - 01/17/2024

clear all; close all; clc; 
%% Load in File Names - files are centered at slicer origin

% labelled xdcr file - labelled file has the elements numbered
% will carry through any transforms applied to the Nifti in Slicer 
fnlabeledxdcr = 'IGT_labeledOrigin_hard.nii.gz'; 

% unlabelled xdcr file w/ focus
fnxdcr = 'IGTorigin_hard.nii.gz'; 

% element locations [mm] - used for phasing array elements
load('IGTelemLocs.mat'); 

% Output filename
fnout = 'pout_-6.3_-1_0_1.4ff'; 
%% Transducer Setting 

freq = 650e3;                      % [Hertz]

%steerLoc is in US-space
steerLoc = [0 0 0];                % [mm] [X Y Z]

%Amplitude - can set manually or calculate 
%Amp = 100000;

% calculate amplitude from getInputAmplitude
ffPressure = 1.4;
Amp = getInputAmplitude(ffPressure,'steer',steerLoc,'max');  % amplitude at element surface

%% Load in the data to set simulation size 

xdcr = niftiread(fnxdcr); 

dim = size(xdcr); 
'Your simpace size'
dim
%checkFactors(min(dim),max(dim)+100)

% pad sim space to prime factors
padsize = [432 432 432]; 

% Get geometric focus before padding 
[~,I] = max(xdcr(:));
[l, m, n]= ind2sub(dim,I);
focus_pos_init = [l, m, n];

%% Set up simulation and define medium 

% Calculate padding
padpre = floor((padsize-dim)/2);
padpost = ceil((padsize-dim)/2);

% Pad arrays
xdcr = padarray(xdcr,padpre,0,'pre');
xdcr = padarray(xdcr,padpost,0,'post');
dim = padsize;

medium.sound_speed = ones(dim)*1500;
medium.density = ones(dim)*1000;
medium.alpha_coeff = ones(dim)*0.2;
medium.alpha_power = 1.1; 

%% Setup simulation grid and input signal

vox = 1e-3.*[0.25 0.25 0.25];            % [mm->m]

kgrid = kWaveGrid(dim(1), vox(1), dim(2), vox(2), dim(3), vox(3));
[kgrid.t_array, ~] = makeTime(kgrid, medium.sound_speed);

% Find Focus from NIFTI File after padding
[~,I] = max(xdcr(:));
[l, m, n]= ind2sub(dim,I);
focus_pos = [l, m, n];

% check to make sure focus_pos is inside dim
if min(dim-focus_pos)<1||min(focus_pos)<1
    error('focus position is not in the grid')
end

%% Scale the element locations to the nifti file 

ROC_m = 72*1e-3;                    % [m] - radius of curvature
diameter_m = 6.6*1e-3;              % [m] - element diameter

ROC = round(ROC_m/vox(3));                          % [grid points]
diameter = round(diameter_m/vox(2));                % [grid points]
if ~mod(diameter,2) 
    diameter = diameter + 1;
end 

%% Define Transducer, Grid, Time 
    
kgrid = kWaveGrid(padsize(1),vox(1),padsize(2),vox(2),padsize(3),vox(3));

[kgrid.t_array, ~] = makeTime(kgrid, medium.sound_speed);
%% 7. Compute Amps and Phases 

kwavenum = 2*pi*freq/1500;

% inputs to function are in meters 
uamp = getFocusedElementVals(kwavenum, elemLocs.*1e-3, steerLoc.*1e-3, 1000);

uamp = uamp./max(abs(uamp)); % normalization bc IGT phases does a normalization
amps = abs(uamp);
phases = -1*angle(uamp); 

% 8. create pulses from ampps/phases 
% make a pulse - vec for each element
Ncyc = 5;                                           % num cycles
t = kgrid.t_array;                                  % time vector 
dt = kgrid.dt;                                      % time step 
nptspulse= round( Ncyc * 1 / (dt*freq) );           % npts in our pulse
tp = dt.*(1:1:nptspulse);                           % time vector for pulse

Ne = length(amps);                                  % num elements 
win = gausswin(nptspulse); 
pvec = zeros([Ne length(t)]);

% pvec = [num elements x pulse length in time] 
for i=1:Ne
    pulse = Amp*amps(i).*sin(2*pi*freq*tp+phases(i)); 
    pvec(i,1:nptspulse)=pulse.*win';
end

%% 8. Assign Pulses to Source Elements 
% lables follows order of the elem_pix 
% amps/phases follows order of the elem_pix mc
% so just match one-one

%pvec = [number of source elements x num of time points] 
% source.p = [source point index (with linear indexing) x time point index] 
% source.p_mask = binary mask of where there is and isn't a point
longLabels = niftiread(fnlabeledxdcr); 
longLabels = longLabels(:);
longLabels(longLabels==0)=[];

source.p = zeros([length(longLabels), length(t)]);
for i = 1 : Ne
    % Ne needs phase from sensor point where labels == i
    % this grabs indexes where labels == i 
    ind = find(longLabels == i);
    % repeats the signal through the length of longLabels 
    source.p(ind,:) = repmat(pvec(i,:),[length(ind),1]); 
 
end

xdcrMask = xdcr;
xdcrMask(xdcrMask>1)=0;
source.p_mask=xdcrMask; 

% 9. Define Sensor 

% adjust focus by steering amaount - focus is grid, steer must also be grid  
focus_pos = focus_pos + round((steerLoc.*1e-3)./vox(1));

% Records sensor data from everywhere or can choose to only record from the
% sensor at the focus location 
sensor.mask = ones(padsize);
% sensor.mask(focus_pos(1),:,:)=1;
% sensor.mask(:,focus_pos(2),:)=1;
% sensor.mask(:,:,focus_pos(3))=1;

sensor.record = {'p_max','p_rms'};
%% Run the Simulation

sensor_data = kspaceFirstOrder3DG(kgrid,medium,source,sensor);

  %% Reshape the sensor data 

pout = zeros([kgrid.Nx kgrid.Ny kgrid.Nz]);
curr=1;
for k=1:kgrid.Nz
    for j=1:kgrid.Ny
        for i=1:kgrid.Nx
            if sensor.mask(i,j,k)==1
                pmax(i,j,k)= sensor_data.p_max(curr);
                prms(i,j,k)=sensor_data.p_rms(curr); 
                curr=curr+1;
            end
        end
    end
end

yz = pmax(focus_pos(1),:,:);
xy = pmax(:,:,focus_pos(3));
xz = pmax(:,focus_pos(2),:);
pdata = [];
pdata.yz = yz;
pdata.xy = xy;
pdata.xz = xz;
%% 11. (Optional) Save the Data as a .mat file 

prms = prms(padpre(1)+1:end-padpost(1),padpre(2)+1:end-padpost(2),...
    padpre(3)+1:end-padpost(3));
pmax = pmax(padpre(1)+1:end-padpost(1),padpre(2)+1:end-padpost(2),...
    padpre(3)+1:end-padpost(3));


save(strcat(fnout,'.mat'),'focus_pos','focus_pos_init','Amp','steerLoc','prms','pmax');
%% Save to nifti

niftiwrite(pmax,fnout,'compressed',true);
%% 12. Extra Step - Visualization at the focus slice
figure;
sgtitle("Steered"); 
subplot(131)
imagesc(rot90(squeeze(pmax(focus_pos(1),:,:))))
xlabel("Y"); ylabel("Z"); 
axis image

subplot(132)
imagesc(rot90(squeeze(pmax(:,focus_pos(2),:))))
xlabel("X"); ylabel("Z");
axis image

subplot(133)
imagesc(rot90(squeeze(pmax(:,:,focus_pos(3)))))
xlabel("X"); ylabel("Y");
axis image












