% Slicer2Kwave - Transcranial Multielement Steering
% Written by Michelle Sigona, Edited by Allison Dockum 
% Last Updated - 01/17/2024

clear all; close all; clc; 
%% Load in File Names - files are outputs from Slicer 

% labelled xdcr  - with number assigned to the transducer elements
% carries through the transforms applied in Slicer 
fnlabeledxdcr = 'IGT_labeled.nii.gz'; 

% unlabelled xdcr file w/ focus
fnxdcr = 'IGT.nii.gz'; 

% element locations [mm] - used for phasing array elements
load('IGTelemLocs.mat'); 

% resampled CT from slicer
fnct = 'myCT_resampled.nii.gz'; 

% Output filename
fnout = 'pout_myCT'; 

 %% Transducer Setting 

freq = 650e3;                          % [Hertz]
ffPressure = 1.68;                     % [MPa] 

%steerLoc is in US-space
steerLoc = [0 0 -5];                   % [mm] - [X Y Z]

% amplitude is at element surface - scaled based on estimate free-field
% pressure 
Amp = getInputAmplitude(ffPressure,'steer',steerLoc,'max');  
%% Read in required files

ct.raw_data = double(niftiread(fnct)); 
ct.info = niftiinfo(fnct); 
xdcr = niftiread(fnxdcr); 
labeled_xdcr = niftiread(fnlabeledxdcr); 

% Rescale values
ct.raw_data = ct.raw_data + ct.info.AdditiveOffset;
%% Create computational grid 

% Resize CT to match desired computational grid 
vox = [0.25 0.25 0.25]*1e-3;        % Make the grid isotropic

% Check Prime Factors
dim = size(xdcr);
'Your simpace size'
dim
%checkFactors(min(dim),max(dim)+100)

[~,I] = max(xdcr(:));
[l, m, n]= ind2sub(dim,I);
focus_pos_init = [l, m, n];
%% Pad CT to prime factors 

padsize = [432,432,432]; 

padpre = round((padsize-dim)/2);
padpre(3) = 0;
padpost = floor((padsize-dim)/2);
padpost(3) = round(padsize(3)-dim(3)); 

ct.data = padarray(ct.raw_data,padpre,0,'pre');
ct.data = padarray(ct.data,padpost,0,'post');

xdcr = padarray(xdcr,padpre,0,'pre'); 
xdcr = padarray(xdcr,padpost,0,'post'); 
%% Create medium 
% Threshold CT
boneLow = 400;            % used to set background and tissue as same value
boneHigh = 2000;          % used to set teeth or any implants as same value 

CT = ct.data; 
CT(CT>boneHigh) = boneHigh; 
CT(CT<boneLow) = 0;

% Calculate bone porosity 
CT_bone = (1-(CT/max(CT(:)))); 
bone_i = CT > boneLow; 

% Calculate density
rho_water = 1000;      % [kg/m^3]
rho_bone =  2200;      % [kg/m^3]
medium.density = rho_water*CT_bone + rho_bone.*(1-CT_bone);

% Calculate speed of sound
c_water = 1500;        % [m/s]
c_bone = 3100;         % [m/s]
medium.sound_speed= c_water*CT_bone + c_bone.*(1-CT_bone);

beta = 0.5;
medium.alpha_power = 1.1;

% from aubry and commonly used              % [dB/MHz/cm]
abs_water = 0;                              % [dB/MHz/cm]
abs_min = 0.2;                              % [dB/MHz/cm]
abs_max = 8;                                % [dB/MHz/cm]
medium.alpha_coeff = abs_water.*ones(size(CT));
medium.alpha_coeff(bone_i) = abs_min + (abs_max-abs_min).* ...
    (CT_bone(bone_i)).^beta;

%% Make kgrid and time

kgrid = kWaveGrid(padsize(1),vox(1),padsize(2),vox(2),padsize(3),vox(3));
[kgrid.t_array, ~] = makeTime(kgrid, medium.sound_speed);


% focus after padding 
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
    % this grabs indexes where labels == i 
    
kgrid = kWaveGrid(padsize(1),vox(1),padsize(2),vox(2),padsize(3),vox(3));

% define time too 
[kgrid.t_array, ~] = makeTime(kgrid, medium.sound_speed);
%% 7. Compute Amps and Phases 

kwavenum = 2*pi*freq/1500;

% inputs to function are in meters 
uamp = getFocusedElementVals(kwavenum, elemLocs.*1e-3, steerLoc.*1e-3, 1000);

% normalize the amplitudes because internal IGT-system does a normalization
% when it computes its pulses 
uamp = uamp./max(abs(uamp)); 
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
% amps/phases follows order of the elem_pix 
% so just match one-one

% pvec = [number of source elements x num of time points] 
% source.p = [source point index (with linear indexing) x time point index] 
% source.p_mask = binary mask of where there is and isn't an element
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

%adjust focus by steering amaount - focus is grid, steer must also be grid  
focus_pos = focus_pos + round((steerLoc.*1e-3)./vox(1));

% Records sensor data from everywhere, can choose to only record from the
% sensor at the focus location 
sensor.mask = ones(padsize);
% sensor.mask(focus_pos(1),:,:)=1;
% sensor.mask(:,focus_pos(2),:)=1;
% sensor.mask(:,:,focus_pos(3))=1;

% can select to record max pressure or rms pressure 
sensor.record = {'p_max'};
%% Run the Simulation

sensor_data = kspaceFirstOrder3DG(kgrid,medium,source,sensor);

%% Reshape the sensor data 

pmax = zeros([kgrid.Nx kgrid.Ny kgrid.Nz]);
curr=1;
for k=1:kgrid.Nz
    for j=1:kgrid.Ny
        for i=1:kgrid.Nx
            if sensor.mask(i,j,k)==1
                pmax(i,j,k)= sensor_data.p_max(curr);
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

pmax = pmax(padpre(1)+1:end-padpost(1),padpre(2)+1:end-padpost(2),...
    padpre(3)+1:end-padpost(3));

save(strcat(fnout,'.mat'),'focus_pos','focus_pos_init','Amp','steerLoc','pmax'); 

%% Unpad the pressure data to fit the original CT Size and save to nifti

pmax = pmax(padpre(1)+1:end-padpost(1),padpre(2)+1:end-padpost(2),...
    padpre(3)+1:end-padpost(3));

niftiwrite(pmax,fnout,'compressed',true);
%% 12. Extra Step - Visualization at the focus slice
figure;
sgtitle("Steered"); 
subplot(131)
p=pmax;
imagesc(rot90(squeeze(p(focus_pos(1),:,:))))
xlabel("Y"); ylabel("Z"); 
axis image

subplot(132)
imagesc(rot90(squeeze(p(:,focus_pos(2),:))))
xlabel("X"); ylabel("Z");
axis image

subplot(133)
imagesc(rot90(squeeze(p(:,:,focus_pos(3)))))
xlabel("X"); ylabel("Y");
axis image












