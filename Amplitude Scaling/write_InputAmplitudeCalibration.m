% Writes the InputAmplitudeCalibration.mat file 
% Outputs:
% coord - steering coordinate [x.x y.y z.z] 
% max_pressmax - pressure max to the corresponding coordinate 
% max_pressrms - max rms pressure to the corresponding coordinate 
% x - amplitude values used 
% written by Allie - 20230105 

clear all; close all; clc; 

files  = dir('pout*.mat');

for i = 1:length(files) 
    data(i) = load(files(i).name); 
end

%% Find and save values we care about

for i = 1:length(data); 
   rms_p = data(i).prms; 
   max_p = data(i).pmax; 
    
   max_maxp(i) = max(max_p(:));
   max_rmsp(i) = max(rms_p(:)); 
   
   amp(i) = data(i).Amp; 
   steerCoords(i,:) = data(i).steerLoc;
end 

%% Find the indices of steerCoords where the rows match (i.e. same steer coordinate)

coord = unique(steerCoords,'rows');
x = unique(amp);

for i = 1:length(coord)
    steerLoc = coord(i,:); 
    ind = find(ismember(steerCoords,steerLoc,'rows')); 
    
    max_pressrms(i,:) = [max_rmsp(ind(1)), max_rmsp(ind(2))]; 
    max_pressmax(i,:) = [max_maxp(ind(1)), max_maxp(ind(2))];      
end 

%% Now save everything into the .mat file 

save('InputAmplitudeCalibration.mat','x','coord','max_pressrms','max_pressmax'); 
