function amp = getInputAmplitude(p,varargin)
%getInputAmplitude
% Inputs:
% p - desired free-field pressure
%
% Outputs:
% amp - input amplitude to use at the transducer surface to get desired
% pressure
%
% OPTIONAL INPUTS: 
%   Optional 'string', value pairs that may be used to modify the default
%   computational settings. 
%   
%   'max'        - use max recorded pressure for calibration instead of RMS
%   (default is RMS pressure)
%   'steer'      - steering coordinates in xdcr space and mm 
% 
% Author: Michelle K. Sigona
% Edited: May 2, 2023

% Load saved values from calibration
load('InputAmplitudeCalibration.mat','x','max_pressmax','max_pressrms', ...
    'coord');

% Determine if using RMS pressure (default) or max pressure
ind = 0;
steer = [0,0,0];
steerFlag = 0;
calibration = max_pressrms; 

if nargin > 2
    for input_index = 1:2:length(varargin) 
        switch varargin{input_index}
            case 'max'
                calibration = max_pressmax;
            case 'steer'
                steerFlag = 1;
                steer = varargin{input_index+1};
        end
    end
end

% Find index for steered or geometric coordinates
if steerFlag
    for ci = 1:length(coord)
        if coord(ci,1) == steer(1) && coord(ci,2) == steer(2) && ...
                coord(ci,3) == steer(3)
            ind = ci; 
        end
    end
elseif ~steerFlag
    for ci = 1:length(coord)
        if coord(ci,1) == 0 && coord(ci,2) == 0 && coord(ci,3) == 0
            ind = ci; 
        end
    end
end

% If steered coordinates do not exist in the calibration 
if ~ind
    disp('Steered coordinates not available'); 
    return;
end

% Calculate linear fit
pf = polyfit(x,calibration(ind,:),1); 

% Calculate input amplitude required for desired pressure 
amp = (p*1e6+pf(2))/pf(1); 

end