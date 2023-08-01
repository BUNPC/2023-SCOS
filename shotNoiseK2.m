function [K2_shot] = shotNoiseK2(intensity,gain)
%correctShotNoise calculates shot noise K2 from intensity and camera gain
%(ADU/e-)
%   Inputs:
%       intensity = nx1 vector
%       gain (ADU/e-) = constant

offset = 1./intensity;
K2_shot = offset*gain;

end

