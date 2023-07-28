function [K2_shot] = shotNoiseK2(intensity,gain)
%correctShotNoise Summary of this function goes here
%   speckleContrastSquared
%   intensity
%   gain (ADU/e-)

offset = 1./intensity;
K2_shot = offset*gain;

end

