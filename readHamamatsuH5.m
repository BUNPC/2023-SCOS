function images = readHamamatsuH5(varargin)
% readBaslerH5 reads Hamamatsu camera information in h5 file format
%   Inputs:
%       h5File = char array specifying h5 file full directory
%       timeInd = 2x1 vector specifying the frames to read. For example
%       [101 500] specifies reading time frame 101 to 500.

h5File = varargin{1};
if nargin > 1
    timeInd = varargin{2};
else
    nt = h5info(h5File,'/Images');
    nt = nt.Dataspace.Size;
    timeInd = [1 nt];
end

images = h5read(h5File,'/Images',[1 1 timeInd(1)],[Inf Inf timeInd(2) - timeInd(1) + 1]);
%images = h5read(h5File,'/Images',timeInd(1),timeInd(2) - timeInd(1) + 1);
images = imrotate(images,-90);
images = flip(images,2);
end

