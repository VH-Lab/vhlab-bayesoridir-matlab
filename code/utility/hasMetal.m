function b = hasMetal()
% hasMetal - is this a Mac with MatlabMetal installed?
%
% b = hasMetal()
%
% Returns true if the computer is a Mac with the function 'MetalCallKernel' available,
% otherwise returns false.
% 

b = logical(ismac() && ~isempty(which('MetalCallKernel')));
