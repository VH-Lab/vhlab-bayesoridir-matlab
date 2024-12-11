function b = angdiff(x)
flip = 0;
if numel(size(x))<3 && size(x,1)>size(x,2),x=x'; flip = 1; end;

b = min(abs([x;x+360;x-360])); 
if flip, b = b'; end;
end