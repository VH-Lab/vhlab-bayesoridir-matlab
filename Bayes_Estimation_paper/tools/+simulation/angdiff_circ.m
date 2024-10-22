function b = angdiff_circ(x)
flip = 0;
if size(x,1)>size(x,2),x=x'; flip = 1; end;

b = min(abs([x;x+2*pi;x-2*pi])); 
if flip, b = b'; end;
end