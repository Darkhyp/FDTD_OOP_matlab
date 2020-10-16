function [Sx,Hy,Ez] = Poynting(obj,X)
%#codegen

if nargin==1
    Sx = real(-obj.snapshotE('z0').*obj.snapshotH('y0'));
else
    Ez = interpn(obj.xEz,obj.Ez,X,'cubic');
    Hy = interpn(obj.xHy,obj.Hy,X,'cubic');
    Sx = real(-Ez.*Hy);
end

end
