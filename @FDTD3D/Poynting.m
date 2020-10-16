function [Sx,Sy,Sz,Ex,Ey,Ez,Hx,Hy,Hz] = Poynting(obj,X,Y,Z)
%#codegen

if nargin==1
	Ex = obj.snapshotE('x0');
	Ey = obj.snapshotE('y0');
	Ez = obj.snapshotE('z0');
	Hx = obj.snapshotH('x0');
	Hy = obj.snapshotH('y0');
	Hz = obj.snapshotH('z0');
else
	Ex = interpn(obj.xEx,obj.yEx,obj.zEx,obj.Ex,X,Y,Z,'cubic');
	Ey = interpn(obj.xEy,obj.yEy,obj.zEy,obj.Ey,X,Y,Z,'cubic');
	Ez = interpn(obj.xEz,obj.yEz,obj.zEz,obj.Ez,X,Y,Z,'cubic');
	Hx = interpn(obj.xHx,obj.yHx,obj.zHx,obj.Hx,X,Y,Z,'cubic');
	Hy = interpn(obj.xHy,obj.yHy,obj.zHy,obj.Hy,X,Y,Z,'cubic');
	Hz = interpn(obj.xHz,obj.yHz,obj.zHz,obj.Hz,X,Y,Z,'cubic');
end
Sx = real(Ey.*Hz-Ez.*Hy);
Sy = real(Ez.*Hx-Ex.*Hz);
Sz = real(Ex.*Hy-Ey.*Hx);

end
