function [Sx,Sy,Ex,Ey,Ez,Hx,Hy,Hz] = Poynting(obj,X,Y)
%#codegen

Ex = [];
Ey = [];
Ez = [];
Hx = [];
Hy = [];
Hz = [];
if nargin==1
    if obj.isTE
        Ez = obj.snapshotE('z0');
        Hx = obj.snapshotH('x0');
        Hy = obj.snapshotH('y0');
        Sx = -real(Ez.*Hy);
        Sy =  real(Ez.*Hx);
    else
        Ex = obj.snapshotE('x0');
        Ey = obj.snapshotE('y0');
        Hz = obj.snapshotH('z0');
        Sx =  real(Ey.*Hz);
        Sy = -real(Ex.*Hz);
    end
else
    if obj.isTE
        Ez = interpn(obj.xEz,obj.yEz,obj.Ez,X,Y,'cubic');
        Hx = interpn(obj.xHx,obj.yHx,obj.Hx,X,Y,'cubic');
        Hy = interpn(obj.xHy,obj.yHy,obj.Hy,X,Y,'cubic');
        Sx = -real(Ez.*Hy);
        Sy =  real(Ez.*Hx);
    else
        Ex = interpn(obj.xEx,obj.yEx,obj.Ex,X,Y,'cubic');
        Ey = interpn(obj.xEy,obj.yEy,obj.Ey,X,Y,'cubic');
        Hz = interpn(obj.xHz,obj.yHz,obj.Hz,X,Y,'cubic');
        Sx =  real(Ey.*Hz);
        Sy = -real(Ex.*Hz);
    end

end
