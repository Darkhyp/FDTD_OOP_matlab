function [Sx,Sy,Sz] = Poynting(obj,X,Y,Z)
    if nargin==1
        Sx = real(obj.snapshotE('y0').*obj.snapshotH('z0')-obj.snapshotE('z0').*obj.snapshotH('y0'));
        Sy = real(obj.snapshotE('z0').*obj.snapshotH('x0')-obj.snapshotE('x0').*obj.snapshotH('z0'));
        Sz = real(obj.snapshotE('x0').*obj.snapshotH('y0')-obj.snapshotE('y0').*obj.snapshotH('x0'));
    else
        tEx = interpn(obj.xEx,obj.yEx,obj.zEx,obj.Ex,X,Y,Z,'cubic');
        tEy = interpn(obj.xEy,obj.yEy,obj.zEy,obj.Ey,X,Y,Z,'cubic');
        tEz = interpn(obj.xEz,obj.yEz,obj.zEz,obj.Ez,X,Y,Z,'cubic');
        tHx = interpn(obj.xHx,obj.yHx,obj.zHx,obj.Hx,X,Y,Z,'cubic');
        tHy = interpn(obj.xHy,obj.yHy,obj.zHy,obj.Hy,X,Y,Z,'cubic');
        tHz = interpn(obj.xHz,obj.yHz,obj.zHz,obj.Hz,X,Y,Z,'cubic');
        Sx = real(tEy.*tHz-tEz.*tHy);
        Sy = real(tEz.*tHx-tEx.*tHz);
        Sz = real(tEx.*tHy-tEy.*tHx);


    end
end
