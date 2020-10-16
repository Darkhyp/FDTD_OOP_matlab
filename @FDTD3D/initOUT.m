%% initialize out field data
function obj = initOUT(obj,Objects, MaxTime)
%#codegen

if ~isempty(Objects)
	% check objects
	obj.OUT.Objects = Objects;
	for n_obj=length(obj.OUT.Objects):-1:1
        if obj.PMLtype>0
            ind = obj.OUT.Objects{n_obj}(1:3,:) + obj.PML.layers;
        else
            ind = obj.OUT.Objects{n_obj}(1:3,:);
        end

		tmp = diff(ind,[],2);

		if any(tmp==0)
			obj.OUT.data{n_obj}.X = obj.x(ind(1,1):ind(1,2),ind(2,1):ind(2,2),ind(3,1):ind(3,2));
			obj.OUT.data{n_obj}.Y = obj.y(ind(1,1):ind(1,2),ind(2,1):ind(2,2),ind(3,1):ind(3,2));
			obj.OUT.data{n_obj}.Z = obj.z(ind(1,1):ind(1,2),ind(2,1):ind(2,2),ind(3,1):ind(3,2));
		end

		switch true
			case tmp(1)==0
				obj.OUT.data{n_obj}.type	= 1;
				obj.OUT.data{n_obj}.Sx      = zeros([tmp(tmp~=0)+1;MaxTime]');
                obj.OUT.data{n_obj}.Hy      = obj.OUT.data{n_obj}.Sx;
                obj.OUT.data{n_obj}.Hz      = obj.OUT.data{n_obj}.Sx;
                obj.OUT.data{n_obj}.Ey      = obj.OUT.data{n_obj}.Sx;
                obj.OUT.data{n_obj}.Ez      = obj.OUT.data{n_obj}.Sx;
			case tmp(2)==0
				obj.OUT.data{n_obj}.type	= 2;
				obj.OUT.data{n_obj}.Sy      = zeros([tmp(tmp~=0)+1;MaxTime]');
                obj.OUT.data{n_obj}.Hx      = obj.OUT.data{n_obj}.Sy;
                obj.OUT.data{n_obj}.Hz      = obj.OUT.data{n_obj}.Sy;
                obj.OUT.data{n_obj}.Ex      = obj.OUT.data{n_obj}.Sy;
                obj.OUT.data{n_obj}.Ez      = obj.OUT.data{n_obj}.Sy;
			case tmp(3)==0
				obj.OUT.data{n_obj}.type	= 3;
				obj.OUT.data{n_obj}.Sz      = zeros([tmp(tmp~=0)+1;MaxTime]');
                obj.OUT.data{n_obj}.Hx      = obj.OUT.data{n_obj}.Sz;
                obj.OUT.data{n_obj}.Hy      = obj.OUT.data{n_obj}.Sz;
                obj.OUT.data{n_obj}.Ex      = obj.OUT.data{n_obj}.Sz;
                obj.OUT.data{n_obj}.Ey      = obj.OUT.data{n_obj}.Sz;
			otherwise
				obj.OUT.data{n_obj}.type	= 4;

				obj.OUT.data{n_obj}.X1x = obj.x(ind(1,1),ind(2,1):ind(2,2),ind(3,1):ind(3,2));
				obj.OUT.data{n_obj}.Y1x = obj.y(ind(1,1),ind(2,1):ind(2,2),ind(3,1):ind(3,2));
				obj.OUT.data{n_obj}.Z1x = obj.z(ind(1,1),ind(2,1):ind(2,2),ind(3,1):ind(3,2));
				obj.OUT.data{n_obj}.X2x = obj.x(ind(1,2),ind(2,1):ind(2,2),ind(3,1):ind(3,2));
				obj.OUT.data{n_obj}.Y2x = obj.y(ind(1,2),ind(2,1):ind(2,2),ind(3,1):ind(3,2));
				obj.OUT.data{n_obj}.Z2x = obj.z(ind(1,2),ind(2,1):ind(2,2),ind(3,1):ind(3,2));
				obj.OUT.data{n_obj}.S1x	= zeros([tmp(2:3)+1;MaxTime]');
				obj.OUT.data{n_obj}.S2x	= zeros([tmp(2:3)+1;MaxTime]');

				obj.OUT.data{n_obj}.X1y = obj.x(ind(1,1):ind(1,2),ind(2,1),ind(3,1):ind(3,2));
				obj.OUT.data{n_obj}.Y1y = obj.y(ind(1,1):ind(1,2),ind(2,1),ind(3,1):ind(3,2));
				obj.OUT.data{n_obj}.Z1y = obj.z(ind(1,1):ind(1,2),ind(2,1),ind(3,1):ind(3,2));
				obj.OUT.data{n_obj}.X2y = obj.x(ind(1,1):ind(1,2),ind(2,2),ind(3,1):ind(3,2));
				obj.OUT.data{n_obj}.Y2y = obj.y(ind(1,1):ind(1,2),ind(2,2),ind(3,1):ind(3,2));
				obj.OUT.data{n_obj}.Z2y = obj.z(ind(1,1):ind(1,2),ind(2,2),ind(3,1):ind(3,2));
				obj.OUT.data{n_obj}.S1y	= zeros([tmp([1,3])+1;MaxTime]');
				obj.OUT.data{n_obj}.S2y	= zeros([tmp([1,3])+1;MaxTime]');

				obj.OUT.data{n_obj}.X1z = obj.x(ind(1,1):ind(1,2),ind(2,1):ind(2,2),ind(3,1));
				obj.OUT.data{n_obj}.Y1z = obj.y(ind(1,1):ind(1,2),ind(2,1):ind(2,2),ind(3,1));
				obj.OUT.data{n_obj}.Z1z = obj.z(ind(1,1):ind(1,2),ind(2,1):ind(2,2),ind(3,1));
				obj.OUT.data{n_obj}.X2z = obj.x(ind(1,1):ind(1,2),ind(2,1):ind(2,2),ind(3,2));
				obj.OUT.data{n_obj}.Y2z = obj.y(ind(1,1):ind(1,2),ind(2,1):ind(2,2),ind(3,2));
				obj.OUT.data{n_obj}.Z2z = obj.z(ind(1,1):ind(1,2),ind(2,1):ind(2,2),ind(3,2));
				obj.OUT.data{n_obj}.S1z	= zeros([tmp(1:2)+1;MaxTime]');
				obj.OUT.data{n_obj}.S2z	= zeros([tmp(1:2)+1;MaxTime]');
		end
	end
%{
	for n_obj=1:length(Objects)
		otmp = Objects{n_obj};
		switch otmp.type
			case 'point'
			case 'plane'
				condition = otmp.normal(1)*(x-otmp.position(1))+otmp.normal(2)*(y-otmp.position(2))+otmp.normal(3)*(z-otmp.position(3))==0;

				Ex = @(X,Y,Z) interp3(obj.xEx,obj.yEx,obj.zEx,obj.Ex,X,Y,Z,'cubic');
				Ey = @(X,Y,Z) interp3(obj.xEy,obj.yEy,obj.zEy,obj.Ey,X,Y,Z,'cubic');
				Ez = @(X,Y,Z) interp3(obj.xEz,obj.yEz,obj.zEz,obj.Ez,X,Y,Z,'cubic');
				Hx = @(X,Y,Z) interp3(obj.xHx,obj.YHx,obj.zHx,obj.Hx,X,Y,Z,'cubic');
				Hy = @(X,Y,Z) interp3(obj.xHy,obj.yHy,obj.zHy,obj.Hy,X,Y,Z,'cubic');
				Hz = @(X,Y,Z) interp3(obj.xHz,obj.yHz,obj.zHz,obj.Hz,X,Y,Z,'cubic');


		end
	end
%}                            
end

end

