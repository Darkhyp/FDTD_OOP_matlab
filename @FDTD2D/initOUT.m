%% initialize out field data
function obj = initOUT(obj,Objects, MaxTime)
%#codegen

if ~isempty(Objects)
	% check objects
	obj.OUT.Objects = Objects;
	for n_obj=length(obj.OUT.Objects):-1:1
        if obj.PMLtype>0
            ind = obj.OUT.Objects{n_obj}(1:2,:) + obj.PML.layers;
        else
            ind = obj.OUT.Objects{n_obj}(1:2,:);
        end

		tmp = diff(ind(1:2,:),[],2);

		if any(tmp==0)
			obj.OUT.data{n_obj}.X = obj.x(ind(1,1):ind(1,2),ind(2,1):ind(2,2));
			obj.OUT.data{n_obj}.Y = obj.y(ind(1,1):ind(1,2),ind(2,1):ind(2,2));
		end

		switch true
			case tmp(1)==0
				obj.OUT.data{n_obj}.type	= 1;
                obj.OUT.data{n_obj}.Sx      = zeros([tmp(tmp~=0)+1;MaxTime]');
                if obj.isTE
                    obj.OUT.data{n_obj}.Hx      = obj.OUT.data{n_obj}.Sx;
                    obj.OUT.data{n_obj}.Hy      = obj.OUT.data{n_obj}.Sx;
                    obj.OUT.data{n_obj}.Ez      = obj.OUT.data{n_obj}.Sx;
                else
                    obj.OUT.data{n_obj}.Hz      = obj.OUT.data{n_obj}.Sx;
                    obj.OUT.data{n_obj}.Ex      = obj.OUT.data{n_obj}.Sx;
                    obj.OUT.data{n_obj}.Ey      = obj.OUT.data{n_obj}.Sx;
                end
			case tmp(2)==0
				obj.OUT.data{n_obj}.type	= 2;
				obj.OUT.data{n_obj}.Sy      = zeros([tmp(tmp~=0)+1;MaxTime]');
                if obj.isTE
                    obj.OUT.data{n_obj}.Hx      = obj.OUT.data{n_obj}.Sy;
                    obj.OUT.data{n_obj}.Hy      = obj.OUT.data{n_obj}.Sy;
                    obj.OUT.data{n_obj}.Ez      = obj.OUT.data{n_obj}.Sy;
                else
                    obj.OUT.data{n_obj}.Hz      = obj.OUT.data{n_obj}.Sy;
                    obj.OUT.data{n_obj}.Ex      = obj.OUT.data{n_obj}.Sy;
                    obj.OUT.data{n_obj}.Ey      = obj.OUT.data{n_obj}.Sy;
                end
			otherwise
				obj.OUT.data{n_obj}.type	= 4;

				obj.OUT.data{n_obj}.X1x = obj.x(ind(1,1),ind(2,1):ind(2,2));
				obj.OUT.data{n_obj}.Y1x = obj.y(ind(1,1),ind(2,1):ind(2,2));
				obj.OUT.data{n_obj}.X2x = obj.x(ind(1,2),ind(2,1):ind(2,2));
				obj.OUT.data{n_obj}.Y2x = obj.y(ind(1,2),ind(2,1):ind(2,2));
				obj.OUT.data{n_obj}.S1x	= zeros([tmp(2)+1;MaxTime]');
				obj.OUT.data{n_obj}.S2x	= zeros([tmp(2)+1;MaxTime]');

				obj.OUT.data{n_obj}.X1y = obj.x(ind(1,1):ind(1,2),ind(2,1));
				obj.OUT.data{n_obj}.Y1y = obj.y(ind(1,1):ind(1,2),ind(2,1));
				obj.OUT.data{n_obj}.X2y = obj.x(ind(1,1):ind(1,2),ind(2,2));
				obj.OUT.data{n_obj}.Y2y = obj.y(ind(1,1):ind(1,2),ind(2,2));
				obj.OUT.data{n_obj}.S1y	= zeros([tmp(1)+1;MaxTime]');
				obj.OUT.data{n_obj}.S2y	= zeros([tmp(1)+1;MaxTime]');
		end
	end
end

end

