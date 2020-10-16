%% Store field data
function obj = updateOUT(obj,time)
%#codegen

if isempty(obj.OUT)
    fprintf('Objects are not defined ...\n')
    return
end

% calculate Poynting
for n_obj=length(obj.OUT.Objects):-1:1
    switch obj.OUT.data{n_obj}.type
        case 1
            Norm = (max(obj.OUT.data{n_obj}.Y(:))-min(obj.OUT.data{n_obj}.Y(:)));
            X = obj.OUT.data{n_obj}.X(1);
            if obj.isTE
                [Sx,~,~,~,Ez,~,Hy,~] = Poynting(obj,obj.OUT.data{n_obj}.X,obj.OUT.data{n_obj}.Y);
                obj.OUT.data{n_obj}.Ez(:,time) = Ez;
                obj.OUT.data{n_obj}.Hy(:,time) = Hy;
                obj.OUT.data{n_obj}.Sx_aver(time) = integral(@(Y)(-interpn(obj.xEz,obj.yEz,obj.Ez,X*ones(size(Y)),Y,'cubic').*interpn(obj.xHy,obj.yHy,obj.Hy,X*ones(size(Y)),Y,'cubic')), ...
                    min(obj.OUT.data{n_obj}.Y(:)),max(obj.OUT.data{n_obj}.Y(:)) )/Norm;
            else
                [Sx,~,~,Ey,~,~,~,Hz] = Poynting(obj,obj.OUT.data{n_obj}.X,obj.OUT.data{n_obj}.Y);
                obj.OUT.data{n_obj}.Hz(:,time) = Hz;
                obj.OUT.data{n_obj}.Ey(:,time) = Ey;
                obj.OUT.data{n_obj}.Sx_aver(time) = integral(@(Y)( interpn(obj.xEy,obj.yEy,obj.Ey,X*ones(size(Y)),Y,'cubic').*interpn(obj.xHz,obj.yHz,obj.Hz,X*ones(size(Y)),Y,'cubic')), ...
                    min(obj.OUT.data{n_obj}.Y(:)),max(obj.OUT.data{n_obj}.Y(:)) )/Norm;
            end
            % plane perpendicular to x-axis (Sx)
            obj.OUT.data{n_obj}.Sx(:,time)	= Sx;
        case 2
            Norm = (max(obj.OUT.data{n_obj}.X(:))-min(obj.OUT.data{n_obj}.X(:)))*(max(obj.OUT.data{n_obj}.Z(:))-min(obj.OUT.data{n_obj}.Z(:)));
            Y = obj.OUT.data{n_obj}.Y(1);
            if obj.isTE
                [~,Sy,~,~,Ez,Hx,~,~] = Poynting(obj,obj.OUT.data{n_obj}.X,obj.OUT.data{n_obj}.Y,obj.OUT.data{n_obj}.Z);
                obj.OUT.data{n_obj}.Ez(:,time) = Ez;
                obj.OUT.data{n_obj}.Hx(:,time) = Hx;
                obj.OUT.data{n_obj}.Sy_aver(time) = integral(@(X)( interpn(obj.xEz,obj.yEz,obj.Ez,X,Y*ones(size(X)),'cubic').*interpn(obj.xHx,obj.yHx,obj.Hx,X,Y*ones(size(X)),'cubic')), ...
                    min(obj.OUT.data{n_obj}.X(:)),max(obj.OUT.data{n_obj}.X(:)) )/Norm;
            else
                [~,Sy,Ex,~,~,~,~,Hz] = Poynting(obj,obj.OUT.data{n_obj}.X,obj.OUT.data{n_obj}.Y,obj.OUT.data{n_obj}.Z);
                obj.OUT.data{n_obj}.Hz(:,time) = Hz;
                obj.OUT.data{n_obj}.Ex(:,time) = Ex;
                obj.OUT.data{n_obj}.Sy_aver(time) = integral(@(X)(-interpn(obj.xEx,obj.yEx,obj.Ex,X,Y*ones(size(X)),'cubic').*interpn(obj.xHz,obj.yHz,obj.Hz,X,Y*ones(size(X)),'cubic')), ...
                    min(obj.OUT.data{n_obj}.X(:)),max(obj.OUT.data{n_obj}.X(:)) )/Norm;
            end
            % plane perpendicular to y-axis (Sy)
            obj.OUT.data{n_obj}.Sy(:,time)	= Sy;
        case 4
            if obj.isTE
                % x-planes
                Sx = Poynting(obj,obj.OUT.data{n_obj}.X1x,obj.OUT.data{n_obj}.Y1x);
                obj.OUT.data{n_obj}.S1x(:,time)	= squeeze(Sx);
                Sx = Poynting(obj,obj.OUT.data{n_obj}.X2x,obj.OUT.data{n_obj}.Y2x);
                obj.OUT.data{n_obj}.S2x(:,time)	= squeeze(Sx);

                Norm = (max(obj.OUT.data{n_obj}.Y1x(:))-min(obj.OUT.data{n_obj}.Y1x(:)));
                X = obj.OUT.data{n_obj}.X1x(1);
                S_aver = integral(@(Y)(-interpn(obj.xEz,obj.yEz,obj.Ez,X*ones(size(Y)),Y,'cubic').*interpn(obj.xHy,obj.yHy,obj.Hy,X*ones(size(Y)),Y,'cubic')), ...
                    min(obj.OUT.data{n_obj}.Y1x(:)),max(obj.OUT.data{n_obj}.Y1x(:)) )/Norm;
                X = obj.OUT.data{n_obj}.X2x(1);
                S_aver = S_aver -integral(@(Y)(-interpn(obj.xEz,obj.yEz,obj.Ez,X*ones(size(Y)),Y,'cubic').*interpn(obj.xHy,obj.yHy,obj.Hy,X*ones(size(Y)),Y,'cubic')), ...
                    min(obj.OUT.data{n_obj}.Y1x(:)),max(obj.OUT.data{n_obj}.Y1x(:)) )/Norm;

                % y-planes
                [~,Sy] = Poynting(obj,obj.OUT.data{n_obj}.X1y,obj.OUT.data{n_obj}.Y1y);
                obj.OUT.data{n_obj}.S1y(:,:,time)	= squeeze(Sy);
                [~,Sy] = Poynting(obj,obj.OUT.data{n_obj}.X2y,obj.OUT.data{n_obj}.Y2y);
                obj.OUT.data{n_obj}.S2y(:,:,time)	= squeeze(Sy);

                Norm = (max(obj.OUT.data{n_obj}.X1y(:))-min(obj.OUT.data{n_obj}.X1y(:)));
                Y = obj.OUT.data{n_obj}.Y1y(1);
                S_aver = S_aver + integral(@(X)(interpn(obj.xEz,obj.yEz,obj.Ez,X,Y*ones(size(X)),'cubic').*interpn(obj.xHx,obj.yHx,obj.Hx,X,Y*ones(size(X)),'cubic')), ...
                    min(obj.OUT.data{n_obj}.X1y(:)),max(obj.OUT.data{n_obj}.X1y(:)) )/Norm;
                Y = obj.OUT.data{n_obj}.Y2y(1);
                S_aver = S_aver - integral(@(X)(interpn(obj.xEz,obj.yEz,obj.Ez,X,Y*ones(size(X)),'cubic').*interpn(obj.xHx,obj.yHx,obj.Hx,X,Y*ones(size(X)),'cubic')), ...
                    min(obj.OUT.data{n_obj}.X1y(:)),max(obj.OUT.data{n_obj}.X1y(:)) )/Norm;
            else
                % x-planes
                Sx = Poynting(obj,obj.OUT.data{n_obj}.X1x,obj.OUT.data{n_obj}.Y1x);
                obj.OUT.data{n_obj}.S1x(:,time)	= squeeze(Sx);
                Sx = Poynting(obj,obj.OUT.data{n_obj}.X2x,obj.OUT.data{n_obj}.Y2x);
                obj.OUT.data{n_obj}.S2x(:,time)	= squeeze(Sx);

                Norm = (max(obj.OUT.data{n_obj}.Y1x(:))-min(obj.OUT.data{n_obj}.Y1x(:)));
                X = obj.OUT.data{n_obj}.X1x(1);
                S_aver = integral(@(Y)(interpn(obj.xEy,obj.yEy,obj.Ey,X*ones(size(Y)),Y,'cubic').*interpn(obj.xHz,obj.yHz,obj.Hz,X*ones(size(Y)),Y,'cubic')), ...
                    min(obj.OUT.data{n_obj}.Y1x(:)),max(obj.OUT.data{n_obj}.Y1x(:)) )/Norm;
                X = obj.OUT.data{n_obj}.X2x(1);
                S_aver = S_aver -integral(@(Y)(interpn(obj.xEy,obj.yEy,obj.Ey,X*ones(size(Y)),Y,'cubic').*interpn(obj.xHz,obj.yHz,obj.Hz,X*ones(size(Y)),Y,'cubic')), ...
                    min(obj.OUT.data{n_obj}.Y1x(:)),max(obj.OUT.data{n_obj}.Y1x(:)) )/Norm;

                % y-planes
                [~,Sy] = Poynting(obj,obj.OUT.data{n_obj}.X1y,obj.OUT.data{n_obj}.Y1y);
                obj.OUT.data{n_obj}.S1y(:,:,time)	= squeeze(Sy);
                [~,Sy] = Poynting(obj,obj.OUT.data{n_obj}.X2y,obj.OUT.data{n_obj}.Y2y);
                obj.OUT.data{n_obj}.S2y(:,:,time)	= squeeze(Sy);

                Norm = (max(obj.OUT.data{n_obj}.X1y(:))-min(obj.OUT.data{n_obj}.X1y(:)));
                Y = obj.OUT.data{n_obj}.Y1y(1);
                S_aver = S_aver + integral(@(X)(-interpn(obj.xEx,obj.yEx,obj.Ex,X,Y*ones(size(X)),'cubic').*interpn(obj.xHz,obj.yHz,obj.Hz,X,Y*ones(size(X)),'cubic')), ...
                    min(obj.OUT.data{n_obj}.X1y(:)),max(obj.OUT.data{n_obj}.X1y(:)) )/Norm;
                Y = obj.OUT.data{n_obj}.Y2y(1);
                S_aver = S_aver - integral(@(X)(-interpn(obj.xEx,obj.yEx,obj.Ex,X,Y*ones(size(X)),'cubic').*interpn(obj.xHz,obj.yHz,obj.Hz,X,Y*ones(size(X)),'cubic')), ...
                    min(obj.OUT.data{n_obj}.X1y(:)),max(obj.OUT.data{n_obj}.X1y(:)) )/Norm;
            end
            obj.OUT.data{n_obj}.S_aver(time) = S_aver;
    end
end

end
