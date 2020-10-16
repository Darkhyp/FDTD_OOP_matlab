%% Store field data
function obj = updateOUT(obj,time)
    if isempty(obj.OUT)
        fprintf('Objects are not defined ...\n')
        return
    end
    % calculate Poynting
    for n_obj=length(obj.OUT.Objects):-1:1
        switch obj.OUT.data{n_obj}.type
            case 1
                Sx = Poynting(obj,obj.OUT.data{n_obj}.X,obj.OUT.data{n_obj}.Y,obj.OUT.data{n_obj}.Z);
                % plane perpendicular to x-axis (Sx)
                obj.OUT.data{n_obj}.Sx(:,:,time)	= squeeze(Sx);

                Norm = (max(obj.OUT.data{n_obj}.Y(:))-min(obj.OUT.data{n_obj}.Y(:)))*(max(obj.OUT.data{n_obj}.Z(:))-min(obj.OUT.data{n_obj}.Z(:)));
                X = obj.OUT.data{n_obj}.X(1);
                obj.OUT.data{n_obj}.Sx_aver(time) = integral2(@(Y,Z)(interpn(obj.xEy,obj.yEy,obj.zEy,obj.Ey,X*ones(size(Y)),Y,Z,'cubic').*interpn(obj.xHz,obj.yHz,obj.zHz,obj.Hz,X*ones(size(Y)),Y,Z,'cubic') ...
                                                                    -interpn(obj.xEz,obj.yEz,obj.zEz,obj.Ez,X*ones(size(Y)),Y,Z,'cubic').*interpn(obj.xHy,obj.yHy,obj.zHy,obj.Hy,X*ones(size(Y)),Y,Z,'cubic')), ...
                    min(obj.OUT.data{n_obj}.Y(:)),max(obj.OUT.data{n_obj}.Y(:)), min(obj.OUT.data{n_obj}.Z(:)),max(obj.OUT.data{n_obj}.Z(:)))/Norm;
            case 2
               [~,Sy] = Poynting(obj,obj.OUT.data{n_obj}.X,obj.OUT.data{n_obj}.Y,obj.OUT.data{n_obj}.Z);
                % plane perpendicular to y-axis (Sy)
                obj.OUT.data{n_obj}.Sy(:,:,time)	= squeeze(Sy);

                Norm = (max(obj.OUT.data{n_obj}.X(:))-min(obj.OUT.data{n_obj}.X(:)))*(max(obj.OUT.data{n_obj}.Z(:))-min(obj.OUT.data{n_obj}.Z(:)));
                Y = obj.OUT.data{n_obj}.Y(1);
                obj.OUT.data{n_obj}.Sy_aver(time) = integral2(@(X,Z)(interpn(obj.xEz,obj.yEz,obj.zEz,obj.Ez,X,Y*ones(size(X)),Z,'cubic').*interpn(obj.xHx,obj.yHx,obj.zHx,obj.Hx,X,Y*ones(size(X)),Z,'cubic') ...
                                                                    -interpn(obj.xEx,obj.yEx,obj.zEx,obj.Ex,X,Y*ones(size(X)),Z,'cubic').*interpn(obj.xHz,obj.yHz,obj.zHz,obj.Hz,X,Y*ones(size(X)),Z,'cubic')), ...
                    min(obj.OUT.data{n_obj}.X(:)),max(obj.OUT.data{n_obj}.X(:)), min(obj.OUT.data{n_obj}.Z(:)),max(obj.OUT.data{n_obj}.Z(:)))/Norm;

            case 3
                [~,~,Sz] = Poynting(obj,obj.OUT.data{n_obj}.X,obj.OUT.data{n_obj}.Y,obj.OUT.data{n_obj}.Z);
                % plane perpendicular to z-axis (Sz)
                obj.OUT.data{n_obj}.Sz(:,:,time)	= Sz;

                Norm = (max(obj.OUT.data{n_obj}.X(:))-min(obj.OUT.data{n_obj}.X(:)))*(max(obj.OUT.data{n_obj}.Y(:))-min(obj.OUT.data{n_obj}.Y(:)));
                Z = obj.OUT.data{n_obj}.Z(1);
                obj.OUT.data{n_obj}.Sz_aver(time) = integral2(@(X,Y)(interpn(obj.xEx,obj.yEx,obj.zEx,obj.Ex,X,Y,Z*ones(size(X)),'cubic').*interpn(obj.xHy,obj.yHy,obj.zHy,obj.Hy,X,Y,Z*ones(size(X)),'cubic') ...
                                                                    -interpn(obj.xEy,obj.yEy,obj.zEy,obj.Ey,X,Y,Z*ones(size(X)),'cubic').*interpn(obj.xHx,obj.yHx,obj.zHx,obj.Hx,X,Y,Z*ones(size(X)),'cubic')), ...
                    min(obj.OUT.data{n_obj}.X(:)),max(obj.OUT.data{n_obj}.X(:)), min(obj.OUT.data{n_obj}.Y(:)),max(obj.OUT.data{n_obj}.Y(:)))/Norm;
            case 4
                % x-planes
                Sx = Poynting(obj,obj.OUT.data{n_obj}.X1x,obj.OUT.data{n_obj}.Y1x,obj.OUT.data{n_obj}.Z1x);
                obj.OUT.data{n_obj}.S1x(:,:,time)	= squeeze(Sx);
                Sx = Poynting(obj,obj.OUT.data{n_obj}.X2x,obj.OUT.data{n_obj}.Y2x,obj.OUT.data{n_obj}.Z2x);
                obj.OUT.data{n_obj}.S2x(:,:,time)	= squeeze(Sx);

                Norm = (max(obj.OUT.data{n_obj}.Y1x(:))-min(obj.OUT.data{n_obj}.Y1x(:)))*(max(obj.OUT.data{n_obj}.Z1x(:))-min(obj.OUT.data{n_obj}.Z1x(:)));
                X = obj.OUT.data{n_obj}.X1x(1);
                S_aver = integral2(@(Y,Z)(interpn(obj.xEy,obj.yEy,obj.zEy,obj.Ey,X*ones(size(Y)),Y,Z,'cubic').*interpn(obj.xHz,obj.yHz,obj.zHz,obj.Hz,X*ones(size(Y)),Y,Z,'cubic') ...
                                         -interpn(obj.xEz,obj.yEz,obj.zEz,obj.Ez,X*ones(size(Y)),Y,Z,'cubic').*interpn(obj.xHy,obj.yHy,obj.zHy,obj.Hy,X*ones(size(Y)),Y,Z,'cubic')), ...
                    min(obj.OUT.data{n_obj}.Y1x(:)),max(obj.OUT.data{n_obj}.Y1x(:)), min(obj.OUT.data{n_obj}.Z1x(:)),max(obj.OUT.data{n_obj}.Z1x(:)))/Norm;
                X = obj.OUT.data{n_obj}.X2x(1);
                S_aver = S_aver -integral2(@(Y,Z)(interpn(obj.xEy,obj.yEy,obj.zEy,obj.Ey,X*ones(size(Y)),Y,Z,'cubic').*interpn(obj.xHz,obj.yHz,obj.zHz,obj.Hz,X*ones(size(Y)),Y,Z,'cubic') ...
                                                 -interpn(obj.xEz,obj.yEz,obj.zEz,obj.Ez,X*ones(size(Y)),Y,Z,'cubic').*interpn(obj.xHy,obj.yHy,obj.zHy,obj.Hy,X*ones(size(Y)),Y,Z,'cubic')), ...
                    min(obj.OUT.data{n_obj}.Y1x(:)),max(obj.OUT.data{n_obj}.Y1x(:)), min(obj.OUT.data{n_obj}.Z1x(:)),max(obj.OUT.data{n_obj}.Z1x(:)))/Norm;

                % y-planes
                [~,Sy] = Poynting(obj,obj.OUT.data{n_obj}.X1y,obj.OUT.data{n_obj}.Y1y,obj.OUT.data{n_obj}.Z1y);
                obj.OUT.data{n_obj}.S1y(:,:,time)	= squeeze(Sy);
                [~,Sy] = Poynting(obj,obj.OUT.data{n_obj}.X2y,obj.OUT.data{n_obj}.Y2y,obj.OUT.data{n_obj}.Z2y);
                obj.OUT.data{n_obj}.S2y(:,:,time)	= squeeze(Sy);

                Norm = (max(obj.OUT.data{n_obj}.X1y(:))-min(obj.OUT.data{n_obj}.X1y(:)))*(max(obj.OUT.data{n_obj}.Z1y(:))-min(obj.OUT.data{n_obj}.Z1y(:)));
                Y = obj.OUT.data{n_obj}.Y1y(1);
                S_aver = S_aver + integral2(@(X,Z)(interpn(obj.xEy,obj.yEy,obj.zEy,obj.Ey,X,Y*ones(size(X)),Z,'cubic').*interpn(obj.xHz,obj.yHz,obj.zHz,obj.Hz,X,Y*ones(size(X)),Z,'cubic') ...
                                                  -interpn(obj.xEz,obj.yEz,obj.zEz,obj.Ez,X,Y*ones(size(X)),Z,'cubic').*interpn(obj.xHy,obj.yHy,obj.zHy,obj.Hy,X,Y*ones(size(X)),Z,'cubic')), ...
                    min(obj.OUT.data{n_obj}.X1y(:)),max(obj.OUT.data{n_obj}.X1y(:)), min(obj.OUT.data{n_obj}.Z1y(:)),max(obj.OUT.data{n_obj}.Z1y(:)))/Norm;
                Y = obj.OUT.data{n_obj}.Y2y(1);
                S_aver = S_aver - integral2(@(X,Z)(interpn(obj.xEy,obj.yEy,obj.zEy,obj.Ey,X,Y*ones(size(X)),Z,'cubic').*interpn(obj.xHz,obj.yHz,obj.zHz,obj.Hz,X,Y*ones(size(X)),Z,'cubic') ...
                                                  -interpn(obj.xEz,obj.yEz,obj.zEz,obj.Ez,X,Y*ones(size(X)),Z,'cubic').*interpn(obj.xHy,obj.yHy,obj.zHy,obj.Hy,X,Y*ones(size(X)),Z,'cubic')), ...
                    min(obj.OUT.data{n_obj}.X1y(:)),max(obj.OUT.data{n_obj}.X1y(:)), min(obj.OUT.data{n_obj}.Z1y(:)),max(obj.OUT.data{n_obj}.Z1y(:)))/Norm;

                % z-planes
                [~,~,Sz] = Poynting(obj,obj.OUT.data{n_obj}.X1z,obj.OUT.data{n_obj}.Y1z,obj.OUT.data{n_obj}.Z1z);
                obj.OUT.data{n_obj}.S1z(:,:,time)	= Sz;
                [~,~,Sz] = Poynting(obj,obj.OUT.data{n_obj}.X2z,obj.OUT.data{n_obj}.Y2z,obj.OUT.data{n_obj}.Z2z);
                obj.OUT.data{n_obj}.S2z(:,:,time)	= Sz;

                Norm = (max(obj.OUT.data{n_obj}.X1z(:))-min(obj.OUT.data{n_obj}.X1z(:)))*(max(obj.OUT.data{n_obj}.Y1z(:))-min(obj.OUT.data{n_obj}.Y1z(:)));
                Z = obj.OUT.data{n_obj}.Z1z(1);
                S_aver = S_aver + integral2(@(X,Y)(interpn(obj.xEy,obj.yEy,obj.zEy,obj.Ey,X,Y,Z*ones(size(X)),'cubic').*interpn(obj.xHz,obj.yHz,obj.zHz,obj.Hz,X,Y,Z*ones(size(X)),'cubic') ...
                                                  -interpn(obj.xEz,obj.yEz,obj.zEz,obj.Ez,X,Y,Z*ones(size(X)),'cubic').*interpn(obj.xHy,obj.yHy,obj.zHy,obj.Hy,X,Y,Z*ones(size(X)),'cubic')), ...
                    min(obj.OUT.data{n_obj}.X1z(:)),max(obj.OUT.data{n_obj}.X1z(:)), min(obj.OUT.data{n_obj}.Y1z(:)),max(obj.OUT.data{n_obj}.Y1z(:)))/Norm;
                Z = obj.OUT.data{n_obj}.Z2z(1);
                S_aver = S_aver - integral2(@(X,Y)(interpn(obj.xEy,obj.yEy,obj.zEy,obj.Ey,X,Y,Z*ones(size(X)),'cubic').*interpn(obj.xHz,obj.yHz,obj.zHz,obj.Hz,X,Y,Z*ones(size(X)),'cubic') ...
                                                  -interpn(obj.xEz,obj.yEz,obj.zEz,obj.Ez,X,Y,Z*ones(size(X)),'cubic').*interpn(obj.xHy,obj.yHy,obj.zHy,obj.Hy,X,Y,Z*ones(size(X)),'cubic')), ...
                    min(obj.OUT.data{n_obj}.X1z(:)),max(obj.OUT.data{n_obj}.X1z(:)), min(obj.OUT.data{n_obj}.Y1z(:)),max(obj.OUT.data{n_obj}.Y1z(:)))/Norm;

                obj.OUT.data{n_obj}.S_aver(time) = S_aver;
        end
    end
end
