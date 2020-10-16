%% update electric field
function obj = updateE(obj)
    global eps0
    switch obj.Dimensionality
        case 1
            if obj.isPML
                DxHy = diff(obj.Hy,[],1)/obj.dx;

                %Update of the PML-Matrices
                obj.PML.QxHy = obj.PML.bxEz.*(obj.PML.QxHy+DxHy)-DxHy;

                obj.Ez(2:end-1) = obj.Ez(2:end-1) ...
                    +(obj.dt/eps0./obj.epsilon_z).*(DxHy+obj.PML.QxHy); % + JEz(2:end,2:end)
            else
                obj.Ez(2:end-1) = obj.coefEz.*obj.Ez(2:end-1) + obj.coefEz_y.*diff(obj.Hy,1);
            end
        case 2
            if obj.isTE  % TM polarization
                if obj.isPML
                    DxHy = diff(obj.Hy(:,2:end-1),[],1)/obj.dx;
                    DyHx = diff(obj.Hx(2:end-1,:),[],2)/obj.dy;

                    %Update of the PML-Matrices
                    obj.PML.QxHy = obj.PML.bxEz.*(obj.PML.QxHy+DxHy)-DxHy;
                    obj.PML.QyHx = obj.PML.byEz.*(obj.PML.QyHx+DyHx)-DyHx;

                    obj.Ez(2:end-1,2:end-1) = obj.Ez(2:end-1,2:end-1) ...
                        +(obj.dt/eps0./obj.epsilon_z).*(DxHy+obj.PML.QxHy -DyHx-obj.PML.QyHx); % + JEz(2:end,2:end)
                else
                    obj.Ez(2:end-1,2:end-1) = obj.coefEz.*obj.Ez(2:end-1,2:end-1) ...
                        + obj.coefEz_y.*diff(obj.Hy(:,2:end-1),[],1) - obj.coefEz_x.*diff(obj.Hx(2:end-1,:),[],2);
                end
            else % TE polarization
                if obj.isPML
                    DyHz = diff(obj.Hz,[],2)/obj.dy;
                    DxHz = diff(obj.Hz,[],1)/obj.dx;

                    %Update of the PML-Matrices
                    obj.PML.QyHz = obj.PML.byEx.*(obj.PML.QyHz+DyHz)-DyHz;
                    obj.PML.QxHz = obj.PML.bxEy.*(obj.PML.QxHz+DxHz)-DxHz;

                    obj.Ex(:,2:end-1,2:end-1) = obj.Ex(:,2:end-1,2:end-1) ...
                        +(obj.dt/eps0./obj.epsilon_x).*( DyHz+obj.PML.QyHz); % + JEz(2:end,2:end)
                    obj.Ey(2:end-1,:,2:end-1) = obj.Ey(2:end-1,:,2:end-1) ...
                        +(obj.dt/eps0./obj.epsilon_y).*(-DxHz-obj.PML.QxHz); % + JEz(2:end,2:end)
                else
                    obj.Ex(:,2:end-1) = obj.coefEx.*obj.Ex(:,2:end-1) ...
                        + obj.coefEx_z.*diff(obj.Hz,[],2);
                    obj.Ey(2:end-1,:) = obj.coefEy.*obj.Ey(2:end-1,:) ...
                        - obj.coefEy_z.*diff(obj.Hz,[],1);
                end
            end
        case 3
            if obj.isPML
                DyHz = diff(obj.Hz(:,:,2:end-1),[],2)/obj.dy;
                DzHy = diff(obj.Hy(:,2:end-1,:),[],3)/obj.dz;
                DxHz = diff(obj.Hz(:,:,2:end-1),[],1)/obj.dx;
                DzHx = diff(obj.Hx(2:end-1,:,:),[],3)/obj.dz;
                DxHy = diff(obj.Hy(:,2:end-1,:),[],1)/obj.dx;
                DyHx = diff(obj.Hx(2:end-1,:,:),[],2)/obj.dy;

                %Update of the PML-Matrices
                obj.PML.QyHz = obj.PML.byEx.*(obj.PML.QyHz+DyHz)-DyHz;
                obj.PML.QzHy = obj.PML.bzEx.*(obj.PML.QzHy+DzHy)-DzHy;
                obj.PML.QxHz = obj.PML.bxEy.*(obj.PML.QxHz+DxHz)-DxHz;
                obj.PML.QzHx = obj.PML.bzEy.*(obj.PML.QzHx+DzHx)-DzHx;
                obj.PML.QxHy = obj.PML.bxEz.*(obj.PML.QxHy+DxHy)-DxHy;
                obj.PML.QyHx = obj.PML.byEz.*(obj.PML.QyHx+DyHx)-DyHx;

                obj.Ex(:,2:end-1,2:end-1) = obj.Ex(:,2:end-1,2:end-1) ...
                    +(obj.dt/eps0./obj.epsilon_x).*(DyHz+obj.PML.QyHz -DzHy-obj.PML.QzHy); % + JEz(2:end,2:end)
                obj.Ey(2:end-1,:,2:end-1) = obj.Ey(2:end-1,:,2:end-1) ...
                    +(obj.dt/eps0./obj.epsilon_y).*(DzHx+obj.PML.QzHx -DxHz-obj.PML.QxHz); % + JEz(2:end,2:end)
                obj.Ez(2:end-1,2:end-1,:) = obj.Ez(2:end-1,2:end-1,:) ...
                    +(obj.dt/eps0./obj.epsilon_z).*(DxHy+obj.PML.QxHy -DyHx-obj.PML.QyHx); % + JEz(2:end,2:end)
            else
                obj.Ex(:,2:end-1,2:end-1) = obj.coefEx.*obj.Ex(:,2:end-1,2:end-1) ...
                    + obj.coefEx_z.*diff(obj.Hz(:,:,2:end-1),[],2) - obj.coefEx_y.*diff(obj.Hy(:,2:end-1,:),[],3);
                obj.Ey(2:end-1,:,2:end-1) = obj.coefEy.*obj.Ey(2:end-1,:,2:end-1) ...
                    + obj.coefEy_x.*diff(obj.Hx(2:end-1,:,:),[],3) - obj.coefEy_z.*diff(obj.Hz(:,:,2:end-1),[],1);
                obj.Ez(2:end-1,2:end-1,:) = obj.coefEz.*obj.Ez(2:end-1,2:end-1,:) ...
                    + obj.coefEz_y.*diff(obj.Hy(:,2:end-1,:),[],1) - obj.coefEz_x.*diff(obj.Hx(2:end-1,:,:),[],2);
            end
    end
end
