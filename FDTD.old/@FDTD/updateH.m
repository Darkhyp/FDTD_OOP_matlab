%% update magnetic field
function obj = updateH(obj)
    global mu0
    switch obj.Dimensionality
        case 1
            if obj.isPML
                DxEz = diff(obj.Ez,[],1)/obj.dx;

                %Update of the PML-Matrices
                obj.PML.QxEz = obj.PML.bxHy.*(obj.PML.QxEz+DxEz)-DxEz;

                obj.Hy(:) = obj.Hy ...
                    + (obj.dt/mu0./obj.mu_y).*(DxEz+obj.PML.QxEz); % + JHy(2:end,1:end-1)
            else
                obj.Hy(:) = obj.coefHy.*obj.Hy + obj.coefHy_z.*diff(obj.Ez,[],1);
            end
        case 2
            if obj.isTE  % TM polarization
                if obj.isPML
                    DyEz = diff(obj.Ez,[],2)/obj.dy;
                    DxEz = diff(obj.Ez,[],1)/obj.dx;

                    %Update of the PML-Matrices
                    obj.PML.QyEz = obj.PML.byHx.*(obj.PML.QyEz+DyEz)-DyEz;
                    obj.PML.QxEz = obj.PML.bxHy.*(obj.PML.QxEz+DxEz)-DxEz;

                    obj.Hx(:,:) = obj.Hx ...
                        + (obj.dt/mu0./obj.mu_x).*(-DyEz-obj.PML.QyEz); % + JHx(1:end-1,2:end) 
                    obj.Hy(:,:) = obj.Hy ...
                        + (obj.dt/mu0./obj.mu_y).*( DxEz+obj.PML.QxEz); % + JHy(2:end,1:end-1)
                else
                    obj.Hx(:,:) = obj.coefHx.*obj.Hx - obj.coefHx_z.*diff(obj.Ez,[],2);
                    obj.Hy(:,:) = obj.coefHy.*obj.Hy + obj.coefHy_z.*diff(obj.Ez,[],1);
                end
            else % TE polarization
                if obj.isPML
                    DxEy = diff(obj.Ey,[],1)/obj.dx;
                    DyEx = diff(obj.Ex,[],2)/obj.dy;

                    %Update of the PML-Matrices
                    obj.PML.QxEy = obj.PML.bxHz.*(obj.PML.QxEy+DxEy)-DxEy;
                    obj.PML.QyEx = obj.PML.byHz.*(obj.PML.QyEx+DyEx)-DyEx;

                    obj.Hz(:,:) = obj.Hz ...
                        + (obj.dt/mu0./obj.mu_z).*(DyEx+obj.PML.QyEx -DxEy-obj.PML.QxEy); % + JHz(2:end,1:end-1)
                else
                    obj.Hz(:,:) = obj.coefHz.*obj.Hz + obj.coefHz_x.*diff(obj.Ex,[],2) - obj.coefHz_y.*diff(obj.Ey,[],1);
                end
            end
        case 3
            if obj.isPML
                DyEz = diff(obj.Ez,[],2)/obj.dy;
                DzEy = diff(obj.Ey,[],3)/obj.dz;
                DxEz = diff(obj.Ez,[],1)/obj.dx;
                DzEx = diff(obj.Ex,[],3)/obj.dz;
                DxEy = diff(obj.Ey,[],1)/obj.dx;
                DyEx = diff(obj.Ex,[],2)/obj.dy;

                %Update of the PML-Matrices
                obj.PML.QyEz = obj.PML.byHx.*(obj.PML.QyEz+DyEz)-DyEz;
                obj.PML.QzEy = obj.PML.bzHx.*(obj.PML.QzEy+DzEy)-DzEy;
                obj.PML.QxEz = obj.PML.bxHy.*(obj.PML.QxEz+DxEz)-DxEz;
                obj.PML.QzEx = obj.PML.bzHy.*(obj.PML.QzEx+DzEx)-DzEx;
                obj.PML.QxEy = obj.PML.bxHz.*(obj.PML.QxEy+DxEy)-DxEy;
                obj.PML.QyEx = obj.PML.byHz.*(obj.PML.QyEx+DyEx)-DyEx;

                obj.Hx(:,:,:) = obj.Hx ...
                    + (obj.dt/mu0./obj.mu_x).*(DzEy+obj.PML.QzEy -DyEz-obj.PML.QyEz); % + JHx(1:end-1,2:end) 
                obj.Hy(:,:,:) = obj.Hy ...
                    + (obj.dt/mu0./obj.mu_y).*(DxEz+obj.PML.QxEz -DzEx-obj.PML.QzEx); % + JHy(2:end,1:end-1)
                obj.Hz(:,:,:) = obj.Hz ...
                    + (obj.dt/mu0./obj.mu_z).*(DyEx+obj.PML.QyEx -DxEy-obj.PML.QxEy); % + JHz(2:end,1:end-1)
            else
                obj.Hx(:,:,:) = obj.coefHx.*obj.Hx + obj.coefHx_y.*diff(obj.Ey,[],3) - obj.coefHx_z.*diff(obj.Ez,[],2);
                obj.Hy(:,:,:) = obj.coefHy.*obj.Hy + obj.coefHy_z.*diff(obj.Ez,[],1) - obj.coefHy_x.*diff(obj.Ex,[],3);
                obj.Hz(:,:,:) = obj.coefHz.*obj.Hz + obj.coefHz_x.*diff(obj.Ex,[],2) - obj.coefHz_y.*diff(obj.Ey,[],1);
            end
    end
end
