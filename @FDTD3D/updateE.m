%% update electric field
function obj = updateE(obj)
%#codegen

global imp0

switch obj.PMLtype
	case 1
		DyHz = diff(obj.Hz(:,:,2:end-1),[],2);
		DzHy = diff(obj.Hy(:,2:end-1,:),[],3);
		DxHz = diff(obj.Hz(:,:,2:end-1),[],1);
		DzHx = diff(obj.Hx(2:end-1,:,:),[],3);
		DxHy = diff(obj.Hy(:,2:end-1,:),[],1);
		DyHx = diff(obj.Hx(2:end-1,:,:),[],2);

		%Update of the PML-Matrices
		obj.PML.QyHz = obj.PML.byEx.*(obj.PML.QyHz+DyHz)-DyHz;
		obj.PML.QzHy = obj.PML.bzEx.*(obj.PML.QzHy+DzHy)-DzHy;
		obj.PML.QxHz = obj.PML.bxEy.*(obj.PML.QxHz+DxHz)-DxHz;
		obj.PML.QzHx = obj.PML.bzEy.*(obj.PML.QzHx+DzHx)-DzHx;
		obj.PML.QxHy = obj.PML.bxEz.*(obj.PML.QxHy+DxHy)-DxHy;
		obj.PML.QyHx = obj.PML.byEz.*(obj.PML.QyHx+DyHx)-DyHx;

		obj.Ex(:,2:end-1,2:end-1) = obj.coefEx.*obj.Ex(:,2:end-1,2:end-1) ...
			+ obj.coefEx_z.*(DyHz+obj.PML.QyHz) - obj.coefEx_y.*(DzHy+obj.PML.QzHy); % + JEz(2:end,2:end)
		obj.Ey(2:end-1,:,2:end-1) = obj.coefEy.*obj.Ey(2:end-1,:,2:end-1) ...
			+ obj.coefEy_x.*(DzHx+obj.PML.QzHx) - obj.coefEy_z.*(DxHz+obj.PML.QxHz); % + JEz(2:end,2:end)
		obj.Ez(2:end-1,2:end-1,:) = obj.coefEz.*obj.Ez(2:end-1,2:end-1,:) ...
			+ obj.coefEz_y.*(DxHy+obj.PML.QxHy) - obj.coefEz_x.*(DyHx+obj.PML.QyHx); % + JEz(2:end,2:end)
	case 2
        %% Calculate Fx -> Gx -> Ex(:,2:end-1,2:end-1)
        test = -1;
        
        Fx0 = obj.PML.Fx;
        Gx0 = obj.PML.Gx;
        switch test
            case -1
                REx0 = obj.CCPR.REx;
                obj.CCPR.REx = obj.CCPR.REx + ((obj.Courant(2)*diff(obj.Hz(:,:,2:end-1),[],2) - obj.Courant(3)*diff(obj.Hy(:,2:end-1,:),[],3))*imp0 ...
                                              - sum(real(obj.CCPR.gammaEx.*obj.CCPR.Jx),4))./obj.epsilon_xx;
                obj.CCPR.Jx = obj.CCPR.alphaEx.*obj.CCPR.Jx + bsxfun(@times,obj.CCPR.betaEx,obj.CCPR.REx-REx0);
                
                obj.PML.Gx = obj.PML.cGx_1.*Gx0 + obj.PML.cGx_2.*obj.CCPR.REx - obj.PML.cGx_3.*REx0;
                obj.Ex(:,2:end-1,2:end-1) = obj.PML.cEx_1.*obj.Ex(:,2:end-1,2:end-1) + obj.PML.cEx_2.*obj.PML.Gx - obj.PML.cEx_3.*Gx0;
            case 0
                obj.PML.Fx = obj.PML.cFx_1.*Fx0 ...
                           + obj.PML.cFx_2.*(obj.Courant(2)*diff(obj.Hz(:,:,2:end-1),[],2) - obj.Courant(3)*diff(obj.Hy(:,2:end-1,:),[],3));
                
                obj.PML.Gx = obj.PML.cGx_1.*Gx0 + obj.PML.cGx_2.*(obj.PML.Fx - Fx0);
                obj.Ex(:,2:end-1,2:end-1) = obj.PML.cEx_1.*obj.Ex(:,2:end-1,2:end-1) + obj.PML.cEx_2.*obj.PML.Gx - obj.PML.cEx_3.*Gx0;
            case 1
                curlEx = obj.PML.cFx_2.*(obj.Courant(2)*diff(obj.Hz(:,:,2:end-1),[],2) - obj.Courant(3)*diff(obj.Hy(:,2:end-1,:),[],3));
                obj.PML.Fx = obj.PML.cFx_1.*Fx0 + curlEx(obj.conditions.ExPML);
                obj.PML.Gx = obj.PML.cGx_1.*Gx0 + obj.PML.cGx_2.*(obj.PML.Fx - Fx0);
                cond_tmp = false(size(obj.Ex)); cond_tmp(:,2:end-1,2:end-1) = obj.conditions.ExPML;
                obj.Ex(cond_tmp) = obj.PML.cEx_1.*obj.Ex(cond_tmp) + obj.PML.cEx_2.*obj.PML.Gx - obj.PML.cEx_3.*Gx0;
                cond_tmp = false(size(obj.Ex)); cond_tmp(:,2:end-1,2:end-1) = ~obj.conditions.ExPML;
                obj.Ex(cond_tmp) = obj.Ex(cond_tmp) + curlEx(~obj.conditions.ExPML);
        end
        %% Calculate Fy -> Gy -> Ey(2:end-1,:,2:end-1)
        Fy0 = obj.PML.Fy;
        Gy0 = obj.PML.Gy;
        switch test
            case -1
                REy0 = obj.CCPR.REy;
                obj.CCPR.REy = obj.CCPR.REy + ((obj.Courant(3)*diff(obj.Hx(2:end-1,:,:),[],3) - obj.Courant(1)*diff(obj.Hz(:,:,2:end-1),[],1))*imp0 ...
                                              - sum(real(obj.CCPR.gammaEy.*obj.CCPR.Jy),4))./obj.epsilon_yy;
                obj.CCPR.Jy = obj.CCPR.alphaEy.*obj.CCPR.Jy + bsxfun(@times,obj.CCPR.betaEy,obj.CCPR.REy-REy0);
                
                obj.PML.Gy = obj.PML.cGy_1.*Gy0 + obj.PML.cGy_2.*obj.CCPR.REy - obj.PML.cGy_3.*REy0;
                obj.Ey(2:end-1,:,2:end-1) = obj.PML.cEy_1.*obj.Ey(2:end-1,:,2:end-1) + obj.PML.cEy_2.*obj.PML.Gy - obj.PML.cEy_3.*Gy0;
            case 0
                obj.PML.Fy = obj.PML.cFy_1.*Fy0 ...
                           + obj.PML.cFy_2.*(obj.Courant(3)*diff(obj.Hx(2:end-1,:,:),[],3) - obj.Courant(1)*diff(obj.Hz(:,:,2:end-1),[],1));
                
                obj.PML.Gy = obj.PML.cGy_1.*Gy0 + obj.PML.cGy_2.*(obj.PML.Fy - Fy0);
                obj.Ey(2:end-1,:,2:end-1) = obj.PML.cEy_1.*obj.Ey(2:end-1,:,2:end-1) + obj.PML.cEy_2.*obj.PML.Gy - obj.PML.cEy_3.*Gy0;
            case 1
                curlEy = obj.PML.cFy_2.*(obj.Courant(3)*diff(obj.Hx(2:end-1,:,:),[],3) - obj.Courant(1)*diff(obj.Hz(:,:,2:end-1),[],1));
                obj.PML.Fy = obj.PML.cFy_1.*Fy0 + curlEy(obj.conditions.EyPML);
                obj.PML.Gy = obj.PML.cGy_1.*Gy0 + obj.PML.cGy_2.*(obj.PML.Fy - Fy0);
                cond_tmp = false(size(obj.Ey)); cond_tmp(2:end-1,:,2:end-1) = obj.conditions.EyPML;
                obj.Ey(cond_tmp) = obj.PML.cEy_1.*obj.Ey(cond_tmp) + obj.PML.cEy_2.*obj.PML.Gy - obj.PML.cEy_3.*Gy0;
                cond_tmp = false(size(obj.Ey)); cond_tmp(2:end-1,:,2:end-1) = ~obj.conditions.EyPML;
                obj.Ey(cond_tmp) = obj.Ey(cond_tmp) + curlEy(~obj.conditions.EyPML);
        end
        %% Calculate Fz -> Gz -> Ez(2:end-1,2:end-1,:)
        Fz0 = obj.PML.Fz;
        Gz0 = obj.PML.Gz;
        switch test
            case -1
                REz0 = obj.CCPR.REz;
                obj.CCPR.REz = obj.CCPR.REz + ((obj.Courant(1)*diff(obj.Hy(:,2:end-1,:),[],1) - obj.Courant(2)*diff(obj.Hx(2:end-1,:,:),[],2))*imp0 ...
                                              - sum(real(obj.CCPR.gammaEz.*obj.CCPR.Jz),4))./obj.epsilon_zz;
                obj.CCPR.Jz = obj.CCPR.alphaEz.*obj.CCPR.Jz + bsxfun(@times, obj.CCPR.betaEz, obj.CCPR.REz-REz0);
                
                obj.PML.Gz = obj.PML.cGz_1.*Gz0 + obj.PML.cGz_2.*obj.CCPR.REz - obj.PML.cGz_3.*REz0;
                obj.Ez(2:end-1,2:end-1,:) = obj.PML.cEz_1.*obj.Ez(2:end-1,2:end-1,:) + obj.PML.cEz_2.*obj.PML.Gz - obj.PML.cEz_3.*Gz0;
            case 0
                obj.PML.Fz = obj.PML.cFz_1.*Fz0 ...
                           + obj.PML.cFz_2.*(obj.Courant(1)*diff(obj.Hy(:,2:end-1,:),[],1) - obj.Courant(2)*diff(obj.Hx(2:end-1,:,:),[],2));
                
                obj.PML.Gz = obj.PML.cGz_1.*Gz0 + obj.PML.cGz_2.*(obj.PML.Fz - Fz0);
                obj.Ez(2:end-1,2:end-1,:) = obj.PML.cEz_1.*obj.Ez(2:end-1,2:end-1,:) + obj.PML.cEz_2.*obj.PML.Gz - obj.PML.cEz_3.*Gz0;
            case 1
                curlEz = obj.PML.cFz_2.*(obj.Courant(1)*diff(obj.Hy(:,2:end-1,:),[],1) - obj.Courant(2)*diff(obj.Hx(2:end-1,:,:),[],2));
                obj.PML.Fz = obj.PML.cFz_1.*Fz0 + curlEz(obj.conditions.EzPML);
                obj.PML.Gz = obj.PML.cGz_1.*Gz0 + obj.PML.cGz_2.*(obj.PML.Fz - Fz0);
                cond_tmp = false(size(obj.Ez)); cond_tmp(2:end-1,2:end-1,:) = obj.conditions.EzPML;
                obj.Ez(cond_tmp) = obj.PML.cEz_1.*obj.Ez(cond_tmp) + obj.PML.cEz_2.*obj.PML.Gz - obj.PML.cEz_3.*Gz0;
                cond_tmp = false(size(obj.Ez)); cond_tmp(2:end-1,2:end-1,:) = ~obj.conditions.EzPML;
                obj.Ez(cond_tmp) = obj.Ez(cond_tmp) + curlEz(~obj.conditions.EzPML);
        end
    otherwise
		obj.Ex(:,2:end-1,2:end-1) = obj.coefEx.*obj.Ex(:,2:end-1,2:end-1) ...
			+ obj.coefEx_z.*diff(obj.Hz(:,:,2:end-1),[],2) - obj.coefEx_y.*diff(obj.Hy(:,2:end-1,:),[],3);
		obj.Ey(2:end-1,:,2:end-1) = obj.coefEy.*obj.Ey(2:end-1,:,2:end-1) ...
			+ obj.coefEy_x.*diff(obj.Hx(2:end-1,:,:),[],3) - obj.coefEy_z.*diff(obj.Hz(:,:,2:end-1),[],1);
		obj.Ez(2:end-1,2:end-1,:) = obj.coefEz.*obj.Ez(2:end-1,2:end-1,:) ...
			+ obj.coefEz_y.*diff(obj.Hy(:,2:end-1,:),[],1) - obj.coefEz_x.*diff(obj.Hx(2:end-1,:,:),[],2);
end

end
