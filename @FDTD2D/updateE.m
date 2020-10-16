%% update electric field
function obj = updateE(obj)
%#codegen

if obj.isTE  % TM polarization
	switch obj.PMLtype
		case 1
			DxHy = diff(obj.Hy(:,2:end-1),[],1);
			DyHx = diff(obj.Hx(2:end-1,:),[],2);

			%Update of the PML-Matrices
			obj.PML.QxHy = obj.PML.bxEz.*(obj.PML.QxHy+DxHy)-DxHy;
			obj.PML.QyHx = obj.PML.byEz.*(obj.PML.QyHx+DyHx)-DyHx;

			obj.Ez(2:end-1,2:end-1) = obj.coefEz.*obj.Ez(2:end-1,2:end-1) ...
				+ obj.coefEz_y.*(DxHy+obj.PML.QxHy) - obj.coefEz_x.*(DyHx+obj.PML.QyHx); % + JEz(2:end,2:end)
		case 2
            %% Calculate Fz -> Gz -> Ez(2:end-1,2:end-1,:)
            Fz0 = obj.PML.Fz;
            Gz0 = obj.PML.Gz;
            obj.PML.Fz = obj.PML.cFz_1.*Fz0 ...
                       + obj.PML.cFz_2.*(obj.Courant(1)*diff(obj.Hy(:,2:end-1),[],1) - obj.Courant(2)*diff(obj.Hx(2:end-1,:),[],2));
            obj.PML.Gz = obj.PML.cGz_1.*Gz0 + obj.PML.cGz_2.*(obj.PML.Fz - Fz0);
            obj.Ez(2:end-1,2:end-1) = obj.PML.cEz_1.*obj.Ez(2:end-1,2:end-1) + obj.PML.cEz_2.*(obj.PML.Gz - Gz0);
		otherwise
			obj.Ez(2:end-1,2:end-1) = obj.coefEz.*obj.Ez(2:end-1,2:end-1) ...
				+ obj.coefEz_y.*diff(obj.Hy(:,2:end-1),[],1) - obj.coefEz_x.*diff(obj.Hx(2:end-1,:),[],2);
	end
else % TE polarization
	switch obj.PMLtype
		case 1
			DyHz = diff(obj.Hz,[],2);
			DxHz = diff(obj.Hz,[],1);

			%Update of the PML-Matrices
			obj.PML.QyHz = obj.PML.byEx.*(obj.PML.QyHz+DyHz)-DyHz;
			obj.PML.QxHz = obj.PML.bxEy.*(obj.PML.QxHz+DxHz)-DxHz;

			obj.Ex(:,2:end-1) = obj.coefEx.*obj.Ex(:,2:end-1) + obj.coefEx_z.*( DyHz+obj.PML.QyHz); % + JEz(2:end,2:end)
			obj.Ey(2:end-1,:) = obj.coefEy.*obj.Ey(2:end-1,:) + obj.coefEy_z.*(-DxHz-obj.PML.QxHz); % + JEz(2:end,2:end)
		case 2
            %% Calculate Fx -> Gx -> Ex(:,2:end-1,2:end-1)
            Fx0 = obj.PML.Fx;
            Gx0 = obj.PML.Gx;
            obj.PML.Fx = obj.PML.cFx_1.*Fx0 + obj.PML.cFx_2.*diff(obj.Hz,[],2);
            obj.PML.Gx = obj.PML.cGx_1.*Gx0 + obj.PML.cGx_2.*(obj.PML.Fx - Fx0);
            obj.Ex(:,2:end-1) = obj.Ex(:,2:end-1) + obj.PML.cEx_2.*obj.PML.Gx - obj.PML.cEx_3.*Gx0;

            %% Calculate Fy -> Gy -> Ey(2:end-1,:,2:end-1)
            Fy0 = obj.PML.Fy;
            Gy0 = obj.PML.Gy;
            obj.PML.Fy = obj.PML.cFy_1.*Fy0 - obj.PML.cFy_2.*diff(obj.Hz,[],1);
            obj.PML.Gy = Gy0 + (obj.PML.Fy - Fy0);
            obj.Ey(2:end-1,:) = obj.PML.cEy_1.*obj.Ey(2:end-1,:) + obj.PML.cEy_2.*obj.PML.Gy - obj.PML.cEy_3.*Gy0;
		otherwise
			obj.Ex(:,2:end-1) = obj.coefEx.*obj.Ex(:,2:end-1) + obj.coefEx_z.*diff(obj.Hz,[],2);
			obj.Ey(2:end-1,:) = obj.coefEy.*obj.Ey(2:end-1,:) - obj.coefEy_z.*diff(obj.Hz,[],1);
	end
end

end
