%% update magnetic field
function obj = updateH(obj)
%#codegen

switch obj.PMLtype
	case 1
		DyEz = diff(obj.Ez,[],2);
		DzEy = diff(obj.Ey,[],3);
		DxEz = diff(obj.Ez,[],1);
		DzEx = diff(obj.Ex,[],3);
		DxEy = diff(obj.Ey,[],1);
		DyEx = diff(obj.Ex,[],2);

		%Update of the PML-Matrices
		obj.PML.QyEz = obj.PML.byHx.*(obj.PML.QyEz+DyEz)-DyEz;
		obj.PML.QzEy = obj.PML.bzHx.*(obj.PML.QzEy+DzEy)-DzEy;
		obj.PML.QxEz = obj.PML.bxHy.*(obj.PML.QxEz+DxEz)-DxEz;
		obj.PML.QzEx = obj.PML.bzHy.*(obj.PML.QzEx+DzEx)-DzEx;
		obj.PML.QxEy = obj.PML.bxHz.*(obj.PML.QxEy+DxEy)-DxEy;
		obj.PML.QyEx = obj.PML.byHz.*(obj.PML.QyEx+DyEx)-DyEx;

		obj.Hx = obj.coefHx.*obj.Hx + obj.coefHx_y.*(DzEy+obj.PML.QzEy) - obj.coefHx_z.*(DyEz+obj.PML.QyEz); % + JHx(1:end-1,2:end) 
		obj.Hy = obj.coefHy.*obj.Hy + obj.coefHy_z.*(DxEz+obj.PML.QxEz) - obj.coefHy_x.*(DzEx+obj.PML.QzEx); % + JHy(2:end,1:end-1)
		obj.Hz = obj.coefHz.*obj.Hz + obj.coefHz_x.*(DyEx+obj.PML.QyEx) - obj.coefHz_y.*(DxEy+obj.PML.QxEy); % + JHz(2:end,1:end-1)
	case 2
        test = true;
        %% Calculate Mx -> Hx
        Mx0 = obj.PML.Mx;
        if test
            obj.PML.Mx = obj.PML.cMx_1.*Mx0 + obj.PML.cMx_2.*(obj.Courant(3)*diff(obj.Ey,[],3) - obj.Courant(2)*diff(obj.Ez,[],2));
            obj.Hx     = obj.PML.cHx_1.*obj.Hx + obj.PML.cHx_2.*obj.PML.Mx - obj.PML.cHx_3.*Mx0;
        else
            curlMx = obj.PML.cMx_2.*(obj.Courant(3)*diff(obj.Ey,[],3) - obj.Courant(2)*diff(obj.Ez,[],2));
            obj.PML.Mx = obj.PML.cMx_1.*Mx0 + curlMx(obj.conditions.HxPML);
            obj.Hx( obj.conditions.HxPML) = obj.PML.cHx_1.*obj.Hx(obj.conditions.HxPML) + obj.PML.cHx_2.*obj.PML.Mx - obj.PML.cHx_3.*Mx0;
            obj.Hx(~obj.conditions.HxPML) = obj.Hx(~obj.conditions.HxPML) + curlMx(~obj.conditions.HxPML);
        end 
        %% Calculate My -> Hy
        My0 = obj.PML.My;
        if test
            obj.PML.My = obj.PML.cMy_1.*My0 + obj.PML.cMy_2.*(obj.Courant(1)*diff(obj.Ez,[],1) - obj.Courant(3)*diff(obj.Ex,[],3));
            obj.Hy     = obj.PML.cHy_1.*obj.Hy + obj.PML.cHy_2.*obj.PML.My - obj.PML.cHy_3.*My0;
        else
            curlMy = obj.PML.cMy_2.*(obj.Courant(1)*diff(obj.Ez,[],1) - obj.Courant(3)*diff(obj.Ex,[],3));
            obj.PML.My = obj.PML.cMy_1.*My0 + curlMy(obj.conditions.HyPML);
            obj.Hy( obj.conditions.HyPML)     = obj.PML.cHy_1.*obj.Hy( obj.conditions.HyPML) + obj.PML.cHy_2.*obj.PML.My - obj.PML.cHy_3.*My0;
            obj.Hy(~obj.conditions.HyPML)     = obj.Hy(~obj.conditions.HxPML) + curlMy(~obj.conditions.HyPML);
        end
        %% Calculate Mz -> Hz
        Mz0 = obj.PML.Mz;
        if test
            obj.PML.Mz = obj.PML.cMz_1.*Mz0 + obj.PML.cMz_2.*(obj.Courant(2)*diff(obj.Ex,[],2) - obj.Courant(1)*diff(obj.Ey,[],1));
            obj.Hz	   = obj.PML.cHz_1.*obj.Hz + obj.PML.cHz_2.*obj.PML.Mz - obj.PML.cHz_3.*Mz0;
        else
            curlMz = obj.PML.cMz_2.*(obj.Courant(2)*diff(obj.Ex,[],2) - obj.Courant(1)*diff(obj.Ey,[],1));
            obj.PML.Mz = obj.PML.cMz_1.*Mz0 + curlMz( obj.conditions.HzPML);
            obj.Hz( obj.conditions.HzPML)	   = obj.PML.cHz_1.*obj.Hz( obj.conditions.HzPML) + obj.PML.cHz_2.*obj.PML.Mz - obj.PML.cHz_3.*Mz0;
            obj.Hz(~obj.conditions.HzPML)	   = obj.Hz(~obj.conditions.HzPML) + curlMz(~obj.conditions.HzPML);
        end
    otherwise
		obj.Hx = obj.coefHx.*obj.Hx + obj.coefHx_y.*diff(obj.Ey,[],3) - obj.coefHx_z.*diff(obj.Ez,[],2);
		obj.Hy = obj.coefHy.*obj.Hy + obj.coefHy_z.*diff(obj.Ez,[],1) - obj.coefHy_x.*diff(obj.Ex,[],3);
		obj.Hz = obj.coefHz.*obj.Hz + obj.coefHz_x.*diff(obj.Ex,[],2) - obj.coefHz_y.*diff(obj.Ey,[],1);
end

end
