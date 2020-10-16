%% update magnetic field
function obj = updateH(obj)
%#codegen

if obj.isTE  % TM polarization
	switch obj.PMLtype
		case 1
			DyEz = diff(obj.Ez,[],2);
			DxEz = diff(obj.Ez,[],1);

			%Update of the PML-Matrices
			obj.PML.QyEz = obj.PML.byHx.*(obj.PML.QyEz+DyEz)-DyEz;
			obj.PML.QxEz = obj.PML.bxHy.*(obj.PML.QxEz+DxEz)-DxEz;

			obj.Hx = obj.coefHx.*obj.Hx - obj.coefHx_z.*(DyEz+obj.PML.QyEz); % + JHx(1:end-1,2:end) 
			obj.Hy = obj.coefHy.*obj.Hy + obj.coefHy_z.*(DxEz+obj.PML.QxEz); % + JHy(2:end,1:end-1)
		case 2
            %% Calculate Mx -> Hx
            Mx0 = obj.PML.Mx;
            obj.PML.Mx = obj.PML.cMx_1.*Mx0 - obj.PML.cMx_2.*diff(obj.Ez,[],2);
            obj.Hx     = obj.Hx + obj.PML.cHx_2.*obj.PML.Mx - obj.PML.cHx_3.*Mx0;

            %% Calculate My -> Hy
            My0 = obj.PML.My;
            obj.PML.My = My0 + obj.PML.cMy_2.*diff(obj.Ez,[],1);
            obj.Hy     = obj.PML.cHy_1.*obj.Hy + obj.PML.cHy_2.*obj.PML.My - obj.PML.cHy_3.*My0;
		otherwise
			obj.Hx = obj.coefHx.*obj.Hx - obj.coefHx_z.*diff(obj.Ez,[],2);
			obj.Hy = obj.coefHy.*obj.Hy + obj.coefHy_z.*diff(obj.Ez,[],1);
	end
else % TE polarization
	switch obj.PMLtype
		case 1
			DyEx = diff(obj.Ex,[],2);
			DxEy = diff(obj.Ey,[],1);

			%Update of the PML-Matrices
			obj.PML.QyEx = obj.PML.byHz.*(obj.PML.QyEx+DyEx)-DyEx;
			obj.PML.QxEy = obj.PML.bxHz.*(obj.PML.QxEy+DxEy)-DxEy;

			obj.Hz = obj.coefHz.*obj.Hz + obj.coefHz_x.*(DyEx+obj.PML.QyEx) - obj.coefHz_y.*(DxEy+obj.PML.QxEy); % + JHz(2:end,1:end-1)
		case 2
            %% Calculate Mz -> Hz
            Mz0 = obj.PML.Mz;
            obj.PML.Mz = obj.PML.cMz_1.*Mz0 + obj.PML.cMz_2.*(obj.Courant(2)*diff(obj.Ex,[],2) - obj.Courant(1)*diff(obj.Ey,[],1));
            obj.Hz	   = obj.PML.cHz_1.*obj.Hz + obj.PML.cHz_2.*(obj.PML.Mz - Mz0);
		otherwise
			obj.Hz = obj.coefHz.*obj.Hz + obj.coefHz_x.*diff(obj.Ex,[],2) - obj.coefHz_y.*diff(obj.Ey,[],1);
	end
end

end
