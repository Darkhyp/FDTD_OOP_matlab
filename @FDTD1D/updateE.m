%% update electric field
function obj = updateE(obj)
%#codegen

switch obj.PMLtype
	case 1
		DxHy = diff(obj.Hy,[],1);

		%Update of the PML-Matrices
		obj.PML.QxHy = obj.PML.bxEz.*(obj.PML.QxHy+DxHy)-DxHy;

		obj.Ez(2:end-1) = obj.coefEz.*obj.Ez(2:end-1) + obj.coefEz_y.*(DxHy+obj.PML.QxHy); % + JEz(2:end,2:end)
	case 2
        %% Calculate Fz -> Gz -> Ez(2:end-1,2:end-1,:)
        Fz0 = obj.PML.Fz;
        Gz0 = obj.PML.Gz;
        obj.PML.Fz = obj.PML.cFz_1.*Fz0 ...
               	   + obj.PML.cFz_2.*diff(obj.Hy,[],1)*obj.Courant(1);
        obj.PML.Gz = obj.PML.cGz_1.*Gz0 + obj.PML.cGz_2.*(obj.PML.Fz - Fz0);
        obj.Ez(2:end-1) = obj.Ez(2:end-1) + obj.PML.Gz - Gz0;
	otherwise
		obj.Ez(2:end-1) = obj.coefEz.*obj.Ez(2:end-1) + obj.coefEz_y.*diff(obj.Hy,1);
end

end
