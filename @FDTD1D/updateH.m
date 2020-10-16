%% update magnetic field
function obj = updateH(obj)
%#codegen

switch obj.PMLtype
	case 1
		DxEz = diff(obj.Ez,[],1);

		%Update of the PML-Matrices
		obj.PML.QxEz = obj.PML.bxHy.*(obj.PML.QxEz+DxEz)-DxEz;

		obj.Hy = obj.coefHy.*obj.Hy + obj.coefHy_z.*(DxEz+obj.PML.QxEz); % + JHy(2:end,1:end-1)
	case 2
        %% Calculate My -> Hy
        obj.Hy = obj.PML.cHy_1.*obj.Hy + obj.PML.cHy_2.*diff(obj.Ez,[],1)*obj.Courant(1);
	otherwise
		obj.Hy = obj.coefHy.*obj.Hy + obj.coefHy_z.*diff(obj.Ez,[],1);
end

end
