%% get grid coordinates
function [X,Y,Z] = getGridXYZ(obj,type)
%#codegen

if nargin==1
	type = '';
end

switch type
	case 'Ex' % Ex grids
		X = linspace(obj.FDTDspace(1,1)+obj.dx/2,obj.FDTDspace(1,2)-obj.dx/2,obj.Nx);
		Y = linspace(obj.FDTDspace(2,1),obj.FDTDspace(2,2),obj.Ny+1);
		Z = linspace(obj.FDTDspace(3,1),obj.FDTDspace(3,2),obj.Nz+1);
	case 'Ey'
		X = linspace(obj.FDTDspace(1,1),obj.FDTDspace(1,2),obj.Nx+1);
		Y = linspace(obj.FDTDspace(2,1)+obj.dx/2,obj.FDTDspace(2,2)-obj.dx/2,obj.Ny);
		Z = linspace(obj.FDTDspace(3,1),obj.FDTDspace(3,2),obj.Nz+1);
	case 'Ez'
		X = linspace(obj.FDTDspace(1,1),obj.FDTDspace(1,2),obj.Nx+1);
		Y = linspace(obj.FDTDspace(2,1),obj.FDTDspace(2,2),obj.Ny+1);
		Z = linspace(obj.FDTDspace(3,1)+obj.dx/2,obj.FDTDspace(3,2)-obj.dx/2,obj.Nz);
	case 'Hx' % Hx grids
		X = linspace(obj.FDTDspace(1,1),obj.FDTDspace(1,2),obj.Nx+1);
		Y = linspace(obj.FDTDspace(2,1)+obj.dx/2,obj.FDTDspace(2,2)-obj.dx/2,obj.Ny);
		Z = linspace(obj.FDTDspace(3,1)+obj.dx/2,obj.FDTDspace(3,2)-obj.dx/2,obj.Nz);
	case 'Hy'
		X = linspace(obj.FDTDspace(1,1)+obj.dx/2,obj.FDTDspace(1,2)-obj.dx/2,obj.Nx);
		Y = linspace(obj.FDTDspace(2,1),obj.FDTDspace(2,2),obj.Ny+1);
		Z = linspace(obj.FDTDspace(3,1)+obj.dx/2,obj.FDTDspace(3,2)-obj.dx/2,obj.Nz);
	case 'Hz'
		X = linspace(obj.FDTDspace(1,1)+obj.dx/2,obj.FDTDspace(1,2)-obj.dx/2,obj.Nx);
		Y = linspace(obj.FDTDspace(2,1)+obj.dx/2,obj.FDTDspace(2,2)-obj.dx/2,obj.Ny);
		Z = linspace(obj.FDTDspace(3,1),obj.FDTDspace(3,2),obj.Nz+1);
	otherwise
		% averaged grid
		X = linspace(obj.FDTDspace(1,1)+obj.dx/2,obj.FDTDspace(1,2)-obj.dx/2,obj.Nx);
		Y = linspace(obj.FDTDspace(2,1)+obj.dy/2,obj.FDTDspace(2,2)-obj.dy/2,obj.Ny);
		Z = linspace(obj.FDTDspace(3,1)+obj.dz/2,obj.FDTDspace(3,2)-obj.dx/2,obj.Nz);
end

end
