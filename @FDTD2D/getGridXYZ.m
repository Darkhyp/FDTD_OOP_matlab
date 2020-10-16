%% get grid coordinates
function [X,Y] = getGridXYZ(obj,type)
%#codegen

if ~exist('type','var')
    type = '';
end

switch type
    case 'Ex' % Ex grids
        X = linspace(obj.FDTDspace(1,1)+obj.dx/2,obj.FDTDspace(1,2)-obj.dx/2,obj.Nx);
        Y = linspace(obj.FDTDspace(2,1),obj.FDTDspace(2,2),obj.Ny+1);
    case 'Ey'
        X = linspace(obj.FDTDspace(1,1),obj.FDTDspace(1,2),obj.Nx+1);
        Y = linspace(obj.FDTDspace(2,1)+obj.dx/2,obj.FDTDspace(2,2)-obj.dx/2,obj.Ny);
    case 'Ez'
        X = linspace(obj.FDTDspace(1,1),obj.FDTDspace(1,2),obj.Nx+1);
        Y = linspace(obj.FDTDspace(2,1),obj.FDTDspace(2,2),obj.Ny+1);
    case 'Hx' % Hx grids
        X = linspace(obj.FDTDspace(1,1),obj.FDTDspace(1,2),obj.Nx+1);
        Y = linspace(obj.FDTDspace(2,1)+obj.dx/2,obj.FDTDspace(2,2)-obj.dx/2,obj.Ny);
    case 'Hy'
        X = linspace(obj.FDTDspace(1,1)+obj.dx/2,obj.FDTDspace(1,2)-obj.dx/2,obj.Nx);
        Y = linspace(obj.FDTDspace(2,1),obj.FDTDspace(2,2),obj.Ny+1);
    case 'Hz'
        X = linspace(obj.FDTDspace(1,1)+obj.dx/2,obj.FDTDspace(1,2)-obj.dx/2,obj.Nx);
        Y = linspace(obj.FDTDspace(2,1)+obj.dx/2,obj.FDTDspace(2,2)-obj.dx/2,obj.Ny);
    otherwise
        % averaged grid
        X = linspace(obj.FDTDspace(1,1)+obj.dx/2,obj.FDTDspace(1,2)-obj.dx/2,obj.Nx);
        Y = linspace(obj.FDTDspace(2,1)+obj.dy/2,obj.FDTDspace(2,2)-obj.dy/2,obj.Ny);
end

end
