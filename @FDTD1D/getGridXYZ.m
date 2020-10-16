%% get grid coordinates
function X = getGridXYZ(obj,type)
%#codegen

if ~exist('type','var')
    type = '';
end

switch type
    case 'Ez'
        X = linspace(obj.FDTDspace(1,1),obj.FDTDspace(1,2),obj.Nx+1)';
    case 'Hy'
        X = linspace(obj.FDTDspace(1,1)+obj.dx/2,obj.FDTDspace(1,2)-obj.dx/2,obj.Nx)';
    otherwise
        % averaged grid
        X = linspace(obj.FDTDspace(1,1)+obj.dx/2,obj.FDTDspace(1,2)-obj.dx/2,obj.Nx)';
end

end
