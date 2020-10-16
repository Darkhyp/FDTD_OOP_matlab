function out = snapshotE(obj,type)
    if ~exist('type','var')
        type = '';
    end
    if strcmpi(type,'intensity')
        switch obj.Dimensionality
            case 1
                out = obj.Ez.^2;
            case 2
                if obj.isTE
                    out = obj.Ez.^2;
                else
                    out = obj.snapshotE('x0').^2 + obj.snapshotE('y0').^2;
                end
            case 3
                out = obj.snapshotE('x0').^2 + obj.snapshotE('y0').^2 + obj.snapshotE('z0').^2;
        end
        return
    end

    switch obj.Dimensionality
        case 1
            out = obj.Ez;
        case 2
            if obj.isTE
                out = obj.Ez;
            else
                switch type
                    % return as it
                    case 'x'
                        out = obj.Ex;
                    case 'y'
                        out = obj.Ey;
                    % fit to space grid
                    case 'x0'
                        out = (obj.Ex(:,1:end-1)+obj.Ex(:,2:end))/2;
                    case 'y0'
                        out = (obj.Ey(1:end-1,:)+obj.Ey(2:end,:))/2;
                    otherwise
                        out = obj.Ex;
                end
            end
        case 3
            switch type
                % return as it
                case 'x'
                    out = obj.Ex;
                case 'y'
                    out = obj.Ey;
                case 'z'
                    out = obj.Ez;
                % fit to space grid
                case 'x0'
                    out = (obj.Ex(:,1:end-1,1:end-1)+obj.Ex(:,1:end-1,2:end)+obj.Ex(:,2:end,1:end-1)+obj.Ex(:,2:end,2:end))/4;
                case 'y0'
                    out = (obj.Ey(1:end-1,:,1:end-1)+obj.Ey(1:end-1,:,2:end)+obj.Ey(2:end,:,1:end-1)+obj.Ey(2:end,:,2:end))/4;
                case 'z0'
                    out = (obj.Ez(1:end-1,1:end-1,:)+obj.Ez(1:end-1,2:end,:)+obj.Ez(2:end,1:end-1,:)+obj.Ez(2:end,2:end,:))/4;
                otherwise
                    out = obj.Ez;
            end
    end
end
