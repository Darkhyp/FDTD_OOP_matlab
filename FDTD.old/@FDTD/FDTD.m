classdef FDTD < handle
%#codegen

    properties (Access = public)
    
        % electro-magnetic field
        Ex,Ey,Ez, Hx,Hy,Hz

        % grid coordinates
        x, y, z
        xEx,yEx,zEx, xEy,yEy,zEy, xEz,yEz,zEz
        xHx,yHx,zHx, xHy,yHy,zHy, xHz,yHz,zHz
        
        % structure for out-field 
        OUT
        
    end
    properties (Access = private)
        
        isTE        = true;
        
        % Perfect-matched layers
        isPML       = false;
        PML
        
        StabilityFactor = 0.99;
        CourantNumber

        Dimensionality
        Nx, Ny, Nz
        FDTDspace
        dx = Inf;
        dy = Inf;
        dz = Inf;
        dt = 1;
        imp0
        
        coefHx,coefHx_y,coefHx_z, coefHy,coefHy_x,coefHy_z, coefHz,coefHz_x,coefHz_y % H-field coefficients
        coefEx,coefEx_y,coefEx_z, coefEy,coefEy_x,coefEy_z, coefEz,coefEz_x,coefEz_y % E-field coefficients
        
        mu_x, mu_y, mu_z
        epsilon_x, epsilon_y, epsilon_z

        % structure for total-field/scattered-field (TFSF) boundary
        TFSF
        
        % structure for scattered-field 
        SCAT

        
        % structure for absorbing boundary conditions
        ABC
        

    end
    
    methods
        
        function dt = get_dt(obj)
            dt = obj.dt;
        end
%         function set.CourantNumber(obj,CourantNumber)
%             obj.CourantNumber = CourantNumber;
%         end
        function CourantNumber = get_CourantNumber(obj)
            CourantNumber = obj.CourantNumber;
        end
        
        function obj = FDTD(DomainSize, Space, PML, isTE)
            global c
            
            obj.StabilityFactor = 0.99;
            
            obj.Dimensionality = length(DomainSize);
            obj.FDTDspace = Space;
            

            if ~isempty(PML)
                obj.isPML   = true;
                obj.PML     = PML;
            end
            
            switch obj.Dimensionality
                case 1 % 1D system
                    obj.Nx = DomainSize(1);

                    obj.dx = diff(Space)/obj.Nx;

                    if obj.isPML
                        obj.Nx = obj.Nx + 2*obj.PML.layers;

                        obj.PML.Wx      = obj.PML.layers*obj.dx;
                        obj.PML.x1      = Space(1,1);
                        obj.PML.x2      = Space(1,2);

                        obj.FDTDspace = Space+[-obj.PML.Wx,obj.PML.Wx];

                        % Ez grid (2:end-1,2:end-1,:)
                        obj.PML.QyHx	= zeros(obj.Nx-1,1);
                        obj.PML.QxHy	= zeros(obj.Nx-1,1);

                        % Hy grid
                        obj.PML.QxEz	= zeros(obj.Nx,1);
                        obj.PML.QzEx	= zeros(obj.Nx,1);
                    end                    

%                     obj.CourantNumber = 1; % Courant number
                    
                    obj.Hy          = zeros(obj.Nx, 1);
                    obj.coefHy      = ones( obj.Nx, 1);
                    obj.coefHy_z	= ones( obj.Nx, 1);
                    obj.mu_y        = ones( obj.Nx, 1);
                    
                    obj.Ez          = zeros(obj.Nx+1, 1);
                    obj.coefEz      = ones( obj.Nx-1, 1);
                    obj.coefEz_y	= ones( obj.Nx-1, 1);
                    obj.epsilon_z	= ones( obj.Nx-1, 1);
                    
%                     obj.Ey = @(obj)  obj.Ez;
%                     obj.Hz = @(obj) -obj.Hy;
                case 2 % 2D system
                    obj.Nx = DomainSize(1);
                    obj.Ny = DomainSize(2);

                    obj.dx = diff(Space(1,:))/obj.Nx;
                    obj.dy = diff(Space(2,:))/obj.Ny;

                    if obj.isPML
                        obj.Nx = obj.Nx + 2*obj.PML.layers;
                        obj.Ny = obj.Ny + 2*obj.PML.layers;
                        
                        obj.PML.Wx      = obj.PML.layers*obj.dx;
                        obj.PML.x1      = Space(1,1);
                        obj.PML.x2      = Space(1,2);

                        obj.PML.Wy      = obj.PML.layers*obj.dy;
                        obj.PML.y1      = Space(2,1);
                        obj.PML.y2      = Space(2,2);

                        obj.FDTDspace = Space+[-obj.PML.Wx,obj.PML.Wx; -obj.PML.Wy,obj.PML.Wy];
                    end                    

%                     obj.CourantNumber = 1/sqrt(2); % Courant number
                    
                    if exist('isTE','var') % check polarization
                        obj.isTE = isTE;
                    end
                    
                    if obj.isTE
                        obj.Hx          = zeros(obj.Nx+1, obj.Ny);
                        obj.coefHx      = ones( obj.Nx+1, obj.Ny);
                        obj.coefHx_z	= ones( obj.Nx+1, obj.Ny);
                        obj.mu_x        = ones( obj.Nx+1, obj.Ny);

                        obj.Hy          = zeros(obj.Nx, obj.Ny+1);
                        obj.coefHy      = ones( obj.Nx, obj.Ny+1);
                        obj.coefHy_z	= ones( obj.Nx, obj.Ny+1);
                        obj.mu_y        = ones( obj.Nx, obj.Ny+1);

                        obj.Ez          = zeros(obj.Nx+1, obj.Ny+1);
                        obj.coefEz      = ones( obj.Nx-1, obj.Ny-1);
                        obj.coefEz_x	= ones( obj.Nx-1, obj.Ny-1);
                        obj.coefEz_y	= ones( obj.Nx-1, obj.Ny-1);
                        obj.epsilon_z	= ones( obj.Nx-1, obj.Ny-1);

                        if obj.isPML
                            % Ez grid (2:end-1,2:end-1,:)
                            obj.PML.QyHx	= zeros(obj.Nx-1, obj.Ny-1);
                            obj.PML.QxHy	= zeros(obj.Nx-1, obj.Ny-1);

                            % Hx grid
                            obj.PML.QyEz	= zeros(obj.Nx+1, obj.Ny);
                            obj.PML.QzEy	= zeros(obj.Nx+1, obj.Ny);
                            % Hy grid
                            obj.PML.QxEz	= zeros(obj.Nx, obj.Ny+1);
                            obj.PML.QzEx	= zeros(obj.Nx, obj.Ny+1);
                        end                    
                    else
                        obj.Hz          = zeros(obj.Nx, obj.Ny);
                        obj.coefHz      = ones( obj.Nx, obj.Ny);
                        obj.coefHz_x	= ones( obj.Nx, obj.Ny);
                        obj.coefHz_y	= ones( obj.Nx, obj.Ny);
                        obj.mu_z        = ones( obj.Nx, obj.Ny);

                        obj.Ex          = zeros(obj.Nx, obj.Ny+1);
                        obj.coefEx      = ones( obj.Nx, obj.Ny-1);
                        obj.coefEx_z	= ones( obj.Nx, obj.Ny-1);
                        obj.epsilon_x	= ones( obj.Nx, obj.Ny-1);

                        obj.Ey          = zeros(obj.Nx+1, obj.Ny);
                        obj.coefEy      = ones( obj.Nx-1, obj.Ny);
                        obj.coefEy_z	= ones( obj.Nx-1, obj.Ny);
                        obj.epsilon_y	= ones( obj.Nx-1, obj.Ny);

                        if obj.isPML

                            % Ex grid (:,2:end-1)
                            obj.PML.QyHz	= zeros(obj.Nx, obj.Ny-1);
                            obj.PML.QzHy	= zeros(obj.Nx, obj.Ny-1);
                            % Ey grid (2:end-1,:,2:end-1)
                            obj.PML.QxHz	= zeros(obj.Nx-1, obj.Ny);
                            obj.PML.QzHx	= zeros(obj.Nx-1, obj.Ny);
                            
                            % Hz grid
                            obj.PML.QxEy	= zeros(obj.Nx, obj.Ny);
                            obj.PML.QyEx	= zeros(obj.Nx, obj.Ny);
                        end                    
                    end
                case 3 % 3D system
                    obj.Nx = DomainSize(1);
                    obj.Ny = DomainSize(2);
                    obj.Nz = DomainSize(3);

                    obj.dx = diff(Space(1,:))/obj.Nx;
                    obj.dy = diff(Space(2,:))/obj.Ny;
                    obj.dz = diff(Space(3,:))/obj.Nz;
                    
                    if obj.isPML
                        obj.Nx = obj.Nx + 2*obj.PML.layers;
                        obj.Ny = obj.Ny + 2*obj.PML.layers;
                        obj.Nz = obj.Nz + 2*obj.PML.layers;

                        obj.PML.Wx      = obj.PML.layers*obj.dx;
                        obj.PML.x1      = Space(1,1);
                        obj.PML.x2      = Space(1,2);
                        
                        obj.PML.Wy      = obj.PML.layers*obj.dy;
                        obj.PML.y1      = Space(2,1);
                        obj.PML.y2      = Space(2,2);
                        
                        obj.PML.Wz      = obj.PML.layers*obj.dz;
                        obj.PML.z1      = Space(3,1);
                        obj.PML.z2      = Space(3,2);

                        obj.FDTDspace = Space+[-obj.PML.Wx,obj.PML.Wx; -obj.PML.Wy,obj.PML.Wy; -obj.PML.Wz,obj.PML.Wz];

                        % Ex grid (:,2:end-1,2:end-1)
                        obj.PML.QyHz	= zeros(obj.Nx, obj.Ny-1, obj.Nz-1);
                        obj.PML.QzHy	= zeros(obj.Nx, obj.Ny-1, obj.Nz-1);
                        % Ey grid (2:end-1,:,2:end-1)
                        obj.PML.QxHz	= zeros(obj.Nx-1, obj.Ny, obj.Nz-1);
                        obj.PML.QzHx	= zeros(obj.Nx-1, obj.Ny, obj.Nz-1);
                        % Ez grid (2:end-1,2:end-1,:)
                        obj.PML.QyHx	= zeros(obj.Nx-1, obj.Ny-1, obj.Nz);
                        obj.PML.QxHy	= zeros(obj.Nx-1, obj.Ny-1, obj.Nz);

                        % Hx grid
                        obj.PML.QyEz	= zeros(obj.Nx+1, obj.Ny, obj.Nz);
                        obj.PML.QzEy	= zeros(obj.Nx+1, obj.Ny, obj.Nz);
                        % Hy grid
                        obj.PML.QxEz	= zeros(obj.Nx, obj.Ny+1, obj.Nz);
                        obj.PML.QzEx	= zeros(obj.Nx, obj.Ny+1, obj.Nz);
                        % Hz grid
                        obj.PML.QxEy	= zeros(obj.Nx, obj.Ny, obj.Nz+1);
                        obj.PML.QyEx	= zeros(obj.Nx, obj.Ny, obj.Nz+1);
                    end                    

%                     obj.CourantNumber = 1/sqrt(3); % Courant number

                    obj.Hx          = zeros(obj.Nx+1, obj.Ny, obj.Nz);
                    obj.coefHx      = ones( obj.Nx+1, obj.Ny, obj.Nz);
                    obj.coefHx_y	= ones( obj.Nx+1, obj.Ny, obj.Nz);
                    obj.coefHx_z	= ones( obj.Nx+1, obj.Ny, obj.Nz);
                    obj.mu_x        = ones( obj.Nx+1, obj.Ny, obj.Nz);

                    obj.Hy          = zeros(obj.Nx, obj.Ny+1, obj.Nz);
                    obj.coefHy      = ones( obj.Nx, obj.Ny+1, obj.Nz);
                    obj.coefHy_x	= ones( obj.Nx, obj.Ny+1, obj.Nz);
                    obj.coefHy_z	= ones( obj.Nx, obj.Ny+1, obj.Nz);
                    obj.mu_y        = ones( obj.Nx, obj.Ny+1, obj.Nz);

                    obj.Hz          = zeros(obj.Nx, obj.Ny, obj.Nz+1);
                    obj.coefHz      = ones( obj.Nx, obj.Ny, obj.Nz+1);
                    obj.coefHz_x	= ones( obj.Nx, obj.Ny, obj.Nz+1);
                    obj.coefHz_y	= ones( obj.Nx, obj.Ny, obj.Nz+1);
                    obj.mu_z        = ones( obj.Nx, obj.Ny, obj.Nz+1);

                    obj.Ex          = zeros(obj.Nx, obj.Ny+1, obj.Nz+1);
                    obj.coefEx      = ones( obj.Nx, obj.Ny-1, obj.Nz-1);
                    obj.coefEx_y	= ones( obj.Nx, obj.Ny-1, obj.Nz-1);
                    obj.coefEx_z	= ones( obj.Nx, obj.Ny-1, obj.Nz-1);
                    obj.epsilon_x	= ones( obj.Nx, obj.Ny-1, obj.Nz-1);

                    obj.Ey          = zeros(obj.Nx+1, obj.Ny, obj.Nz+1);
                    obj.coefEy      = ones( obj.Nx-1, obj.Ny, obj.Nz-1);
                    obj.coefEy_x	= ones( obj.Nx-1, obj.Ny, obj.Nz-1);
                    obj.coefEy_z	= ones( obj.Nx-1, obj.Ny, obj.Nz-1);
                    obj.epsilon_y	= ones( obj.Nx-1, obj.Ny, obj.Nz-1);

                    obj.Ez          = zeros(obj.Nx+1, obj.Ny+1, obj.Nz);
                    obj.coefEz      = ones( obj.Nx-1, obj.Ny-1, obj.Nz);
                    obj.coefEz_x	= ones( obj.Nx-1, obj.Ny-1, obj.Nz);
                    obj.coefEz_y	= ones( obj.Nx-1, obj.Ny-1, obj.Nz);
                    obj.epsilon_z	= ones( obj.Nx-1, obj.Ny-1, obj.Nz);
            end
            
            obj.dt = obj.StabilityFactor/(sqrt(1/obj.dx^2+1/obj.dy^2+1/obj.dz^2)*c);
            obj.CourantNumber = obj.dt*c./[obj.dx,obj.dy,obj.dz]; % Courant number

        end
        
    end
    
    methods (Access = public)

        % update magnetic field
        obj = initGrid(obj, Objects, CourantNumber)
        
        % update magnetic field
        obj = updateH(obj)

        % update electric field
        obj = updateE(obj)
        
        % Absorbing boundary conditioms
        obj = initABC(obj,type)
        obj = applyABC(obj)
        
        % TFSF
        obj = initTFSF(obj,TFSFregion,TFSFsource)
        obj = updateTFSF(obj,time)
        
        % out field data
        obj = initOUT(obj,Objects, MaxTime)
        obj = updateOUT(obj,time)
        
        % calculate source function at given time and location
        obj = add_Ez_inc(obj,Position,time,sourceStr)
        obj = add_Hy_inc(obj,Position,time,sourceStr)
        
        % field components outputs
        out = snapshotE(obj,type)
        out = snapshotH(obj,type)
        
        [Sx,Sy,Sz] = Poynting(obj,X,Y,Z)

        %% get grid coordinates
        [X,Y,Z] = getGridXYZ(obj,type)
%{        
        function out = getGridX(obj,type)
            if ~exist('type','var')
                type = '';
                out = (1/2:1:obj.Nx+1/2)*obj.dx;
            end
            switch type
                case 'Ex'
                    out = (0:obj.Nx-1)*obj.dx;
                case 'Ey'
                    out = (0:obj.Nx)*obj.dy;
                case 'Ez'
                    out = (0:obj.Nx)*obj.dz;
            end
        end
        function out = getGridY(obj)
            out = (1/2:1:obj.Ny+1/2)*obj.dy;
        end
        function out = getGridZ(obj)
            out = (1/2:1:obj.Nz+1/2)*obj.dz;
        end
%}        
    end

    methods(Access = protected)
        
        coef = coefABC(obj,coefEz_y,coefHy_z,order)
        
    end
end




