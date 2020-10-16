classdef FDTD < handle
% %#codegen

    properties (Access = public)
    
        Ex,Ey,Ez,Hx,Hy,Hz; % electro-magnetic field

        %grid coordinates
        x, y, z
        xEx, yEx, zEx, xEy, yEy, zEy, xEz, yEz, zEz
        xHx, yHx, zHx, xHy, yHy, zHy, xHz, yHz, zHz
        
        % structure for out-field 
        OUT
        
    end
    properties (Access = private)
        
        isTE        = true;
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
        
        m2, n2, p2

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
        
    end
    
    methods (Access = public)
        
        function obj = FDTD(DomainSize, Space, PML, isTE)
            global eps0 mu0 c imp0
%             eps0 = 8.8541878176209821e-12;
%             mu0  = 1.2566370614359173e-06;
%             c    = 299792458; % 1/sqrt(eps0*mu0)
            imp0 = sqrt(mu0/eps0);
%             imp0 = 377.0; 
            
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
                    obj.m2 = 2:obj.Nx;

%                     obj.CourantNumber = 1; % Courant number
                    
                    obj.Hy          = zeros(obj.Nx, 1);
                    obj.coefHy      = ones( obj.Nx, 1);
                    obj.coefHy_z	= ones( obj.Nx, 1);
                    obj.mu_y        = ones( obj.Nx, 1);
                    
                    obj.Ez          = zeros(obj.Nx+1, 1);
                    obj.coefEz      = ones( obj.Nx+1, 1);
                    obj.coefEz_y	= ones( obj.Nx+1, 1);
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
                    obj.m2 = 2:obj.Nx;
                    obj.n2 = 2:obj.Ny;

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
                        obj.coefEz      = ones( obj.Nx+1, obj.Ny+1);
                        obj.coefEz_x	= ones( obj.Nx+1, obj.Ny+1);
                        obj.coefEz_y	= ones( obj.Nx+1, obj.Ny+1);
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
                        obj.coefEx      = ones( obj.Nx, obj.Ny+1);
                        obj.coefEx_z	= ones( obj.Nx, obj.Ny+1);
                        obj.epsilon_x	= ones( obj.Nx, obj.Ny-1);

                        obj.Ey          = zeros(obj.Nx+1, obj.Ny);
                        obj.coefEy      = ones( obj.Nx+1, obj.Ny);
                        obj.coefEy_z	= ones( obj.Nx+1, obj.Ny);
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
                    obj.m2 = 2:obj.Nx; % 2:end-1
                    obj.n2 = 2:obj.Ny; % 2:end-1
                    obj.p2 = 2:obj.Nz; % 2:end-1

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
                    obj.coefEx      = ones( obj.Nx, obj.Ny+1, obj.Nz+1);
                    obj.coefEx_y	= ones( obj.Nx, obj.Ny+1, obj.Nz+1);
                    obj.coefEx_z	= ones( obj.Nx, obj.Ny+1, obj.Nz+1);
                    obj.epsilon_x	= ones( obj.Nx, obj.Ny-1, obj.Nz-1);

                    obj.Ey          = zeros(obj.Nx+1, obj.Ny, obj.Nz+1);
                    obj.coefEy      = ones( obj.Nx+1, obj.Ny, obj.Nz+1);
                    obj.coefEy_x	= ones( obj.Nx+1, obj.Ny, obj.Nz+1);
                    obj.coefEy_z	= ones( obj.Nx+1, obj.Ny, obj.Nz+1);
                    obj.epsilon_y	= ones( obj.Nx-1, obj.Ny, obj.Nz-1);

                    obj.Ez          = zeros(obj.Nx+1, obj.Ny+1, obj.Nz);
                    obj.coefEz      = ones( obj.Nx+1, obj.Ny+1, obj.Nz);
                    obj.coefEz_x	= ones( obj.Nx+1, obj.Ny+1, obj.Nz);
                    obj.coefEz_y	= ones( obj.Nx+1, obj.Ny+1, obj.Nz);
                    obj.epsilon_z	= ones( obj.Nx-1, obj.Ny-1, obj.Nz);
            end
            
            obj.dt = obj.StabilityFactor/(sqrt(1/obj.dx^2+1/obj.dy^2+1/obj.dz^2)*c);
            obj.CourantNumber = obj.dt*c./[obj.dx,obj.dy,obj.dz]; % Courant number

        end
        
        function obj = initGrid(obj, Objects, CourantNumber)
            global eps0 imp0 c
            if ~isempty(CourantNumber)
                obj.CourantNumber = CourantNumber;
            end
            
            switch obj.Dimensionality
                case 1 % 1D system
                    if isempty(Objects)
                        % set electric-field update coefficients
                        obj.coefEz_y(:) = obj.CourantNumber(1) * imp0;
                        % set magnetic-field update coefficients */
                        obj.coefHy_z(:) = obj.CourantNumber(1) / imp0;

                        %grid coordinates
                        obj.xEz = obj.getGridXYZ('Ez');
                        obj.xHy = obj.getGridXYZ('Hy');
                    end
                case 2 % 2D system
                    if obj.isTE
                        % set electric-field update coefficients
                        obj.coefEz_y(:,:) = obj.CourantNumber(1) * imp0;
                        obj.coefEz_x(:,:) = obj.CourantNumber(2) * imp0;
                        % set magnetic-field update coefficients */
                        obj.coefHy_z(:,:) = obj.CourantNumber(1) / imp0;
                        obj.coefHx_z(:,:) = obj.CourantNumber(2) / imp0;
                        
                        %grid coordinates
                        [obj.xEz, obj.yEz] = obj.getGridXYZ('Ez');
                        [obj.xHx, obj.yHx] = obj.getGridXYZ('Hx');
                        [obj.xHy, obj.yHy] = obj.getGridXYZ('Hy');
 
                    else
                        % set electric-field update coefficients
                        obj.coefEy_z(:,:) = obj.CourantNumber(1) * imp0;
                        obj.coefEx_z(:,:) = obj.CourantNumber(2) * imp0;
                        % set magnetic-field update coefficients */
                        obj.coefHz_y(:,:) = obj.CourantNumber(1) / imp0;
                        obj.coefHz_x(:,:) = obj.CourantNumber(2) / imp0;
                        
                        %grid coordinates
                        [obj.xEx, obj.yEx] = obj.getGridXYZ('Ex');
                        [obj.xEy, obj.yEy] = obj.getGridXYZ('Ey');
                        [obj.xHz, obj.yHz] = obj.getGridXYZ('Hz');
                    end
                    [obj.x,   obj.y]   = obj.getGridXYZ('');
            case 3 % 3D system
                % set electric-field update coefficients
                obj.coefEx_y(:,:,:) = obj.CourantNumber(3) * imp0;
                obj.coefEx_z(:,:,:) = obj.CourantNumber(2) * imp0;
                obj.coefEy_x(:,:,:) = obj.CourantNumber(3) * imp0;
                obj.coefEy_z(:,:,:) = obj.CourantNumber(1) * imp0;
                obj.coefEz_x(:,:,:) = obj.CourantNumber(2) * imp0;
                obj.coefEz_y(:,:,:) = obj.CourantNumber(1) * imp0;
                
                % set magnetic-field update coefficients */
                obj.coefHx_y(:,:,:) = obj.CourantNumber(3) / imp0;
                obj.coefHx_z(:,:,:) = obj.CourantNumber(2) / imp0;
                obj.coefHy_x(:,:,:) = obj.CourantNumber(3) / imp0;
                obj.coefHy_z(:,:,:) = obj.CourantNumber(1) / imp0;
                obj.coefHz_x(:,:,:) = obj.CourantNumber(2) / imp0;
                obj.coefHz_y(:,:,:) = obj.CourantNumber(1) / imp0;

                %grid coordinates
                [X,Y,Z]   = obj.getGridXYZ('');
                [obj.x,   obj.y,   obj.z] = ndgrid(X,Y,Z);
                
                [X,Y,Z]   = obj.getGridXYZ('Ex');
                [obj.xEx, obj.yEx, obj.zEx] = ndgrid(X,Y,Z);
                [X,Y,Z]   = obj.getGridXYZ('Ey');
                [obj.xEy, obj.yEy, obj.zEy] = ndgrid(X,Y,Z);
                [X,Y,Z]   = obj.getGridXYZ('Ez');
                [obj.xEz, obj.yEz, obj.zEz] = ndgrid(X,Y,Z);
                [X,Y,Z]   = obj.getGridXYZ('Hx');
                [obj.xHx, obj.yHx, obj.zHx] = ndgrid(X,Y,Z);
                [X,Y,Z]   = obj.getGridXYZ('Hy');
                [obj.xHy, obj.yHy, obj.zHy] = ndgrid(X,Y,Z);
                [X,Y,Z]   = obj.getGridXYZ('Hz');
                [obj.xHz, obj.yHz, obj.zHz] = ndgrid(X,Y,Z);

            end
            
            if ~isempty(Objects)
                for n_obj=1:length(Objects)
                    otmp = Objects{n_obj};
                    switch otmp.type
                        case 'sphere'
                            switch obj.Dimensionality
                                case 2
                                    test_func = @(X,Y) ((X-otmp.position(1))/otmp.radius(1)).^2 ...
                                        + ((Y-otmp.position(2))/otmp.radius(2)).^2 <= 1;

                                    % surrounding Ey nodes
                                    if obj.isTE 
                                        % surrounding Ez nodes
                                        [X,Y] = ndgrid(obj.xEz,obj.yEz);
                                        condition = test_func(X,Y);
                                        % lossy dielectric
                                        obj.coefEz(condition)   = obj.coefEz(condition)*(1-otmp.loss)/(1+otmp.loss);
                                        obj.coefEz_y(condition) = obj.coefEz_y(condition)/(otmp.RI)^2/(1+otmp.loss);
                                        if obj.isPML
                                            condition = test_func(X(2:end-1,2:end-1),Y(2:end-1,2:end-1));
                                            obj.epsilon_z(condition) = otmp.RI^2;
                                        end
                                    else
                                        % surrounding Ex nodes
                                        [X,Y] = ndgrid(obj.xEx,obj.yEx);
                                        condition = test_func(X,Y);
                                        % lossy dielectric
                                        obj.coefEx(condition)   = obj.coefEx(condition)*(1-otmp.loss)/(1+otmp.loss);
                                        obj.coefEx_y(condition) = obj.coefEx_y(condition)/(otmp.RI)^2/(1+otmp.loss);
                                        if obj.isPML
                                            condition = test_func(X(:,2:end-1),Y(:,2:end-1));
                                            obj.epsilon_x(condition) = otmp.RI^2;
                                        end

                                        % surrounding Ey nodes
                                        [X,Y] = ndgrid(obj.xEy,obj.yEy);
                                        condition = test_func(X,Y);
                                        % lossy dielectric
                                        obj.coefEy(condition)   = obj.coefEy(condition)*(1-otmp.loss)/(1+otmp.loss);
                                        obj.coefEy_x(condition) = obj.coefEy_x(condition)/(otmp.RI)^2/(1+otmp.loss);
                                        if obj.isPML
                                            condition = test_func(X(2:end-1,:),Y(2:end-1,:));
                                            obj.epsilon_y(condition) = otmp.RI^2;
                                        end
                                    end
                                case 3
                                    test_func = @(X,Y,Z) ((X-otmp.position(1))/otmp.radius(1)).^2 ...
                                        + ((Y-otmp.position(2))/otmp.radius(2)).^2 ...
                                        + ((Z-otmp.position(3))/otmp.radius(3)).^2 <= 1;
                                    % surrounding Ex nodes
                                    [X,Y,Z] = ndgrid(obj.xEx,obj.yEx,obj.zEx);
                                    condition = test_func(X,Y,Z);
                                    % lossy dielectric
                                    obj.coefEx(condition)   = obj.coefEx(condition)*(1-otmp.loss)/(1+otmp.loss);
                                    obj.coefEx_y(condition) = obj.coefEx_y(condition)/(otmp.RI)^2/(1+otmp.loss);
                                    obj.coefEx_z(condition) = obj.coefEx_z(condition)/(otmp.RI)^2/(1+otmp.loss);
                                    if obj.isPML
                                        condition = test_func(X(:,2:end-1,2:end-1),Y(:,2:end-1,2:end-1),Z(:,2:end-1,2:end-1));
                                        obj.epsilon_x(condition) = otmp.RI^2;
                                    end

                                    % surrounding Ey nodes
                                    [X,Y,Z] = ndgrid(obj.xEy,obj.yEy,obj.zEy);
                                    condition = test_func(X,Y,Z);
                                    % lossy dielectric
                                    obj.coefEy(condition)   = obj.coefEy(condition)*(1-otmp.loss)/(1+otmp.loss);
                                    obj.coefEy_x(condition) = obj.coefEy_x(condition)/(otmp.RI)^2/(1+otmp.loss);
                                    obj.coefEy_z(condition) = obj.coefEy_z(condition)/(otmp.RI)^2/(1+otmp.loss);
                                    if obj.isPML
                                        condition = test_func(X(2:end-1,:,2:end-1),Y(2:end-1,:,2:end-1),Z(2:end-1,:,2:end-1));
                                        obj.epsilon_y(condition) = otmp.RI^2;
                                    end

                                    % surrounding Ez nodes
                                    [X,Y,Z] = ndgrid(obj.xEz,obj.yEz,obj.zEz);
                                    condition = test_func(X,Y,Z);
                                    % lossy dielectric
                                    obj.coefEz(condition)   = obj.coefEz(condition)*(1-otmp.loss)/(1+otmp.loss);
                                    obj.coefEz_x(condition) = obj.coefEz_x(condition)/(otmp.RI)^2/(1+otmp.loss);
                                    obj.coefEz_y(condition) = obj.coefEz_y(condition)/(otmp.RI)^2/(1+otmp.loss);
                                    if obj.isPML
                                        condition = test_func(X(2:end-1,2:end-1,:),Y(2:end-1,2:end-1,:),Z(2:end-1,2:end-1,:));
                                        obj.epsilon_z(condition) = otmp.RI^2;
                                    end
                            end
                    end
                end
            end
            
            % initializations for PML
            if obj.isPML
                switch obj.Dimensionality
                    case 1
                        % Hy grid
                        sigmaxHy = zeros(size(obj.Hy));
                        condition = obj.PML.x1>obj.xHy;
                        sigmaxHy(condition) = ((obj.PML.x1-obj.xHy(condition))/obj.PML.Wx).^obj.PML.m./1; % /RI
                        condition = obj.PML.x2<obj.xHy;
                        sigmaxHy(condition) = ((obj.xHy(condition)-obj.PML.x2)/obj.PML.Wx).^obj.PML.m./1; % /RI
                        obj.PML.bxHy = exp(-sigmaxHy*(obj.PML.m+1)*0.8*c*obj.dt/obj.dx);

                        % Ez grid (2:end-1)
                        sigmaxEz = zeros(obj.Nx-1,1);
                        x_tmp = obj.xEz(2:end-1); condition = obj.PML.x1>x_tmp;
                        sigmaxEz(condition) = ((obj.PML.x1-x_tmp(condition))/obj.PML.Wx).^obj.PML.m./1; % /RI
                        condition = obj.PML.x2<x_tmp;
                        sigmaxEz(condition) = ((x_tmp(condition)-obj.PML.x2)/obj.PML.Wx).^obj.PML.m./1; % /RI
                        obj.PML.bxEz = exp(-sigmaxEz*(obj.PML.m+1)*0.8*c*obj.dt/obj.dx);
                    case 2
                        if obj.isTE
                            % Hx grid
                            sigmayHx = zeros(size(obj.Hx));
                            condition = obj.PML.y1>obj.yHx;
                            sigmayHx(condition) = ((obj.PML.y1-obj.yHx(condition))/obj.PML.Wy).^obj.PML.m./1; % /RI
                            condition = obj.PML.y2<obj.yHx;
                            sigmayHx(condition) = ((obj.yHy(condition)-obj.PML.y2)/obj.PML.Wy).^obj.PML.m./1; % /RI
                            obj.PML.byHx = exp(-sigmayHx*(obj.PML.m+1)*0.8*c*obj.dt/obj.dy);

                            % Hy grid
                            sigmaxHy = zeros(size(obj.Hy));
                            condition = obj.PML.x1>obj.xHy;
                            sigmaxHy(condition) = ((obj.PML.x1-obj.xHy(condition))/obj.PML.Wx).^obj.PML.m./1; % /RI
                            condition = obj.PML.x2<obj.xHy;
                            sigmaxHy(condition) = ((obj.xHy(condition)-obj.PML.x2)/obj.PML.Wx).^obj.PML.m./1; % /RI
                            obj.PML.bxHy = exp(-sigmaxHy*(obj.PML.m+1)*0.8*c*obj.dt/obj.dx);

                            % Ez grid (2:end-1,2:end-1)
                            sigmaxEz = zeros(obj.Nx-1,obj.Nx-1);
                            x_tmp = obj.xEz(2:end-1,2:end-1); condition = obj.PML.x1>x_tmp;
                            sigmaxEz(condition) = ((obj.PML.x1-x_tmp(condition))/obj.PML.Wx).^obj.PML.m./1; % /RI
                            condition = obj.PML.x2<x_tmp;
                            sigmaxEz(condition) = ((x_tmp(condition)-obj.PML.x2)/obj.PML.Wx).^obj.PML.m./1; % /RI
                            obj.PML.bxEz = exp(-sigmaxEz*(obj.PML.m+1)*0.8*c*obj.dt/obj.dx);

                            sigmayEz = zeros(obj.Nx-1,obj.Nx-1);
                            y_tmp = obj.yEz(2:end-1,2:end-1); condition = obj.PML.y1>y_tmp;
                            sigmayEz(condition) = ((obj.PML.y1-y_tmp(condition))/obj.PML.Wy).^obj.PML.m./1; % /RI
                            condition = obj.PML.y2<y_tmp;
                            sigmayEz(condition) = ((y_tmp(condition)-obj.PML.y2)/obj.PML.Wy).^obj.PML.m./1; % /RI
                            obj.PML.byEz = exp(-sigmayEz*(obj.PML.m+1)*0.8*c*obj.dt/obj.dy);
                        else
                            % Hz grid
                            sigmaxHz = zeros(size(obj.Hz));
                            condition = obj.PML.x1>obj.xHz;
                            sigmaxHz(condition) = ((obj.PML.x1-obj.xHz(condition))/obj.PML.Wx).^obj.PML.m./1; % /RI
                            condition = obj.PML.x2<obj.xHy;
                            sigmaxHz(condition) = ((obj.xHz(condition)-obj.PML.x2)/obj.PML.Wx).^obj.PML.m./1; % /RI
                            obj.PML.bxHz = exp(-sigmaxHz*(obj.PML.m+1)*0.8*c*obj.dt/obj.dx);

                            sigmayHz = zeros(size(obj.Hz));
                            condition = obj.PML.y1>obj.yHz;
                            sigmayHz(condition) = ((obj.PML.y1-obj.yHz(condition))/obj.PML.Wy).^obj.PML.m./1; % /RI
                            condition = obj.PML.y2<obj.yHz;
                            sigmayHz(condition) = ((obj.yHz(condition)-obj.PML.y2)/obj.PML.Wy).^obj.PML.m./1; % /RI
                            obj.PML.byHz = exp(-sigmayHz*(obj.PML.m+1)*0.8*c*obj.dt/obj.dy);

                            % Ex grid (:,2:end-1)
                            sigmayEx = zeros(obj.Nx,obj.Nx-1);
                            y_tmp = obj.yEx(:,2:end-1); condition = obj.PML.y1>y_tmp;
                            sigmayEx(condition) = ((obj.PML.y1-y_tmp(condition))/obj.PML.Wy).^obj.PML.m./1; % /RI
                            condition = obj.PML.y2<y_tmp;
                            sigmayEx(condition) = ((y_tmp(condition)-obj.PML.y2)/obj.PML.Wy).^obj.PML.m./1; % /RI
                            obj.PML.byEx = exp(-sigmayEx*(obj.PML.m+1)*0.8*c*obj.dt/obj.dy);
                            
                            % Ey grid (2:end-1,:)
                            sigmaxEy = zeros(obj.Nx-1,obj.Nx);
                            x_tmp = obj.xEy(2:end-1,:); condition = obj.PML.x1>x_tmp;
                            sigmaxEy(condition) = ((obj.PML.x1-x_tmp(condition))/obj.PML.Wx).^obj.PML.m./1; % /RI
                            condition = obj.PML.x2<x_tmp;
                            sigmaxEy(condition) = ((x_tmp(condition)-obj.PML.x2)/obj.PML.Wx).^obj.PML.m./1; % /RI
                            obj.PML.bxEy = exp(-sigmaxEy*(obj.PML.m+1)*0.8*c*obj.dt/obj.dx);
                        end
                    case 3
                        % Hx grid
                        sigmayHx = zeros(size(obj.Hx));
                        condition = obj.PML.y1>obj.yHx;
                        sigmayHx(condition) = ((obj.PML.y1-obj.yHx(condition))/obj.PML.Wy).^obj.PML.m./1; % /RI
                        condition = obj.PML.y2<obj.yHx;
                        sigmayHx(condition) = ((obj.yHy(condition)-obj.PML.y2)/obj.PML.Wy).^obj.PML.m./1; % /RI
                        obj.PML.byHx = exp(-sigmayHx*(obj.PML.m+1)*0.8*c*obj.dt/obj.dy);

                        sigmazHx = zeros(size(obj.Hx));
                        condition = obj.PML.z1>obj.zHx;
                        sigmazHx(condition) = ((obj.PML.z1-obj.zHx(condition))/obj.PML.Wz).^obj.PML.m./1; % /RI
                        condition = obj.PML.z2<obj.zHx;
                        sigmazHx(condition) = ((obj.zHx(condition)-obj.PML.z2)/obj.PML.Wz).^obj.PML.m./1; % /RI
                        obj.PML.bzHx = exp(-sigmazHx*(obj.PML.m+1)*0.8*c*obj.dt/obj.dz);

                        % Hy grid
                        sigmaxHy = zeros(size(obj.Hy));
                        condition = obj.PML.x1>obj.xHy;
                        sigmaxHy(condition) = ((obj.PML.x1-obj.xHy(condition))/obj.PML.Wx).^obj.PML.m./1; % /RI
                        condition = obj.PML.x2<obj.xHy;
                        sigmaxHy(condition) = ((obj.xHy(condition)-obj.PML.x2)/obj.PML.Wx).^obj.PML.m./1; % /RI
                        obj.PML.bxHy = exp(-sigmaxHy*(obj.PML.m+1)*0.8*c*obj.dt/obj.dx);

                        sigmazHy = zeros(size(obj.Hy));
                        condition = obj.PML.z1>obj.zHy;
                        sigmazHy(condition) = ((obj.PML.z1-obj.zHy(condition))/obj.PML.Wz).^obj.PML.m./1; % /RI
                        condition = obj.PML.z2<obj.zHy;
                        sigmazHy(condition) = ((obj.zHy(condition)-obj.PML.z2)/obj.PML.Wz).^obj.PML.m./1; % /RI
                        obj.PML.bzHy = exp(-sigmazHy*(obj.PML.m+1)*0.8*c*obj.dt/obj.dz);

                        % Hz grid
                        sigmaxHz = zeros(size(obj.Hz));
                        condition = obj.PML.x1>obj.xHz;
                        sigmaxHz(condition) = ((obj.PML.x1-obj.xHz(condition))/obj.PML.Wx).^obj.PML.m./1; % /RI
                        condition = obj.PML.x2<obj.xHy;
                        sigmaxHz(condition) = ((obj.xHz(condition)-obj.PML.x2)/obj.PML.Wx).^obj.PML.m./1; % /RI
                        obj.PML.bxHz = exp(-sigmaxHz*(obj.PML.m+1)*0.8*c*obj.dt/obj.dx);

                        sigmayHz = zeros(size(obj.Hz));
                        condition = obj.PML.y1>obj.yHz;
                        sigmayHz(condition) = ((obj.PML.y1-obj.yHz(condition))/obj.PML.Wy).^obj.PML.m./1; % /RI
                        condition = obj.PML.y2<obj.yHz;
                        sigmayHz(condition) = ((obj.yHz(condition)-obj.PML.y2)/obj.PML.Wy).^obj.PML.m./1; % /RI
                        obj.PML.byHz = exp(-sigmayHz*(obj.PML.m+1)*0.8*c*obj.dt/obj.dy);

                        % Ex grid (:,2:end-1,2:end-1)
                        sigmayEx = zeros(obj.Nx,obj.Nx-1,obj.Nx-1);
                        y_tmp = obj.yEx(:,2:end-1,2:end-1); condition = obj.PML.y1>y_tmp;
                        sigmayEx(condition) = ((obj.PML.y1-y_tmp(condition))/obj.PML.Wy).^obj.PML.m./1; % /RI
                        condition = obj.PML.y2<y_tmp;
                        sigmayEx(condition) = ((y_tmp(condition)-obj.PML.y2)/obj.PML.Wy).^obj.PML.m./1; % /RI
                        obj.PML.byEx = exp(-sigmayEx*(obj.PML.m+1)*0.8*c*obj.dt/obj.dy);

                        sigmazEx = zeros(obj.Nx,obj.Nx-1,obj.Nx-1);
                        z_tmp = obj.zEx(:,2:end-1,2:end-1); condition = obj.PML.z1>z_tmp;
                        sigmazEx(condition) = ((obj.PML.z1-z_tmp(condition))/obj.PML.Wz).^obj.PML.m./1; % /RI
                        condition = obj.PML.z2<z_tmp;
                        sigmazEx(condition) = ((z_tmp(condition)-obj.PML.z2)/obj.PML.Wz).^obj.PML.m./1; % /RI
                        obj.PML.bzEx = exp(-sigmazEx*(obj.PML.m+1)*0.8*c*obj.dt/obj.dz);

                        % Ey grid (2:end-1,:,2:end-1)
                        sigmaxEy = zeros(obj.Nx-1,obj.Nx,obj.Nx-1);
                        x_tmp = obj.xEy(2:end-1,:,2:end-1); condition = obj.PML.x1>x_tmp;
                        sigmaxEy(condition) = ((obj.PML.x1-x_tmp(condition))/obj.PML.Wx).^obj.PML.m./1; % /RI
                        condition = obj.PML.x2<x_tmp;
                        sigmaxEy(condition) = ((x_tmp(condition)-obj.PML.x2)/obj.PML.Wx).^obj.PML.m./1; % /RI
                        obj.PML.bxEy = exp(-sigmaxEy*(obj.PML.m+1)*0.8*c*obj.dt/obj.dx);

                        sigmazEy = zeros(obj.Nx-1,obj.Nx,obj.Nx-1);
                        z_tmp = obj.xEy(2:end-1,:,2:end-1); condition = obj.PML.z1>z_tmp;
                        sigmazEy(condition) = ((obj.PML.z1-z_tmp(condition))/obj.PML.Wz).^obj.PML.m./1; % /RI
                        condition = obj.PML.z2<z_tmp;
                        sigmazEy(condition) = ((z_tmp(condition)-obj.PML.z2)/obj.PML.Wz).^obj.PML.m./1; % /RI
                        obj.PML.bzEy = exp(-sigmazEy*(obj.PML.m+1)*0.8*c*obj.dt/obj.dz);

                        % Ez grid (2:end-1,2:end-1,:)
                        sigmaxEz = zeros(obj.Nx-1,obj.Nx-1,obj.Nx);
                        x_tmp = obj.xEz(2:end-1,2:end-1,:); condition = obj.PML.x1>x_tmp;
                        sigmaxEz(condition) = ((obj.PML.x1-x_tmp(condition))/obj.PML.Wx).^obj.PML.m./1; % /RI
                        condition = obj.PML.x2<x_tmp;
                        sigmaxEz(condition) = ((x_tmp(condition)-obj.PML.x2)/obj.PML.Wx).^obj.PML.m./1; % /RI
                        obj.PML.bxEz = exp(-sigmaxEz*(obj.PML.m+1)*0.8*c*obj.dt/obj.dx);

                        sigmayEz = zeros(obj.Nx-1,obj.Nx-1,obj.Nx);
                        y_tmp = obj.yEz(2:end-1,2:end-1,:); condition = obj.PML.y1>y_tmp;
                        sigmayEz(condition) = ((obj.PML.y1-y_tmp(condition))/obj.PML.Wy).^obj.PML.m./1; % /RI
                        condition = obj.PML.y2<y_tmp;
                        sigmayEz(condition) = ((y_tmp(condition)-obj.PML.y2)/obj.PML.Wy).^obj.PML.m./1; % /RI
                        obj.PML.byEz = exp(-sigmayEz*(obj.PML.m+1)*0.8*c*obj.dt/obj.dy);
                end
            end
        
        end
        
        %% update magnetic field
        function obj = updateH(obj)
            global mu0
            switch obj.Dimensionality
                case 1
                    if obj.isPML
                        DxEz = diff(obj.Ez,[],1)/obj.dx;

                        %Update of the PML-Matrices
                        obj.PML.QxEz = obj.PML.bxHy.*(obj.PML.QxEz+DxEz)-DxEz;

                        obj.Hy(:) = obj.Hy ...
                            + (obj.dt/mu0./obj.mu_y).*(-DxEz-obj.PML.QxEz); % + JHy(2:end,1:end-1)
                    else
                        obj.Hy(:) = obj.coefHy.*obj.Hy + obj.coefHy_z.*diff(obj.Ez,[],1);
                    end
                case 2
                    if obj.isTE  % TM polarization
                        if obj.isPML
                            DyEz = diff(obj.Ez,[],2)/obj.dy;
                            DxEz = diff(obj.Ez,[],1)/obj.dx;

                            %Update of the PML-Matrices
                            obj.PML.QyEz = obj.PML.byHx.*(obj.PML.QyEz+DyEz)-DyEz;
                            obj.PML.QxEz = obj.PML.bxHy.*(obj.PML.QxEz+DxEz)-DxEz;

                            obj.Hx(:,:) = obj.Hx ...
                                + (obj.dt/mu0./obj.mu_x).*(DyEz+obj.PML.QyEz); % + JHx(1:end-1,2:end) 
                            obj.Hy(:,:) = obj.Hy ...
                                + (obj.dt/mu0./obj.mu_y).*(-DxEz-obj.PML.QxEz); % + JHy(2:end,1:end-1)
                        else
                            obj.Hx(:,:) = obj.coefHx.*obj.Hx - obj.coefHx_z.*diff(obj.Ez,[],2);
                            obj.Hy(:,:) = obj.coefHy.*obj.Hy + obj.coefHy_z.*diff(obj.Ez,[],1);
                        end
                    else % TE polarization
                        if obj.isPML
                            DxEy = diff(obj.Ey,[],1)/obj.dx;
                            DyEx = diff(obj.Ex,[],2)/obj.dy;

                            %Update of the PML-Matrices
                            obj.PML.QxEy = obj.PML.bxHz.*(obj.PML.QxEy+DxEy)-DxEy;
                            obj.PML.QyEx = obj.PML.byHz.*(obj.PML.QyEx+DyEx)-DyEx;

                            obj.Hz(:,:) = obj.Hz ...
                                + (obj.dt/mu0./obj.mu_z).*(DyEx+obj.PML.QyEx - DxEy-obj.PML.QxEy); % + JHz(2:end,1:end-1)

                        else
                            obj.Hz(:,:) = obj.coefHz.*obj.Hz + obj.coefHz_x.*diff(obj.Ex,[],2) - obj.coefHz_y.*diff(obj.Ey,[],1);
                        end
                    end
                case 3
                    if obj.isPML
                        DyEz = diff(obj.Ez,[],2)/obj.dy;
                        DzEy = diff(obj.Ey,[],3)/obj.dz;
                        DxEz = diff(obj.Ez,[],1)/obj.dx;
                        DzEx = diff(obj.Ex,[],3)/obj.dz;
                        DxEy = diff(obj.Ey,[],1)/obj.dx;
                        DyEx = diff(obj.Ex,[],2)/obj.dy;

                        %Update of the PML-Matrices
                        obj.PML.QyEz = obj.PML.byHx.*(obj.PML.QyEz+DyEz)-DyEz;
                        obj.PML.QzEy = obj.PML.bzHx.*(obj.PML.QzEy+DzEy)-DzEy;
                        obj.PML.QxEz = obj.PML.bxHy.*(obj.PML.QxEz+DxEz)-DxEz;
                        obj.PML.QzEx = obj.PML.bzHy.*(obj.PML.QzEx+DzEx)-DzEx;
                        obj.PML.QxEy = obj.PML.bxHz.*(obj.PML.QxEy+DxEy)-DxEy;
                        obj.PML.QyEx = obj.PML.byHz.*(obj.PML.QyEx+DyEx)-DyEx;
                        
                        obj.Hx(:,:,:) = obj.Hx ...
                            + (obj.dt/mu0./obj.mu_x).*(DyEz+obj.PML.QyEz - DzEy-obj.PML.QzEy); % + JHx(1:end-1,2:end) 
                        obj.Hy(:,:,:) = obj.Hy ...
                            + (obj.dt/mu0./obj.mu_y).*(DzEx+obj.PML.QzEx - DxEz-obj.PML.QxEz); % + JHy(2:end,1:end-1)
                        obj.Hz(:,:,:) = obj.Hz ...
                            + (obj.dt/mu0./obj.mu_z).*(DyEx+obj.PML.QyEx - DxEy-obj.PML.QxEy); % + JHz(2:end,1:end-1)
                        
                    else
                        obj.Hx(:,:,:) = obj.coefHx.*obj.Hx + obj.coefHx_y.*diff(obj.Ey,[],3) - obj.coefHx_z.*diff(obj.Ez,[],2);
                        obj.Hy(:,:,:) = obj.coefHy.*obj.Hy + obj.coefHy_z.*diff(obj.Ez,[],1) - obj.coefHy_x.*diff(obj.Ex,[],3);
                        obj.Hz(:,:,:) = obj.coefHz.*obj.Hz + obj.coefHz_x.*diff(obj.Ex,[],2) - obj.coefHz_y.*diff(obj.Ey,[],1);
                    end
            end
        end

        %% update electric field
        function obj = updateE(obj)
            global eps0
            switch obj.Dimensionality
                case 1
                    if obj.isPML
                        DxHy = diff(obj.Hy,[],1)/obj.dx;

                        %Update of the PML-Matrices
                        obj.PML.QxHy = obj.PML.bxEz.*(obj.PML.QxHy+DxHy)-DxHy;

                        obj.Ez(obj.m2) = obj.Ez(obj.m2) ...
                            +(obj.dt/eps0./obj.epsilon_z).*(-DxHy-obj.PML.QxHy); % + JEz(2:end,2:end)

                    else
                        obj.Ez(obj.m2) = obj.coefEz(obj.m2).*obj.Ez(obj.m2) + obj.coefEz_y(obj.m2).*diff(obj.Hy,1);
                    end
                case 2
                    if obj.isTE  % TM polarization
                        if obj.isPML
                            DxHy = diff(obj.Hy(:,obj.n2),[],1)/obj.dx;
                            DyHx = diff(obj.Hx(obj.m2,:),[],2)/obj.dy;

                            %Update of the PML-Matrices
                            obj.PML.QxHy = obj.PML.bxEz.*(obj.PML.QxHy+DxHy)-DxHy;
                            obj.PML.QyHx = obj.PML.byEz.*(obj.PML.QyHx+DyHx)-DyHx;

                            obj.Ez(obj.m2,obj.n2) = obj.Ez(obj.m2,obj.n2,:) ...
                                +(obj.dt/eps0./obj.epsilon_z).*(DyHx+obj.PML.QyHx - DxHy-obj.PML.QxHy); % + JEz(2:end,2:end)
                        else
                            obj.Ez(obj.m2,obj.n2) = obj.coefEz(obj.m2,obj.n2).*obj.Ez(obj.m2,obj.n2) ...
                                + obj.coefEz_y(obj.m2,obj.n2).*diff(obj.Hy(:,obj.n2),[],1)  ...
                                - obj.coefEz_x(obj.m2,obj.n2).*diff(obj.Hx(obj.m2,:),[],2);
                        end
                    else % TE polarization
                        if obj.isPML
                            DyHz = diff(obj.Hz,[],2)/obj.dy;
                            DxHz = diff(obj.Hz,[],1)/obj.dx;

                            %Update of the PML-Matrices
                            obj.PML.QyHz = obj.PML.byEx.*(obj.PML.QyHz+DyHz)-DyHz;
                            obj.PML.QxHz = obj.PML.bxEy.*(obj.PML.QxHz+DxHz)-DxHz;

                            obj.Ex(:,obj.n2,obj.p2) = obj.Ex(:,obj.n2,obj.p2) ...
                                +(obj.dt/eps0./obj.epsilon_x).*(-DyHz-obj.PML.QyHz); % + JEz(2:end,2:end)
                            obj.Ey(obj.m2,:,obj.p2) = obj.Ey(obj.m2,:,obj.p2) ...
                                +(obj.dt/eps0./obj.epsilon_y).*(DxHz+obj.PML.QxHz); % + JEz(2:end,2:end)
                        else
                            obj.Ex(:,obj.n2) = obj.coefEx(:,obj.n2).*obj.Ex(:,obj.n2) ...
                                + obj.coefEx_z(:,obj.n2).*diff(obj.Hz,[],2);
                            obj.Ey(obj.m2,:) = obj.coefEy(obj.m2,:).*obj.Ey(obj.m2,:) ...
                                - obj.coefEy_z(obj.m2,:).*diff(obj.Hz,[],1);
                        end
                    end
                case 3
                    if obj.isPML
                        DyHz = diff(obj.Hz(:,:,obj.p2),[],2)/obj.dy;
                        DzHy = diff(obj.Hy(:,obj.n2,:),[],3)/obj.dz;
                        DxHz = diff(obj.Hz(:,:,obj.p2),[],1)/obj.dx;
                        DzHx = diff(obj.Hx(obj.m2,:,:),[],3)/obj.dz;
                        DxHy = diff(obj.Hy(:,obj.n2,:),[],1)/obj.dx;
                        DyHx = diff(obj.Hx(obj.m2,:,:),[],2)/obj.dy;

                        %Update of the PML-Matrices
                        obj.PML.QyHz = obj.PML.byEx.*(obj.PML.QyHz+DyHz)-DyHz;
                        obj.PML.QzHy = obj.PML.bzEx.*(obj.PML.QzHy+DzHy)-DzHy;
                        obj.PML.QxHz = obj.PML.bxEy.*(obj.PML.QxHz+DxHz)-DxHz;
                        obj.PML.QzHx = obj.PML.bzEy.*(obj.PML.QzHx+DzHx)-DzHx;
                        obj.PML.QxHy = obj.PML.bxEz.*(obj.PML.QxHy+DxHy)-DxHy;
                        obj.PML.QyHx = obj.PML.byEz.*(obj.PML.QyHx+DyHx)-DyHx;
    
                        obj.Ex(:,obj.n2,obj.p2) = obj.Ex(:,obj.n2,obj.p2) ...
                            +(obj.dt/eps0./obj.epsilon_x).*(DzHy+obj.PML.QzHy - DyHz-obj.PML.QyHz); % + JEz(2:end,2:end)
                        obj.Ey(obj.m2,:,obj.p2) = obj.Ey(obj.m2,:,obj.p2) ...
                            +(obj.dt/eps0./obj.epsilon_y).*(DxHz+obj.PML.QxHz - DzHx-obj.PML.QzHx); % + JEz(2:end,2:end)
                        obj.Ez(obj.m2,obj.n2,:) = obj.Ez(obj.m2,obj.n2,:) ...
                            +(obj.dt/eps0./obj.epsilon_z).*(DyHx+obj.PML.QyHx - DxHy-obj.PML.QxHy); % + JEz(2:end,2:end)
                        
                    else
                        obj.Ex(:,obj.n2,obj.p2) = obj.coefEx(:,obj.n2,obj.p2).*obj.Ex(:,obj.n2,obj.p2) ...
                            + obj.coefEx_z(:,obj.n2,obj.p2).*diff(obj.Hz(:,:,obj.p2),[],2) ...
                            - obj.coefEx_y(:,obj.n2,obj.p2).*diff(obj.Hy(:,obj.n2,:),[],3);
                        obj.Ey(obj.m2,:,obj.p2) = obj.coefEy(obj.m2,:,obj.p2).*obj.Ey(obj.m2,:,obj.p2) ...
                            + obj.coefEy_x(obj.m2,:,obj.p2).*diff(obj.Hx(obj.m2,:,:),[],3) ...
                            - obj.coefEy_z(obj.m2,:,obj.p2).*diff(obj.Hz(:,:,obj.p2),[],1);
                        obj.Ez(obj.m2,obj.n2,:) = obj.coefEz(obj.m2,obj.n2,:).*obj.Ez(obj.m2,obj.n2,:) ...
                            + obj.coefEz_y(obj.m2,obj.n2,:).*diff(obj.Hy(:,obj.n2,:),[],1) ...
                            - obj.coefEz_x(obj.m2,obj.n2,:).*diff(obj.Hx(obj.m2,:,:),[],2);
                    end
            end
        end
        
        %% Absorbing boundary conditioms
        function obj = initABC(obj,type) % initialize ABC
            switch obj.Dimensionality
                case 1
                    switch type
                        case '1abc' % first order ABC
                            obj.ABC.type = 1;
                            obj.ABC.oldEz_x	= [0, 0];
                        otherwise % second order ABC
                            obj.ABC.type = 2;
                            obj.ABC.oldEz_x = zeros(3,2,2);
                    end                   
%                     obj.ABC.coef = (obj.CourantNumber-1)/(obj.CourantNumber+1);
                    % calculate coefficient on left-end and right-end of grid
%                     obj.ABC.coef = coefABC(obj,'coefEz_y','coefHy_z',obj.ABC.type); % (3,start:end)
                    obj.ABC.coef = coefABC2(obj.CourantNumber(1),obj.ABC.type); % (3,start:end)
                case 2
                    % allocate memory for ABC arrays (second order)
                    switch type
                        case '1abc' % first order ABC
                            obj.ABC.type = 1;
                            if obj.isTE
                                obj.ABC.oldEz_x	= zeros(2, obj.Ny+1); % (m(start:end),n(all)) <- along x-axis
                                obj.ABC.oldEz_y	= zeros(obj.Nx+1, 2); % (m(all),n(start:end)) <- along y-axis
                                % calculate coefficient on left-end and right-end of grid
                                %{
                                obj.ABC.coef(:,2) = coefABC(obj,'coefEz_x','coefHx_z',obj.ABC.type); % ((start:end),y) <- along y-axis
                                obj.ABC.coef(:,1) = coefABC(obj,'coefEz_y','coefHy_z',obj.ABC.type); % ((start:end),x) <- along x-axis
                                %}
                                obj.ABC.coef(:,2) = coefABC2(obj.CourantNumber(2),obj.ABC.type); % ((start:end),y) <- along y-axis
                                obj.ABC.coef(:,1) = coefABC2(obj.CourantNumber(1),obj.ABC.type); % ((start:end),x) <- along x-axis
                            else
                                obj.ABC.oldEy_x	= zeros(2, obj.Ny); % (m(start:end),n(all)) <- along x-axis
                                obj.ABC.oldEx_y	= zeros(obj.Nx, 2); % (m(all),n(start:end)) <- along y-axis
                                % calculate coefficient on left-end and right-end of grid
                                %{
                                obj.ABC.coef(:,2) = coefABC(obj,'coefEx_z','coefHz_x',obj.ABC.type); % ((start:end),y) <- along y-axis
                                obj.ABC.coef(:,1) = coefABC(obj,'coefEy_z','coefHz_y',obj.ABC.type); % ((start:end),x) <- along x-axis
                                %}
                                obj.ABC.coef(:,2) = coefABC2(obj.CourantNumber(2),obj.ABC.type); % ((start:end),y) <- along y-axis
                                obj.ABC.coef(:,1) = coefABC2(obj.CourantNumber(1),obj.ABC.type); % ((start:end),x) <- along x-axis
                            end
                        otherwise % second order ABC
                            obj.ABC.type = 2;
                            if obj.isTE
                                obj.ABC.oldEz_x	= zeros(3,obj.Ny+1,2,2); % (m(1,2,3),n(all),prev:next,start:end) <- along x-axis
                                obj.ABC.oldEz_y	= zeros(obj.Nx+1,3,2,2); % (m(all),n(1,2,3),prev:next,start:end) <- along x-axis
                                % calculate coefficient on left-end and right-end of grid
                                %{
                                obj.ABC.coef(:,:,2) = coefABC(obj,'coefEz_x','coefHx_z',obj.ABC.type); % (3,start:end,y) <- along y-axis
                                obj.ABC.coef(:,:,1) = coefABC(obj,'coefEz_y','coefHy_z',obj.ABC.type); % (3,start:end,x) <- along x-axis
                                %}
                                obj.ABC.coef(:,:,2) = coefABC2(obj.CourantNumber(2),obj.ABC.type); % (3,start:end,y) <- along y-axis
                                obj.ABC.coef(:,:,1) = coefABC2(obj.CourantNumber(1),obj.ABC.type); % (3,start:end,x) <- along x-axis
                            else
                                obj.ABC.oldEy_x	= zeros(3,obj.Ny,2,2); % (m(1,2,3),n(all),prev:next,start:end) <- along x-axis
                                obj.ABC.oldEx_y	= zeros(obj.Nx,3,2,2); % (m(all),n(1,2,3),prev:next,start:end) <- along x-axis
                                % calculate coefficient on left-end and right-end of grid
                                %{
                                obj.ABC.coef(:,:,2) = coefABC(obj,'coefEx_z','coefHz_x',obj.ABC.type); % (3,start:end) <- along y-axis
                                obj.ABC.coef(:,:,1) = coefABC(obj,'coefEy_z','coefHz_y',obj.ABC.type); % (3,start:end) <- along x-axis
                                %}
                                obj.ABC.coef(:,:,2) = coefABC2(obj.CourantNumber(2),obj.ABC.type); % (3,start:end) <- along y-axis
                                obj.ABC.coef(:,:,1) = coefABC2(obj.CourantNumber(1),obj.ABC.type); % (3,start:end) <- along x-axis
                            end
                    end
                case 3
                    switch type
                        case '1abc' % first order ABC
                            obj.ABC.type = 1;
                            % allocate memory for ABC arrays
                            obj.ABC.oldEy_x = zeros(2,obj.Ny, obj.Nz+1); % (m(start:end,n,p)
                            obj.ABC.oldEz_x = zeros(2,obj.Ny+1,   obj.Nz);

                            obj.ABC.oldEx_y = zeros(obj.Nx, 2,obj.Nz+1);
                            obj.ABC.oldEz_y = zeros(obj.Nx+1,   2,obj.Nz);

                            obj.ABC.oldEx_z = zeros(obj.Nx, obj.Ny+1,  2);
                            obj.ABC.oldEy_z = zeros(obj.Nx+1,   obj.Ny,2);
                            % calculate ABC coefficients (first order)
                            %{
                            obj.ABC.coef(:,6) = coefABC(obj,'coefEy_x','coefHx_y',obj.ABC.type); % (2,start:end,z=3) <- Ey along z-axis
                            obj.ABC.coef(:,5) = coefABC(obj,'coefEx_y','coefHy_x',obj.ABC.type); % (2,start:end,z=3) <- Ex along z-axis
                            obj.ABC.coef(:,4) = coefABC(obj,'coefEz_x','coefHx_z',obj.ABC.type); % (2,start:end,y=2) <- Ez along y-axis
                            obj.ABC.coef(:,3) = coefABC(obj,'coefEx_z','coefHz_x',obj.ABC.type); % (2,start:end,y=2) <- Ex along y-axis
                            obj.ABC.coef(:,2) = coefABC(obj,'coefEz_y','coefHy_z',obj.ABC.type); % (2,start:end,x=1) <- Ez along x-axis
                            obj.ABC.coef(:,1) = coefABC(obj,'coefEy_z','coefHz_y',obj.ABC.type); % (2,start:end,x=1) <- Ey along x-axis
                            %}
                            obj.ABC.coef(:,6) = coefABC2(obj.CourantNumber(3),obj.ABC.type); % (2,start:end,z=3) <- Ey along z-axis
                            obj.ABC.coef(:,5) = coefABC2(obj.CourantNumber(3),obj.ABC.type); % (2,start:end,z=3) <- Ex along z-axis
                            obj.ABC.coef(:,4) = coefABC2(obj.CourantNumsber(2),obj.ABC.type); % (2,start:end,y=2) <- Ez along y-axis
                            obj.ABC.coef(:,3) = coefABC2(obj.CourantNumber(2),obj.ABC.type); % (2,start:end,y=2) <- Ex along y-axis
                            obj.ABC.coef(:,2) = coefABC2(obj.CourantNumber(1),obj.ABC.type); % (2,start:end,x=1) <- Ez along x-axis
                            obj.ABC.coef(:,1) = coefABC2(obj.CourantNumber(1),obj.ABC.type); % (2,start:end,x=1) <- Ey along x-axis
%                             obj.ABC.coef(2,1:2) = (obj.CourantNumber(1:2)-1)./(obj.CourantNumber(1:2)+1);
%                             obj.ABC.coef(1,1:2) = obj.ABC.coef(2,1:3); % (start:end),(x,y,z)
                        otherwise % second order ABC
                            obj.ABC.type = 2;
                            % allocate memory for ABC arrays
                            obj.ABC.oldEy_x = zeros(3,obj.Ny,obj.Nz+1,  2,2); % (m,n,p,prev:next,start:end)
                            obj.ABC.oldEz_x = zeros(3,obj.Ny+1,  obj.Nz,2,2);

                            obj.ABC.oldEx_y = zeros(obj.Nx, 3,obj.Nz+1,  2,2);
                            obj.ABC.oldEz_y = zeros(obj.Nx+1,   3,obj.Nz,2,2);

                            obj.ABC.oldEx_z = zeros(obj.Nx, obj.Ny+1,  3,2,2);
                            obj.ABC.oldEy_z = zeros(obj.Nx+1,   obj.Ny,3,2,2);
                            % calculate ABC coefficients (second order)
                            %{
                            obj.ABC.coef(:,:,6) = coefABC(obj,'coefEy_x','coefHx_y',obj.ABC.type); % (3,start:end,z=3) <- Ey along z-axis
                            obj.ABC.coef(:,:,5) = coefABC(obj,'coefEx_y','coefHy_x',obj.ABC.type); % (3,start:end,z=3) <- Ex along z-axis
                            obj.ABC.coef(:,:,4) = coefABC(obj,'coefEz_x','coefHx_z',obj.ABC.type); % (3,start:end,y=2) <- Ez along y-axis
                            obj.ABC.coef(:,:,3) = coefABC(obj,'coefEx_z','coefHz_x',obj.ABC.type); % (3,start:end,y=2) <- Ex along y-axis
                            obj.ABC.coef(:,:,2) = coefABC(obj,'coefEz_y','coefHy_z',obj.ABC.type); % (3,start:end,x=1) <- Ez along x-axis
                            obj.ABC.coef(:,:,1) = coefABC(obj,'coefEy_z','coefHz_y',obj.ABC.type); % (3,start:end,x=1) <- Ey along x-axis
                            %}

                            obj.ABC.coef(:,:,6) = coefABC2(obj.CourantNumber(3),obj.ABC.type); % (3,start:end,z=3) <- Ey along z-axis
                            obj.ABC.coef(:,:,5) = coefABC2(obj.CourantNumber(3),obj.ABC.type); % (3,start:end,z=3) <- Ex along z-axis
                            obj.ABC.coef(:,:,4) = coefABC2(obj.CourantNumber(2),obj.ABC.type); % (3,start:end,y=2) <- Ez along y-axis
                            obj.ABC.coef(:,:,3) = coefABC2(obj.CourantNumber(2),obj.ABC.type); % (3,start:end,y=2) <- Ex along y-axis
                            obj.ABC.coef(:,:,2) = coefABC2(obj.CourantNumber(1),obj.ABC.type); % (3,start:end,x=1) <- Ez along x-axis
                            obj.ABC.coef(:,:,1) = coefABC2(obj.CourantNumber(1),obj.ABC.type); % (3,start:end,x=1) <- Ey along x-axis

%                             obj.ABC.coef(1:3,2,1:3) = [-(1./obj.CourantNumber-2+obj.CourantNumber);
%                                                          -2*(obj.CourantNumber-1./obj.CourantNumber);
%                                                           4*(obj.CourantNumber+1./obj.CourantNumber)]./(1./obj.CourantNumber([1 1 1],:)+2+obj.CourantNumber([1 1 1],:));
%                             obj.ABC.coef(:,1,1:3) = obj.ABC.coef(:,2,1:3); % m,(start:end),(x,y,z)
                    end
            end
        end
        
        function obj = applyABC(obj)
            if isempty(obj.ABC) % check if initABC() has been called
                fprintf('\nUse initABC() method before applyABC()!...')
%                 exit(-1);
            end
            switch obj.Dimensionality
                case 1
                    switch obj.ABC.type
                        case 0    % simple ABC for left side of grid
                            obj.Ez(1) = obj.Ez(2);
                        case 1
                            % ABC for left and right sides of grid
                            obj.Ez([1,obj.Nx]) = obj.ABC.oldEz_x + obj.ABC.coef.*(obj.Ez([2,obj.Nx+1-2]) - obj.Ez([1,obj.Nx]));
                            % update stored fields
                            obj.ABC.oldEz_x = obj.Ez([2,obj.Nx-1]);
                        case 2 % Second-order ABC.
                            % ABC for left side of grid
                            obj.Ez(1) = ...
                                obj.ABC.coef(1,1) * (obj.ABC.oldEz_x(1,2,1) +  obj.Ez(3)) ...
                              + obj.ABC.coef(2,1) * (obj.ABC.oldEz_x(1,1,1) + obj.ABC.oldEz_x(3,1,1) ...
                                                   - obj.ABC.oldEz_x(2,2,1) -  obj.Ez(2)) ...
                              + obj.ABC.coef(3,1) *  obj.ABC.oldEz_x(2,1,1) - obj.ABC.oldEz_x(3,2,1);
                            % update stored fields
                            obj.ABC.oldEz_x1(:,2) = obj.ABC.oldEz_x(:,1,1);
                            obj.ABC.oldEz_x1(:,1) = obj.Ez(1:3);
                            
                            % ABC for right side of grid
                            m = obj.Nx;
                            obj.Ez(m) = ...
                                obj.ABC.coef(1,2) * (obj.ABC.oldEz_x(1,2,2) +  obj.Ez(m-2)) ...
                              + obj.ABC.coef(2,2) * (obj.ABC.oldEz_x(1,1,2) + obj.ABC.oldEz_x(3,1,2) ...
                                                   - obj.ABC.oldEz_x(2,2,2) -  obj.Ez(m-1)) ...
                              + obj.ABC.coef(3,2) *  obj.ABC.oldEz_x(2,1,2) - obj.ABC.oldEz_x(3,2,2);
                            % update stored fields
                            obj.ABC.oldEz_x(:,2,2) = obj.ABC.oldEz_x(:,1,2);
                            obj.ABC.oldEz_x(:,1,2) = obj.Ez(m-(0:2));
                    end
                case 2
                    switch obj.ABC.type
                        case 1
                            % first-order ABC.
                            if obj.isTE
                                % ABC at x1 (start) and x2 (end)
                                m_1 = 1;
                                m_2 = obj.Nx;
                                obj.Ez(m_1,:)	= obj.ABC.oldEz_x(1,:) + obj.ABC.coef(1)*(obj.Ez(m_1+1,:) - obj.Ez(m_1,:));
                                obj.Ez(m_2,:)	= obj.ABC.oldEz_x(2,:) + obj.ABC.coef(2)*(obj.Ez(m_2-1,:) - obj.Ez(m_2,:));
                                % update stored fields
                                obj.ABC.oldEz_x = obj.Ez([2,obj.Nx-1],:);

                                % ABC at y1 (start) and y2 (end)
                                n_1 = 1;
                                n_2 = obj.Ny;
                                obj.Ez(:,n_1)	= obj.ABC.oldEz_y(:,1) + obj.ABC.coef(1)*(obj.Ez(:,n_1+1) - obj.Ez(:,n_1));
                                obj.Ez(:,n_2)	= obj.ABC.oldEz_y(:,2) + obj.ABC.coef(2)*(obj.Ez(:,n_2-1) - obj.Ez(:,n_2));
                                % update stored fields
                                obj.ABC.oldEz_y = obj.Ez(:,[2,obj.Ny-1]);
                            else
                                % ABC at x1 (start) and x2 (end)
                                m_1 = 1;
                                m_2 = obj.Nx;
                                obj.Ey(m_1,:)	= obj.ABC.oldEy_x(1,:) + obj.ABC.coef(1)*(obj.Ey(m_1+1,:) - obj.Ey(m_1,:));
                                obj.Ey(m_2,:)	= obj.ABC.oldEy_x(2,:) + obj.ABC.coef(2)*(obj.Ey(m_2-1,:) - obj.Ey(m_2,:));
                                % update stored fields
                                obj.ABC.oldEy_x = obj.Ey([2,obj.Nx-1],:);

                                % ABC at y1 (start) and y2 (end)
                                n_1 = 1;
                                n_2 = obj.Ny;
                                obj.Ex(:,n_1)	= obj.ABC.oldEx_y(:,1) + obj.ABC.coef(1)*(obj.Ex(:,n_1+1) - obj.Ex(:,n_1));
                                obj.Ex(:,n_2)	= obj.ABC.oldEx_y(:,2) + obj.ABC.coef(2)*(obj.Ex(:,n_2-1) - obj.Ex(:,n_2));
                                % update stored fields
                                obj.ABC.oldEx_y = obj.Ex(:,[2,obj.Ny-1]);
                            end
                        case 2
                            % Second-order ABC.
                            if obj.isTE
                                % ABC at x1 side of grid
                                obj.Ez(1,:) = ...
                                    obj.ABC.coef(1,1)*(obj.ABC.oldEz_x(1,:,2,1) +  obj.Ez(3,:)) ...
                                  + obj.ABC.coef(2,1)*(obj.ABC.oldEz_x(1,:,1,1) + obj.ABC.oldEz_x(3,:,1,1) ...
                                                     - obj.ABC.oldEz_x(2,:,2,1) -  obj.Ez(2,:)) ...
                                  + obj.ABC.coef(3,1)* obj.ABC.oldEz_x(2,:,1,1) - obj.ABC.oldEz_x(3,:,2,1);
                                % memorize old fields
                                obj.ABC.oldEz_x(:,:,2,1) = obj.ABC.oldEz_x(:,:,1,1);
                                obj.ABC.oldEz_x(:,:,1,1) = obj.Ez(1:3,:);

                                % ABC at x2 side of grid
                                m = obj.Nx;
                                obj.Ez(m,:) = ...
                                    obj.ABC.coef(1,2)*(obj.ABC.oldEz_x(1,:,2,2) +  obj.Ez(m-2,:)) ...
                                  + obj.ABC.coef(2,2)*(obj.ABC.oldEz_x(1,:,1,2) + obj.ABC.oldEz_x(3,:,1,2) ...
                                                     - obj.ABC.oldEz_x(2,:,2,2) -  obj.Ez(m-1,:)) ...
                                  + obj.ABC.coef(3,2)* obj.ABC.oldEz_x(2,:,1,2) - obj.ABC.oldEz_x(3,:,2,2);
                                % memorize old fields
                                obj.ABC.oldEz_x(:,:,2,2) = obj.ABC.oldEz_x(:,:,1,2);
                                obj.ABC.oldEz_x(:,:,1,2) = obj.Ez(m-(0:2),:);

                                % ABC at y1 of grid
                                obj.Ez(:,1) = ...
                                    obj.ABC.coef(1,1)*(obj.ABC.oldEz_y(:,1,2,1) +  obj.Ez(:,3)) ...
                                  + obj.ABC.coef(2,1)*(obj.ABC.oldEz_y(:,1,1,1) + obj.ABC.oldEz_y(:,3,1,1) ...
                                                     - obj.ABC.oldEz_y(:,2,2,1) -  obj.Ez(:,2)) ...
                                  + obj.ABC.coef(3,1)* obj.ABC.oldEz_y(:,2,1,1) - obj.ABC.oldEz_y(:,3,2,1);
                                % memorize old fields
                                obj.ABC.oldEz_y(:,:,2,1) = obj.ABC.oldEz_y(:,:,1,1);
                                obj.ABC.oldEz_y(:,:,1,1) = obj.Ez(:,1:3);

                                % ABC at y2 of grid
                                n = obj.Ny;
                                obj.Ez(:,n) = ...
                                    obj.ABC.coef(1,2)*(obj.ABC.oldEz_y(:,1,2,2) +  obj.Ez(:,n-2)) ...
                                  + obj.ABC.coef(2,2)*(obj.ABC.oldEz_y(:,1,1,2) + obj.ABC.oldEz_y(:,3,1,2) ...
                                                     - obj.ABC.oldEz_y(:,2,2,2) -  obj.Ez(:,n-1)) ...
                                  + obj.ABC.coef(3,2)* obj.ABC.oldEz_y(:,2,1,2) - obj.ABC.oldEz_y(:,3,2,2);
                                % memorize old fields
                                obj.ABC.oldEz_y(:,:,2,2) = obj.ABC.oldEz_y(:,:,1,2);
                                obj.ABC.oldEz_y(:,:,1,2) = obj.Ez(:,n-(0:2));
                            else
                                % ABC at x1 of grid
                                obj.Ey(1,:) = ...
                                    obj.ABC.coef(1,1)*(obj.ABC.oldEy_x(1,:,2,1) +  obj.Ey(3,:)) ...
                                  + obj.ABC.coef(2,1)*(obj.ABC.oldEy_x(1,:,1,1) + obj.ABC.oldEy_x(3,:,1,1) ...
                                                     - obj.ABC.oldEy_x(2,:,2,1) -  obj.Ey(2,:)) ...
                                  + obj.ABC.coef(3,1)* obj.ABC.oldEy_x(2,:,1,1) - obj.ABC.oldEy_x(3,:,2,1);
                                % memorize old fields
                                obj.ABC.oldEy_x(:,:,2,1) = obj.ABC.oldEy_x(:,:,1,1);
                                obj.ABC.oldEy_x(:,:,1,1) = obj.Ey(1:3,:);

                                % ABC at x2 of grid
                                m = obj.Nx;
                                obj.Ey(m,:) = ...
                                    obj.ABC.coef(1,2)*(obj.ABC.oldEy_x(1,:,2,2) +  obj.Ey(m-2,:)) ...
                                  + obj.ABC.coef(2,2)*(obj.ABC.oldEy_x(1,:,1,2) + obj.ABC.oldEy_x(3,:,1,2) ...
                                                     - obj.ABC.oldEy_x(2,:,2,2) -  obj.Ey(m-1,:)) ...
                                  + obj.ABC.coef(3,2)* obj.ABC.oldEy_x(2,:,1,2) - obj.ABC.oldEy_x(3,:,2,2);
                                % memorize old fields
                                obj.ABC.oldEy_x(:,:,2,2) = obj.ABC.oldEy_x(:,:,1,2);
                                obj.ABC.oldEy_x(:,:,1,2) = obj.Ey(m-(0:2),:);

                                % ABC at y1 of grid
                                obj.Ex(:,1) = ...
                                    obj.ABC.coef(1,1)*(obj.ABC.oldEx_y(:,1,2,1) +  obj.Ex(:,3)) ...
                                  + obj.ABC.coef(2,1)*(obj.ABC.oldEx_y(:,1,1,1) + obj.ABC.oldEx_y(:,3,1,1) ...
                                                     - obj.ABC.oldEx_y(:,2,2,1) -  obj.Ex(:,2)) ...
                                  + obj.ABC.coef(3,1)* obj.ABC.oldEx_y(:,2,1,1) - obj.ABC.oldEx_y(:,3,2,1);
                                % memorize old fields
                                obj.ABC.oldEx_y(:,:,2,1) = obj.ABC.oldEx_y(:,:,1,1);
                                obj.ABC.oldEx_y(:,:,1,1) = obj.Ex(:,1:3);

                                % ABC at y2 of grid
                                n = obj.Ny;
                                obj.Ex(:,n) = ...
                                    obj.ABC.coef(1,2)*(obj.ABC.oldEx_y(:,1,2,2) +  obj.Ex(:,n-2)) ...
                                  + obj.ABC.coef(2,2)*(obj.ABC.oldEx_y(:,1,1,2) + obj.ABC.oldEx_y(:,3,1,2) ...
                                                     - obj.ABC.oldEx_y(:,2,2,2) -  obj.Ex(:,n-1)) ...
                                  + obj.ABC.coef(3,2)* obj.ABC.oldEx_y(:,2,1,2) - obj.ABC.oldEx_y(:,3,2,2);
                                % memorize old fields
                                obj.ABC.oldEx_y(:,:,2,2) = obj.ABC.oldEx_y(:,:,1,2);
                                obj.ABC.oldEx_y(:,:,1,2) = obj.Ex(:,n-(0:2));
                            end
                    end
                case 3
                    switch obj.ABC.type
                        case 1  % first-order ABC.
                            % ABC at x1 (start) and x2 (end)
                            m = 1;
                            obj.Ey(m,:,:) = obj.ABC.oldEy_x(1,:,:) + obj.ABC.coef(1,1)*(obj.Ey(m+1,:,:) - obj.Ey(m,:,:));
                            obj.Ez(m,:,:) = obj.ABC.oldEz_x(1,:,:) + obj.ABC.coef(1,2)*(obj.Ez(m+1,:,:) - obj.Ez(m,:,:));
                            m = obj.Nx;
                            obj.Ey(m,:,:) = obj.ABC.oldEy_x(2,:,:) + obj.ABC.coef(2,1)*(obj.Ey(m-1,:,:) - obj.Ey(m,:,:));
                            obj.Ez(m,:,:) = obj.ABC.oldEz_x(2,:,:) + obj.ABC.coef(2,2)*(obj.Ez(m-1,:,:) - obj.Ez(m,:,:));
                            obj.ABC.oldEy_x(:,:,:) = obj.Ey([2,obj.Nx-1],:,:);
                            obj.ABC.oldEz_x(:,:,:) = obj.Ez([2,obj.Nx-1],:,:);

                            % ABC at y1 (start) and y2 (end)
                            n = 1;
                            obj.Ez(:,n,:) = obj.ABC.oldEz_y(:,1,:) + obj.ABC.coef(1,3)*(obj.Ez(:,n+1,:) - obj.Ez(:,n,:));
                            obj.Ex(:,n,:) = obj.ABC.oldEx_y(:,1,:) + obj.ABC.coef(1,4)*(obj.Ex(:,n+1,:) - obj.Ex(:,n,:));
                            n = obj.Ny;
                            obj.Ez(:,n,:) = obj.ABC.oldEz_y(:,2,:) + obj.ABC.coef(2,3)*(obj.Ez(:,n-1,:) - obj.Ez(:,n,:));
                            obj.Ex(:,n,:) = obj.ABC.oldEx_y(:,2,:) + obj.ABC.coef(2,4)*(obj.Ex(:,n-1,:) - obj.Ex(:,n,:));
                            obj.ABC.oldEz_y(:,:,:) = obj.Ez(:,[2, obj.Ny-1],:);
                            obj.ABC.oldEx_y(:,:,:) = obj.Ex(:,[2, obj.Ny-1],:);

                            % ABC at z1 (start) and z2 (end)
                            p = 1;
                            obj.Ex(:,:,p) = obj.ABC.oldEx_z(:,:,1) + obj.ABC.coef(1,5)*(obj.Ex(:,:,p+1) - obj.Ex(:,:,p));
                            obj.Ey(:,:,p) = obj.ABC.oldEy_z(:,:,1) + obj.ABC.coef(1,6)*(obj.Ey(:,:,p+1) - obj.Ey(:,:,p));
                            p = obj.Nz;
                            obj.Ex(:,:,p) = obj.ABC.oldEx_z(:,:,2) + obj.ABC.coef(2,5)*(obj.Ex(:,:,p-1) - obj.Ex(:,:,p));
                            obj.Ey(:,:,p) = obj.ABC.oldEy_z(:,:,2) + obj.ABC.coef(2,6)*(obj.Ey(:,:,p-1) - obj.Ey(:,:,p));
                            obj.ABC.oldEx_z(:,:,:) = obj.Ex(:,:,[2, obj.Nz-1]);
                            obj.ABC.oldEy_z(:,:,:) = obj.Ey(:,:,[2, obj.Nz-1]);
                        case 2 % Second-order ABC.
                            % ABC at x1 side of grid
                            obj.Ey(1,:,:) = ...
                                obj.ABC.coef(1,1,1)*(obj.ABC.oldEy_x(1,:,:,2,1) +  obj.Ey(3,:,:)) ...
                              + obj.ABC.coef(2,1,1)*(obj.ABC.oldEy_x(1,:,:,1,1) + obj.ABC.oldEy_x(3,:,:,1,1) ...
                                                   - obj.ABC.oldEy_x(2,:,:,2,1) -  obj.Ey(2,:,:)) ...
                              + obj.ABC.coef(3,1,1)* obj.ABC.oldEy_x(2,:,:,1,1) - obj.ABC.oldEy_x(3,:,:,2,1);
                            obj.Ez(1,:,:) = ...
                                obj.ABC.coef(1,1,2)*(obj.ABC.oldEz_x(1,:,:,2,1) +  obj.Ez(3,:,:)) ...
                              + obj.ABC.coef(2,1,2)*(obj.ABC.oldEz_x(1,:,:,1,1) + obj.ABC.oldEz_x(3,:,:,1,1) ...
                                                   - obj.ABC.oldEz_x(2,:,:,2,1) -  obj.Ez(2,:,:)) ...
                              + obj.ABC.coef(3,1,2)* obj.ABC.oldEz_x(2,:,:,1,1) - obj.ABC.oldEz_x(3,:,:,2,1);
                            % memorize old fields
                            obj.ABC.oldEy_x(:,:,:,2,1) = obj.ABC.oldEy_x(:,:,:,1,1);
                            obj.ABC.oldEy_x(:,:,:,1,1) = obj.Ey(1:3,:,:);
                            obj.ABC.oldEz_x(:,:,:,2,1) = obj.ABC.oldEz_x(:,:,:,1,1);
                            obj.ABC.oldEz_x(:,:,:,1,1) = obj.Ez(1:3,:,:);

                            % ABC at x2 side of grid
                            m = obj.Nx;
                            obj.Ey(m,:,:) = ...
                                obj.ABC.coef(1,2,1)*(obj.ABC.oldEy_x(1,:,:,2,2) +  obj.Ey(m-2,:,:)) ...
                              + obj.ABC.coef(2,2,1)*(obj.ABC.oldEy_x(1,:,:,1,2) + obj.ABC.oldEy_x(3,:,:,1,2) ...
                                                   - obj.ABC.oldEy_x(2,:,:,2,2) -  obj.Ey(m-1,:,:)) ...
                              + obj.ABC.coef(3,2,1)* obj.ABC.oldEy_x(2,:,:,1,2) - obj.ABC.oldEy_x(3,:,:,2,2);
                            obj.Ez(m,:,:) = ...
                                obj.ABC.coef(1,2,2)*(obj.ABC.oldEz_x(1,:,:,2,2) +  obj.Ez(m-2,:,:)) ...
                              + obj.ABC.coef(2,2,2)*(obj.ABC.oldEz_x(1,:,:,1,2) + obj.ABC.oldEz_x(3,:,:,1,2) ...
                                                   - obj.ABC.oldEz_x(2,:,:,2,2) -  obj.Ez(m-1,:,:)) ...
                              + obj.ABC.coef(3,2,2)* obj.ABC.oldEz_x(2,:,:,1,2) - obj.ABC.oldEz_x(3,:,:,2,2);
                            % memorize old fields
                            obj.ABC.oldEy_x(:,:,:,2,2) = obj.ABC.oldEy_x(:,:,:,1,2);
                            obj.ABC.oldEy_x(:,:,:,1,2) = obj.Ey(m-(0:2),:,:);
                            obj.ABC.oldEz_x(:,:,:,2,2) = obj.ABC.oldEz_x(:,:,:,1,2);
                            obj.ABC.oldEz_x(:,:,:,1,2) = obj.Ez(m-(0:2),:,:);

                            % ABC at y1 of grid
                            obj.Ex(:,1,:) = ...
                                obj.ABC.coef(1,1,3)*(obj.ABC.oldEx_y(:,1,:,2,1) +  obj.Ex(:,3,:)) ...
                              + obj.ABC.coef(2,1,3)*(obj.ABC.oldEx_y(:,1,:,1,1) + obj.ABC.oldEx_y(:,3,:,1,1) ...
                                                   - obj.ABC.oldEx_y(:,2,:,2,1) -  obj.Ex(:,2,:)) ...
                              + obj.ABC.coef(3,1,3)* obj.ABC.oldEx_y(:,2,:,1,1) - obj.ABC.oldEx_y(:,3,:,2,1);
                            obj.Ez(:,1,:) = ...
                                obj.ABC.coef(1,1,4)*(obj.ABC.oldEz_y(:,1,:,2,1) +  obj.Ez(:,3,:)) ...
                              + obj.ABC.coef(2,1,4)*(obj.ABC.oldEz_y(:,1,:,1,1) + obj.ABC.oldEz_y(:,3,:,1,1) ...
                                                   - obj.ABC.oldEz_y(:,2,:,2,1) -  obj.Ez(:,2,:)) ...
                              + obj.ABC.coef(3,1,4)* obj.ABC.oldEz_y(:,2,:,1,1) - obj.ABC.oldEz_y(:,3,:,2,1);
                            % memorize old fields
                            obj.ABC.oldEx_y(:,:,:,2,1) = obj.ABC.oldEx_y(:,:,:,1,1);
                            obj.ABC.oldEx_y(:,:,:,1,1) = obj.Ex(:,1:3,:);
                            obj.ABC.oldEz_y(:,:,:,2,1) = obj.ABC.oldEz_y(:,:,:,1,1);
                            obj.ABC.oldEz_y(:,:,:,1,1) = obj.Ez(:,1:3,:);

                            % ABC at y2 of grid
                            n = obj.Ny;
                            obj.Ex(:,n,:) = ...
                                obj.ABC.coef(1,2,3)*(obj.ABC.oldEx_y(:,1,:,2,2) +  obj.Ex(:,n-2,:)) ...
                              + obj.ABC.coef(2,2,3)*(obj.ABC.oldEx_y(:,1,:,1,2) + obj.ABC.oldEx_y(:,3,:,1,2) ...
                                                   - obj.ABC.oldEx_y(:,2,:,2,2) -  obj.Ex(:,n-1,:)) ...
                              + obj.ABC.coef(3,2,3)* obj.ABC.oldEx_y(:,2,:,1,2) - obj.ABC.oldEx_y(:,3,:,2,2);
                            obj.Ez(:,n,:) = ...
                                obj.ABC.coef(1,2,4)*(obj.ABC.oldEz_y(:,1,:,2,2) +  obj.Ez(:,n-2,:)) ...
                              + obj.ABC.coef(2,2,4)*(obj.ABC.oldEz_y(:,1,:,1,2) + obj.ABC.oldEz_y(:,3,:,1,2) ...
                                                   - obj.ABC.oldEz_y(:,2,:,2,2) -  obj.Ez(:,n-1,:)) ...
                              + obj.ABC.coef(3,2,4)* obj.ABC.oldEz_y(:,2,:,1,2) - obj.ABC.oldEz_y(:,3,:,2,2);
                            % memorize old fields
                            obj.ABC.oldEx_y(:,:,:,2,2) = obj.ABC.oldEx_y(:,:,:,1,2);
                            obj.ABC.oldEx_y(:,:,:,1,2) = obj.Ex(:,n-(0:2),:);
                            obj.ABC.oldEz_y(:,:,:,2,2) = obj.ABC.oldEz_y(:,:,:,1,2);
                            obj.ABC.oldEz_y(:,:,:,1,2) = obj.Ez(:,n-(0:2),:);

                            % ABC at z1 of grid
                            obj.Ex(:,:,1) = ...
                                obj.ABC.coef(1,1,5)*(obj.ABC.oldEx_z(:,:,1,2,1) +  obj.Ex(:,:,3)) ...
                              + obj.ABC.coef(2,1,5)*(obj.ABC.oldEx_z(:,:,1,1,1) + obj.ABC.oldEx_z(:,:,3,1,1) ...
                                                   - obj.ABC.oldEx_z(:,:,2,2,1) -  obj.Ex(:,:,2)) ...
                              + obj.ABC.coef(3,1,5)* obj.ABC.oldEx_z(:,:,2,1,1) - obj.ABC.oldEx_z(:,:,3,2,1);
                            obj.Ey(:,:,1) = ...
                                obj.ABC.coef(1,1,6)*(obj.ABC.oldEy_z(:,:,1,2,1) +  obj.Ey(:,:,3)) ...
                              + obj.ABC.coef(2,1,6)*(obj.ABC.oldEy_z(:,:,1,1,1) + obj.ABC.oldEy_z(:,:,3,1,1) ...
                                                   - obj.ABC.oldEy_z(:,:,2,2,1) -  obj.Ey(:,:,2)) ...
                              + obj.ABC.coef(3,1,6)* obj.ABC.oldEy_z(:,:,2,1,1) - obj.ABC.oldEy_z(:,:,3,2,1);
                            % memorize old fields
                            obj.ABC.oldEx_z(:,:,:,2,1) = obj.ABC.oldEx_z(:,:,:,1,1);
                            obj.ABC.oldEx_z(:,:,:,1,1) = obj.Ex(:,:,1:3);
                            obj.ABC.oldEy_z(:,:,:,2,1) = obj.ABC.oldEy_z(:,:,:,1,1);
                            obj.ABC.oldEy_z(:,:,:,1,1) = obj.Ey(:,:,1:3);

                            % ABC at z2 of grid
                            p = obj.Nz;
                            obj.Ex(:,:,p) = ...
                                obj.ABC.coef(1,2,5)*(obj.ABC.oldEx_z(:,:,1,2,2) +  obj.Ex(:,:,p-2)) ...
                              + obj.ABC.coef(2,2,5)*(obj.ABC.oldEx_z(:,:,1,1,2) + obj.ABC.oldEx_z(:,:,3,1,2) ...
                                                   - obj.ABC.oldEx_z(:,:,2,2,2) -  obj.Ex(:,:,p-1)) ...
                              + obj.ABC.coef(3,2,5)* obj.ABC.oldEx_z(:,:,2,1,2) - obj.ABC.oldEx_z(:,:,3,2,2);
                            obj.Ey(:,:,p) = ...
                                obj.ABC.coef(1,2,6)*(obj.ABC.oldEy_z(:,:,1,2,2) +  obj.Ey(:,:,p-2)) ...
                              + obj.ABC.coef(2,2,6)*(obj.ABC.oldEy_z(:,:,1,1,2) + obj.ABC.oldEy_z(:,:,3,1,2) ...
                                                   - obj.ABC.oldEy_z(:,:,2,2,2) -  obj.Ey(:,:,p-1)) ...
                              + obj.ABC.coef(3,2,6)* obj.ABC.oldEy_z(:,:,2,1,2) - obj.ABC.oldEy_z(:,:,3,2,2);
                            % memorize old fields
                            obj.ABC.oldEx_z(:,:,:,2,2) = obj.ABC.oldEx_z(:,:,:,1,2);
                            obj.ABC.oldEx_z(:,:,:,1,2) = obj.Ex(:,:,p-(0:2));
                            obj.ABC.oldEy_z(:,:,:,2,2) = obj.ABC.oldEy_z(:,:,:,1,2);
                            obj.ABC.oldEy_z(:,:,:,1,2) = obj.Ey(:,:,p-(0:2));
                    end
            end
        end
        
        %% TFSF
        function obj = initTFSF(obj,TFSFregion,TFSFsource) % initialize TFSF boundary
            obj.TFSF.region	= TFSFregion;
            obj.TFSF.source	= TFSFsource;
            obj.TFSF.source.omega = TFSFsource.omega*obj.dt;
            switch obj.Dimensionality
                case {2,3}
                    % for TE Ey -> -Ez, Hz -> Hy
                    obj.TFSF.fields1D = FDTD(obj.Nx+1,obj.FDTDspace(1,:));                 % create the 1D FDTD grid
                    obj.TFSF.fields1D.initGrid([],obj.CourantNumber);	% initialize the 1D FDTD grid with Courant number for 2d/3D model
                    obj.TFSF.fields1D.initTFSF(1,obj.TFSF.source);
                    obj.TFSF.fields1D.initABC('2abc');                   % initialize ABC

            end
        end
        
        function obj = updateTFSF(obj,time)
            if isempty(obj.TFSF) % check if initTFSF() has been called
                fprintf('\nUse initTFSF() method before updateTFSF()!...')
%                 exit(-1);
            end

            switch obj.Dimensionality
                case 1
                    % correct Hy adjacent to TFSF boundary
                    obj.Hy(obj.TFSF.region) = obj.Hy(obj.TFSF.region) ...
                        - obj.coefHy_z(obj.TFSF.region)*SourceFunction(time, obj.TFSF.source,obj.CourantNumber(1));
                    % correct Ez adjacent to TFSF boundary (time+0.5-(-0.5))
                    obj.Ez(obj.TFSF.region+1) = obj.Ez(obj.TFSF.region+1) ...
                        + SourceFunction(time+1,obj.TFSF.source,obj.CourantNumber(1));
                case 2
                    if obj.isTE
                        % correct H adjacent to TFSF boundary
                        
                        %****** constant x faces (scattered-field nodes) ******
                        m_1 = obj.TFSF.region(1,1)-1;
                        m_2 = obj.TFSF.region(1,2);
                        n = obj.TFSF.region(2,1):obj.TFSF.region(2,2);
                        % correct Hy at firstX-1/2 by subtracting Ez_inc
                        obj.Hy(m_1,n) = obj.Hy(m_1,n) - obj.coefHy_z(m_1,n).*obj.TFSF.fields1D.Ez(m_1+1);
                        % correct Hy at lastX + 1/2 by adding Ez_inc
                        obj.Hy(m_2,n) = obj.Hy(m_2,n) + obj.coefHy_z(m_2,n).*obj.TFSF.fields1D.Ez(m_2);

                        %**** constant y faces (scattered-field nodes) ****
                        m = obj.TFSF.region(1,1):obj.TFSF.region(1,2);
                        n_1 = obj.TFSF.region(2,1)-1;
                        n_2 = obj.TFSF.region(2,2);
                        % correct Hx at firstY-1/2 by adding Ez_inc
                        obj.Hx(m,n_1) = obj.Hx(m,n_1) + obj.coefHx_z(m,n_1).*obj.TFSF.fields1D.Ez(m);
                        % correct Hx at lastY+1/2 by subtracting Ez_inc
                        obj.Hx(m,n_2) = obj.Hx(m,n_2) - obj.coefHx_z(m,n_2).*obj.TFSF.fields1D.Ez(m);
                    
                        % update 1D
                        obj.TFSF.fields1D.updateH();         % update 1D magnetic field
                        obj.TFSF.fields1D.updateTFSF(time);	% in 1D, apply TFSF
                        obj.TFSF.fields1D.updateE();         % update 1D electric field
%                         obj.TFSF.fields1D.Ez(3) = SourceFunction(time, obj.TFSF.source,obj.CourantNumber); % set source node
                        obj.TFSF.fields1D.applyABC();        % apply ABC

                        % correct E adjacent to TFSF boundary

                        %**** constant x faces (total-field nodes) ****
                        m_1 = obj.TFSF.region(1,1);
                        m_2 = obj.TFSF.region(1,2);
                        n = obj.TFSF.region(2,1):obj.TFSF.region(2,2);
                        % correct Ez at firstX face by subtracting Hy_inc
                        obj.Ez(m_1,n) = obj.Ez(m_1,n) - obj.coefEz_x(m_1,n).*obj.TFSF.fields1D.Hy(m_1-1);
                        % correct Ez at lastX face by adding Hy_inc
                        obj.Ez(m_2,n) = obj.Ez(m_2,n) + obj.coefEz_x(m_2,n).*obj.TFSF.fields1D.Hy(m_2);

                        % no need to correct Ez along top and bottom since
                        % incident Hx is zero

                    else
                        % correct Hz along x-axis (start)
                        m = obj.TFSF.region(1,1)-1;
                        n = obj.TFSF.region(2,1):obj.TFSF.region(2,2)-1;
                        obj.Hz(m,n) = obj.Hz(m,n) + obj.coefHz_y(m,n).*(-obj.TFSF.fields1D.Ez(m+1)); % Ez <- (-Ey)
                        % correct Hz along x-axis (end)
                        m = obj.TFSF.region(1,2);
                        obj.Hz(m,n) = obj.Hz(m,n) - obj.coefHz_y(m,n).*(-obj.TFSF.fields1D.Ez(m)); % Ez <- (-Ey)
                        
                        % update 1D
                        obj.TFSF.fields1D.updateH();         % update 1D magnetic field
                        obj.TFSF.fields1D.updateTFSF(time);	% in 1D, apply TFSF
                        obj.TFSF.fields1D.updateE();         % update 1D electric field
%                         obj.TFSF.fields1D.Ey(1) = SourceFunction(time,obj.TFSF.source,obj.CourantNumber); %<-Ey in 1D, set source node
                        obj.TFSF.fields1D.applyABC();        % apply ABC

                        % correct Ex along y-axis (start)
                        m = obj.TFSF.region(1,1):obj.TFSF.region(1,2)-1;
                        n = obj.TFSF.region(2,1);
                        obj.Ex(m,n) = obj.Ex(m,n) - obj.coefEx_z(m,n).*obj.TFSF.fields1D.Hy(m); % Hy <- Hz
                        % correct Ex along y-axis (end)
                        n = obj.TFSF.region(2,2);
                        obj.Ex(m,n) = obj.Ex(m,n) + obj.coefEx_z(m,n).*obj.TFSF.fields1D.Hy(m); % Hy <- Hz

                        % correct Ey field along x-axis (start)
                        m = obj.TFSF.region(1,1);
                        n = obj.TFSF.region(2,1):obj.TFSF.region(2,2)-1;
                        obj.Ey(m,n) = obj.Ey(m,n) + obj.coefEy_z(m,n).*obj.TFSF.fields1D.Hy(m-1); % Hy <- Hz
                        % correct Ey field along y-axis (end)
                        m = obj.TFSF.region(1,2);
                        obj.Ey(m,n) = obj.Ey(m,n) - obj.coefEy_z(m,n).*obj.TFSF.fields1D.Hy(m); % Hy <- Hz

                        % no need to correct Ex along top and bottom since
                        % incident Ex is zero

                    end
                case 3
                    % correct H adjacent to TFSF boundary
                    
                    %****** constant x faces (scattered-field nodes) ******
                    m_1 = obj.TFSF.region(1,1)-1;
                    m_2 = obj.TFSF.region(1,2);
                    n = obj.TFSF.region(2,1):obj.TFSF.region(2,2);
                    p = obj.TFSF.region(3,1):obj.TFSF.region(3,2)-1;
                    % correct Hy at firstX-1/2 by subtracting Ez_inc
                    obj.Hy(m_1,n,p) = obj.Hy(m_1,n,p) - obj.coefHy_x(m_1,n,p).*obj.TFSF.fields1D.Ez(m_1+1);
                    % correct Hy at lastX + 1/2 by adding Ez_inc
                    obj.Hy(m_2,n,p) = obj.Hy(m_2,n,p) + obj.coefHy_x(m_2,n,p).*obj.TFSF.fields1D.Ez(m_2);
                    
                    %**** constant y faces (scattered-field nodes) ****
                    m = obj.TFSF.region(1,1):obj.TFSF.region(1,2);
                    n_1 = obj.TFSF.region(2,1)-1;
                    n_2 = obj.TFSF.region(2,2);
                    for p = obj.TFSF.region(3,1):obj.TFSF.region(3,2)-1;
                        % correct Hx at firstY-1/2 by adding Ez_inc
                        obj.Hx(m,n_1,p) = obj.Hx(m,n_1,p) + obj.coefHx_y(m,n_1,p).*obj.TFSF.fields1D.Ez(m);
                        % correct Hx at lastY+1/2 by subtracting Ez_inc
                        obj.Hx(m,n_2,p) = obj.Hx(m,n_2,p) - obj.coefHx_y(m,n_2,p).*obj.TFSF.fields1D.Ez(m);
                    end
                    
                    %**** constant z faces (scattered-field nodes) ****
                    % nothing to correct on this face
                    
                    %**** update the fields in the auxiliary 1D grid ****
                    obj.TFSF.fields1D.updateH();         % update 1D magnetic field
                    obj.TFSF.fields1D.updateTFSF(time);	% in 1D, apply TFSF
                    obj.TFSF.fields1D.updateE();         % update 1D electric field
%                     obj.TFSF.fields1D.Ez(1) = SourceFunction(time,obj.TFSF.source,obj.CourantNumber); %<-Ey in 1D, set source node
                    obj.TFSF.fields1D.applyABC();        % apply ABC

                    % correct E adjacent to TFSF boundary
                    
                    %**** constant x faces (total-field nodes) ****
                    m_1 = obj.TFSF.region(1,1);
                    m_2 = obj.TFSF.region(1,2);
                    n = obj.TFSF.region(2,1):obj.TFSF.region(2,2);
                    p = obj.TFSF.region(3,1):obj.TFSF.region(3,2)-1;
                    % correct Ez at firstX face by subtracting Hy_inc
                    obj.Ez(m_1,n,p) = obj.Ez(m_1,n,p) - obj.coefEz_x(m_1,n,p).*obj.TFSF.fields1D.Hy(m_1-1);
                    % correct Ez at lastX face by adding Hy_inc
                    obj.Ez(m_2,n,p) = obj.Ez(m_2,n,p) + obj.coefEz_x(m_2,n,p).*obj.TFSF.fields1D.Hy(m_2);

                    %**** constant y faces (total-field nodes) ****
                    % nothing to correct on this face
                    
                    %**** constant z faces (total-field nodes) ****
                    m = obj.TFSF.region(1,1):obj.TFSF.region(1,2)-1;
                    p_1 = obj.TFSF.region(3,1);
                    p_2 = obj.TFSF.region(3,2);
                    for n = obj.TFSF.region(2,1):obj.TFSF.region(2,2);
                        % correct Ex at firstZ face by adding Hy_inc
                        obj.Ex(m,n,p_1) = obj.Ex(m,n,p_1) + obj.coefEx_y(m,n,p_1).*obj.TFSF.fields1D.Hy(m);
                        % correct Ex at lastZ face by subtracting Hy_inc
                        obj.Ex(m,n,p_2) = obj.Ex(m,n,p_2) - obj.coefEx_y(m,n,p_2).*obj.TFSF.fields1D.Hy(m);
                    end
            end
        end
        
        %% Scattering field
        function obj = initSCAT(obj,SCATregion) % initialize TFSF boundary
            obj.SCAT.region	= SCATregion;
        end

        function obj = updateSCAT(obj)
            if isempty(obj.TFSF) % check if initTFSF() has been called
                fprintf('\nUse initTFSF() method before updateTFSF()!...')
%                 exit(-1);
            end

            switch obj.Dimensionality
                case 1
                    % correct Hy adjacent to TFSF boundary
                    obj.Hy(obj.TFSF.region) = obj.Hy(obj.TFSF.region) ...
                        - obj.coefHy_z(obj.TFSF.region)*SourceFunction(time, obj.TFSF.source,obj.CourantNumber(1));
                    % correct Ez adjacent to TFSF boundary (time+0.5-(-0.5))
                    obj.Ez(obj.TFSF.region+1) = obj.Ez(obj.TFSF.region+1) ...
                        + SourceFunction(time+1,obj.TFSF.source,obj.CourantNumber(1));
                case 2
                    if obj.isTE
                        % correct H adjacent to TFSF boundary
                        
                        %****** constant x faces (scattered-field nodes) ******
                        m_1 = obj.TFSF.region(1,1)-1;
                        m_2 = obj.TFSF.region(2,1);
                        n = obj.TFSF.region(1,2):obj.TFSF.region(2,2);
                        % correct Hy at firstX-1/2 by subtracting Ez_inc
                        obj.Hy(m_1,n) = obj.Hy(m_1,n) - obj.coefHy_z(m_1,n).*obj.TFSF.fields1D.Ez(m_1+1);
                        % correct Hy at lastX + 1/2 by adding Ez_inc
                        obj.Hy(m_2,n) = obj.Hy(m_2,n) + obj.coefHy_z(m_2,n).*obj.TFSF.fields1D.Ez(m_2);

                        %**** constant y faces (scattered-field nodes) ****
                        m = obj.TFSF.region(1,1):obj.TFSF.region(2,1);
                        n_1 = obj.TFSF.region(1,2)-1;
                        n_2 = obj.TFSF.region(2,2);
                        % correct Hx at firstY-1/2 by adding Ez_inc
                        obj.Hx(m,n_1) = obj.Hx(m,n_1) + obj.coefHx_z(m,n_1).*obj.TFSF.fields1D.Ez(m);
                        % correct Hx at lastY+1/2 by subtracting Ez_inc
                        obj.Hx(m,n_2) = obj.Hx(m,n_2) - obj.coefHx_z(m,n_2).*obj.TFSF.fields1D.Ez(m);
                    
                        % update 1D
                        obj.TFSF.fields1D.updateH();         % update 1D magnetic field
                        obj.TFSF.fields1D.updateTFSF(time);	% in 1D, apply TFSF
                        obj.TFSF.fields1D.updateE();         % update 1D electric field
%                         obj.TFSF.fields1D.Ez(3) = SourceFunction(time, obj.TFSF.source,obj.CourantNumber); % set source node
                        obj.TFSF.fields1D.applyABC();        % apply ABC

                        % correct E adjacent to TFSF boundary

                        %**** constant x faces (total-field nodes) ****
                        m_1 = obj.TFSF.region(1,1);
                        m_2 = obj.TFSF.region(2,1);
                        n = obj.TFSF.region(1,2):obj.TFSF.region(2,2);
                        % correct Ez at firstX face by subtracting Hy_inc
                        obj.Ez(m_1,n) = obj.Ez(m_1,n) - obj.coefEz_x(m_1,n).*obj.TFSF.fields1D.Hy(m_1-1);
                        % correct Ez at lastX face by adding Hy_inc
                        obj.Ez(m_2,n) = obj.Ez(m_2,n) + obj.coefEz_x(m_2,n).*obj.TFSF.fields1D.Hy(m_2);

                        % no need to correct Ez along top and bottom since
                        % incident Hx is zero

                    else
                        % correct Hz along x-axis (start)
                        m = obj.TFSF.region(1,1)-1;
                        n = obj.TFSF.region(1,2):obj.TFSF.region(2,2)-1;
                        obj.Hz(m,n) = obj.Hz(m,n) + obj.coefHz_y(m,n).*(-obj.TFSF.fields1D.Ez(m+1)); % Ez <- (-Ey)
                        % correct Hz along x-axis (end)
                        m = obj.TFSF.region(2,1);
                        obj.Hz(m,n) = obj.Hz(m,n) - obj.coefHz_y(m,n).*(-obj.TFSF.fields1D.Ez(m)); % Ez <- (-Ey)
                        
                        % update 1D
                        obj.TFSF.fields1D.updateH();         % update 1D magnetic field
                        obj.TFSF.fields1D.updateTFSF(time);	% in 1D, apply TFSF
                        obj.TFSF.fields1D.updateE();         % update 1D electric field
%                         obj.TFSF.fields1D.Ey(1) = SourceFunction(time,obj.TFSF.source,obj.CourantNumber); %<-Ey in 1D, set source node
                        obj.TFSF.fields1D.applyABC();        % apply ABC

                        % correct Ex along y-axis (start)
                        m = obj.TFSF.region(1,1):obj.TFSF.region(2,1)-1;
                        n = obj.TFSF.region(1,2);
                        obj.Ex(m,n) = obj.Ex(m,n) - obj.coefEx_z(m,n).*obj.TFSF.fields1D.Hy(m); % Hy <- Hz
                        % correct Ex along y-axis (end)
                        n = obj.TFSF.region(2,2);
                        obj.Ex(m,n) = obj.Ex(m,n) + obj.coefEx_z(m,n).*obj.TFSF.fields1D.Hy(m); % Hy <- Hz

                        % correct Ey field along x-axis (start)
                        m = obj.TFSF.region(1,1);
                        n = obj.TFSF.region(1,2):obj.TFSF.region(2,2)-1;
                        obj.Ey(m,n) = obj.Ey(m,n) + obj.coefEy_z(m,n).*obj.TFSF.fields1D.Hy(m-1); % Hy <- Hz
                        % correct Ey field along y-axis (end)
                        m = obj.TFSF.region(2,1);
                        obj.Ey(m,n) = obj.Ey(m,n) - obj.coefEy_z(m,n).*obj.TFSF.fields1D.Hy(m); % Hy <- Hz

                        % no need to correct Ex along top and bottom since
                        % incident Ex is zero

                    end
                case 3
                    % correct H adjacent to TFSF boundary
                    
                    %****** constant x faces (scattered-field nodes) ******
                    m_1 = obj.TFSF.region(1,1)-1;
                    m_2 = obj.TFSF.region(2,1);
                    n = obj.TFSF.region(1,2):obj.TFSF.region(2,2);
                    p = obj.TFSF.region(1,3):obj.TFSF.region(2,3)-1;
                    % correct Hy at firstX-1/2 by subtracting Ez_inc
                    obj.Hy(m_1,n,p) = obj.Hy(m_1,n,p) - obj.coefHy_x(m_1,n,p).*obj.TFSF.fields1D.Ez(m_1+1);
                    % correct Hy at lastX + 1/2 by adding Ez_inc
                    obj.Hy(m_2,n,p) = obj.Hy(m_2,n,p) + obj.coefHy_x(m_2,n,p).*obj.TFSF.fields1D.Ez(m_2);
                    
                    %**** constant y faces (scattered-field nodes) ****
                    m = obj.TFSF.region(1,1):obj.TFSF.region(2,1);
                    n_1 = obj.TFSF.region(1,2)-1;
                    n_2 = obj.TFSF.region(2,2);
                    for p = obj.TFSF.region(1,3):obj.TFSF.region(2,3)-1;
                        % correct Hx at firstY-1/2 by adding Ez_inc
                        obj.Hx(m,n_1,p) = obj.Hx(m,n_1,p) + obj.coefHx_y(m,n_1,p).*obj.TFSF.fields1D.Ez(m);
                        % correct Hx at lastY+1/2 by subtracting Ez_inc
                        obj.Hx(m,n_2,p) = obj.Hx(m,n_2,p) - obj.coefHx_y(m,n_2,p).*obj.TFSF.fields1D.Ez(m);
                    end
                    
                    %**** constant z faces (scattered-field nodes) ****
                    % nothing to correct on this face
                    
                    %**** update the fields in the auxiliary 1D grid ****
                    obj.TFSF.fields1D.updateH();         % update 1D magnetic field
                    obj.TFSF.fields1D.updateTFSF(time);	% in 1D, apply TFSF
                    obj.TFSF.fields1D.updateE();         % update 1D electric field
%                     obj.TFSF.fields1D.Ez(1) = SourceFunction(time,obj.TFSF.source,obj.CourantNumber); %<-Ey in 1D, set source node
                    obj.TFSF.fields1D.applyABC();        % apply ABC

                    % correct E adjacent to TFSF boundary
                    
                    %**** constant x faces (total-field nodes) ****
                    m_1 = obj.TFSF.region(1,1);
                    m_2 = obj.TFSF.region(2,1);
                    n = obj.TFSF.region(1,2):obj.TFSF.region(2,2);
                    p = obj.TFSF.region(1,3):obj.TFSF.region(2,3)-1;
                    % correct Ez at firstX face by subtracting Hy_inc
                    obj.Ez(m_1,n,p) = obj.Ez(m_1,n,p) - obj.coefEz_x(m_1,n,p).*obj.TFSF.fields1D.Hy(m_1-1);
                    % correct Ez at lastX face by adding Hy_inc
                    obj.Ez(m_2,n,p) = obj.Ez(m_2,n,p) + obj.coefEz_x(m_2,n,p).*obj.TFSF.fields1D.Hy(m_2);

                    %**** constant y faces (total-field nodes) ****
                    % nothing to correct on this face
                    
                    %**** constant z faces (total-field nodes) ****
                    m = obj.TFSF.region(1,1):obj.TFSF.region(2,1)-1;
                    p_1 = obj.TFSF.region(1,3);
                    p_2 = obj.TFSF.region(2,3);
                    for n = obj.TFSF.region(1,2):obj.TFSF.region(2,2);
                        % correct Ex at firstZ face by adding Hy_inc
                        obj.Ex(m,n,p_1) = obj.Ex(m,n,p_1) + obj.coefEx_y(m,n,p_1).*obj.TFSF.fields1D.Hy(m);
                        % correct Ex at lastZ face by subtracting Hy_inc
                        obj.Ex(m,n,p_2) = obj.Ex(m,n,p_2) - obj.coefEx_y(m,n,p_2).*obj.TFSF.fields1D.Hy(m);
                    end
            end
        end
        
        %% out field data
        function obj = initOUT(obj,Objects, MaxTime)
            if ~isempty(Objects)
                % check objects
                obj.OUT.Objects = Objects;
                for n_obj=length(obj.OUT.Objects):-1:1
                    ind = obj.OUT.Objects{n_obj};

                    tmp = diff(ind,[],2);

                    if any(tmp==0)
                        obj.OUT.data{n_obj}.X = obj.x(ind(1,1):ind(1,2),ind(2,1):ind(2,2),ind(3,1):ind(3,2));
                        obj.OUT.data{n_obj}.Y = obj.y(ind(1,1):ind(1,2),ind(2,1):ind(2,2),ind(3,1):ind(3,2));
                        obj.OUT.data{n_obj}.Z = obj.z(ind(1,1):ind(1,2),ind(2,1):ind(2,2),ind(3,1):ind(3,2));
                    end
                    
                    switch true
                        case tmp(1)==0
                            obj.OUT.data{n_obj}.type	= 1;
                            obj.OUT.data{n_obj}.Sx      = zeros([tmp(tmp~=0)+1;MaxTime]');
                        case tmp(2)==0
                            obj.OUT.data{n_obj}.type	= 2;
                            obj.OUT.data{n_obj}.Sy      = zeros([tmp(tmp~=0)+1;MaxTime]');
                        case tmp(3)==0
                            obj.OUT.data{n_obj}.type	= 3;
                            obj.OUT.data{n_obj}.Sz      = zeros([tmp(tmp~=0)+1;MaxTime]');
                        otherwise
                            obj.OUT.data{n_obj}.type	= 4;

                            obj.OUT.data{n_obj}.X1x = obj.x(ind(1,1),ind(2,1):ind(2,2),ind(3,1):ind(3,2));
                            obj.OUT.data{n_obj}.Y1x = obj.y(ind(1,1),ind(2,1):ind(2,2),ind(3,1):ind(3,2));
                            obj.OUT.data{n_obj}.Z1x = obj.z(ind(1,1),ind(2,1):ind(2,2),ind(3,1):ind(3,2));
                            obj.OUT.data{n_obj}.X2x = obj.x(ind(1,2),ind(2,1):ind(2,2),ind(3,1):ind(3,2));
                            obj.OUT.data{n_obj}.Y2x = obj.y(ind(1,2),ind(2,1):ind(2,2),ind(3,1):ind(3,2));
                            obj.OUT.data{n_obj}.Z2x = obj.z(ind(1,2),ind(2,1):ind(2,2),ind(3,1):ind(3,2));
                            obj.OUT.data{n_obj}.S1x(:,:,time)	= zeros([tmp(2:3)+1;MaxTime]');
                            obj.OUT.data{n_obj}.S2x(:,:,time)	= zeros([tmp(2:3)+1;MaxTime]');

                            obj.OUT.data{n_obj}.X1y = obj.x(ind(1,1):ind(1,2),ind(2,1),ind(3,1):ind(3,2));
                            obj.OUT.data{n_obj}.Y1y = obj.y(ind(1,1):ind(1,2),ind(2,1),ind(3,1):ind(3,2));
                            obj.OUT.data{n_obj}.Z1y = obj.z(ind(1,1):ind(1,2),ind(2,1),ind(3,1):ind(3,2));
                            obj.OUT.data{n_obj}.X2y = obj.x(ind(1,1):ind(1,2),ind(2,2),ind(3,1):ind(3,2));
                            obj.OUT.data{n_obj}.Y2y = obj.y(ind(1,1):ind(1,2),ind(2,2),ind(3,1):ind(3,2));
                            obj.OUT.data{n_obj}.Z2y = obj.z(ind(1,1):ind(1,2),ind(2,2),ind(3,1):ind(3,2));
                            obj.OUT.data{n_obj}.S1y(:,:,time)	= zeros([tmp([1,3])+1;MaxTime]');
                            obj.OUT.data{n_obj}.S2y(:,:,time)	= zeros([tmp([1,3])+1;MaxTime]');

                            obj.OUT.data{n_obj}.X1z = obj.x(ind(1,1):ind(1,2),ind(2,1):ind(2,2),ind(3,1));
                            obj.OUT.data{n_obj}.Y1z = obj.y(ind(1,1):ind(1,2),ind(2,1):ind(2,2),ind(3,1));
                            obj.OUT.data{n_obj}.Z1z = obj.z(ind(1,1):ind(1,2),ind(2,1):ind(2,2),ind(3,1));
                            obj.OUT.data{n_obj}.X2z = obj.x(ind(1,1):ind(1,2),ind(2,1):ind(2,2),ind(3,2));
                            obj.OUT.data{n_obj}.Y2z = obj.y(ind(1,1):ind(1,2),ind(2,1):ind(2,2),ind(3,2));
                            obj.OUT.data{n_obj}.Z2z = obj.z(ind(1,1):ind(1,2),ind(2,1):ind(2,2),ind(3,2));
                            obj.OUT.data{n_obj}.S1z(:,:,time)	= zeros([tmp(1:2)+1;MaxTime]');
                            obj.OUT.data{n_obj}.S2z(:,:,time)	= zeros([tmp(1:2)+1;MaxTime]');
                    end
                end
%{
                for n_obj=1:length(Objects)
                    otmp = Objects{n_obj};
                    switch otmp.type
                        case 'point'
                        case 'plane'
                            condition = otmp.normal(1)*(x-otmp.position(1))+otmp.normal(2)*(y-otmp.position(2))+otmp.normal(3)*(z-otmp.position(3))==0;
                            
                            Ex = @(X,Y,Z) interp3(obj.xEx,obj.yEx,obj.zEx,obj.Ex,X,Y,Z,'cubic');
                            Ey = @(X,Y,Z) interp3(obj.xEy,obj.yEy,obj.zEy,obj.Ey,X,Y,Z,'cubic');
                            Ez = @(X,Y,Z) interp3(obj.xEz,obj.yEz,obj.zEz,obj.Ez,X,Y,Z,'cubic');
                            Hx = @(X,Y,Z) interp3(obj.xHx,obj.YHx,obj.zHx,obj.Hx,X,Y,Z,'cubic');
                            Hy = @(X,Y,Z) interp3(obj.xHy,obj.yHy,obj.zHy,obj.Hy,X,Y,Z,'cubic');
                            Hz = @(X,Y,Z) interp3(obj.xHz,obj.yHz,obj.zHz,obj.Hz,X,Y,Z,'cubic');
                            

                    end
                end
%}                            
            end
        end

        function obj = updateOUT(obj,time)
            if isempty(obj.OUT)
                fprintf('Objects are not defined ...\n')
                return
            end
            % calculate Poynting
            for n_obj=length(obj.OUT.Objects):-1:1
                switch obj.OUT.data{n_obj}.type
                    case 1
                        Sx = Poynting(obj,obj.OUT.data{n_obj}.X,obj.OUT.data{n_obj}.Y,obj.OUT.data{n_obj}.Z);
                        % plane perpendicular to x-axis (Sx)
                        obj.OUT.data{n_obj}.Sx(:,:,time)	= squeeze(Sx);
                
                        Norm = (max(obj.OUT.data{n_obj}.Y(:))-min(obj.OUT.data{n_obj}.Y(:)))*(max(obj.OUT.data{n_obj}.Z(:))-min(obj.OUT.data{n_obj}.Z(:)));
                        X = obj.OUT.data{n_obj}.X(1);
                        obj.OUT.data{n_obj}.Sx_aver(time) = integral2(@(Y,Z)(interpn(obj.xEy,obj.yEy,obj.zEy,obj.Ey,X*ones(size(Y)),Y,Z,'cubic').*interpn(obj.xHz,obj.yHz,obj.zHz,obj.Hz,X*ones(size(Y)),Y,Z,'cubic') ...
                                                                            -interpn(obj.xEz,obj.yEz,obj.zEz,obj.Ez,X*ones(size(Y)),Y,Z,'cubic').*interpn(obj.xHy,obj.yHy,obj.zHy,obj.Hy,X*ones(size(Y)),Y,Z,'cubic')), ...
                            min(obj.OUT.data{n_obj}.Y(:)),max(obj.OUT.data{n_obj}.Y(:)), min(obj.OUT.data{n_obj}.Z(:)),max(obj.OUT.data{n_obj}.Z(:)))/Norm;
                    case 2
                       [~,Sy] = Poynting(obj,obj.OUT.data{n_obj}.X,obj.OUT.data{n_obj}.Y,obj.OUT.data{n_obj}.Z);
                        % plane perpendicular to y-axis (Sy)
                        obj.OUT.data{n_obj}.Sy(:,:,time)	= squeeze(Sy);
                        
                        Norm = (max(obj.OUT.data{n_obj}.X(:))-min(obj.OUT.data{n_obj}.X(:)))*(max(obj.OUT.data{n_obj}.Z(:))-min(obj.OUT.data{n_obj}.Z(:)));
                        Y = obj.OUT.data{n_obj}.Y(1);
                        obj.OUT.data{n_obj}.Sy_aver(time) = integral2(@(X,Z)(interpn(obj.xEz,obj.yEz,obj.zEz,obj.Ez,X,Y*ones(size(X)),Z,'cubic').*interpn(obj.xHx,obj.yHx,obj.zHx,obj.Hx,X,Y*ones(size(X)),Z,'cubic') ...
                                                                            -interpn(obj.xEx,obj.yEx,obj.zEx,obj.Ex,X,Y*ones(size(X)),Z,'cubic').*interpn(obj.xHz,obj.yHz,obj.zHz,obj.Hz,X,Y*ones(size(X)),Z,'cubic')), ...
                            min(obj.OUT.data{n_obj}.X(:)),max(obj.OUT.data{n_obj}.X(:)), min(obj.OUT.data{n_obj}.Z(:)),max(obj.OUT.data{n_obj}.Z(:)))/Norm;

                    case 3
                        [~,~,Sz] = Poynting(obj,obj.OUT.data{n_obj}.X,obj.OUT.data{n_obj}.Y,obj.OUT.data{n_obj}.Z);
                        % plane perpendicular to z-axis (Sz)
                        obj.OUT.data{n_obj}.Sz(:,:,time)	= Sz;

                        Norm = (max(obj.OUT.data{n_obj}.X(:))-min(obj.OUT.data{n_obj}.X(:)))*(max(obj.OUT.data{n_obj}.Y(:))-min(obj.OUT.data{n_obj}.Y(:)));
                        Z = obj.OUT.data{n_obj}.Z(1);
                        obj.OUT.data{n_obj}.Sz_aver(time) = integral2(@(X,Y)(interpn(obj.xEx,obj.yEx,obj.zEx,obj.Ex,X,Y,Z*ones(size(X)),'cubic').*interpn(obj.xHy,obj.yHy,obj.zHy,obj.Hy,X,Y,Z*ones(size(X)),'cubic') ...
                                                                            -interpn(obj.xEy,obj.yEy,obj.zEy,obj.Ey,X,Y,Z*ones(size(X)),'cubic').*interpn(obj.xHx,obj.yHx,obj.zHx,obj.Hx,X,Y,Z*ones(size(X)),'cubic')), ...
                            min(obj.OUT.data{n_obj}.X(:)),max(obj.OUT.data{n_obj}.X(:)), min(obj.OUT.data{n_obj}.Y(:)),max(obj.OUT.data{n_obj}.Y(:)))/Norm;
                    case 4
                        % x-planes
                        Sx = Poynting(obj,obj.OUT.data{n_obj}.X1x,obj.OUT.data{n_obj}.Y1x,obj.OUT.data{n_obj}.Z1x);
                        obj.OUT.data{n_obj}.S1x(:,:,time)	= squeeze(Sx);
                        Sx = Poynting(obj,obj.OUT.data{n_obj}.X2x,obj.OUT.data{n_obj}.Y2x,obj.OUT.data{n_obj}.Z2x);
                        obj.OUT.data{n_obj}.S2x(:,:,time)	= squeeze(Sx);

                        Norm = (max(obj.OUT.data{n_obj}.Y1x(:))-min(obj.OUT.data{n_obj}.Y1x(:)))*(max(obj.OUT.data{n_obj}.Z1x(:))-min(obj.OUT.data{n_obj}.Z1x(:)));
                        X = obj.OUT.data{n_obj}.X1x(1);
                        S_aver = integral2(@(Y,Z)(interpn(obj.xEy,obj.yEy,obj.zEy,obj.Ey,X*ones(size(Y)),Y,Z,'cubic').*interpn(obj.xHz,obj.yHz,obj.zHz,obj.Hz,X*ones(size(Y)),Y,Z,'cubic') ...
                                                 -interpn(obj.xEz,obj.yEz,obj.zEz,obj.Ez,X*ones(size(Y)),Y,Z,'cubic').*interpn(obj.xHy,obj.yHy,obj.zHy,obj.Hy,X*ones(size(Y)),Y,Z,'cubic')), ...
                            min(obj.OUT.data{n_obj}.Y1x(:)),max(obj.OUT.data{n_obj}.Y1x(:)), min(obj.OUT.data{n_obj}.Z1x(:)),max(obj.OUT.data{n_obj}.Z1x(:)))/Norm;
                        X = obj.OUT.data{n_obj}.X2x(1);
                        S_aver = S_aver -integral2(@(Y,Z)(interpn(obj.xEy,obj.yEy,obj.zEy,obj.Ey,X*ones(size(Y)),Y,Z,'cubic').*interpn(obj.xHz,obj.yHz,obj.zHz,obj.Hz,X*ones(size(Y)),Y,Z,'cubic') ...
                                                         -interpn(obj.xEz,obj.yEz,obj.zEz,obj.Ez,X*ones(size(Y)),Y,Z,'cubic').*interpn(obj.xHy,obj.yHy,obj.zHy,obj.Hy,X*ones(size(Y)),Y,Z,'cubic')), ...
                            min(obj.OUT.data{n_obj}.Y1x(:)),max(obj.OUT.data{n_obj}.Y1x(:)), min(obj.OUT.data{n_obj}.Z1x(:)),max(obj.OUT.data{n_obj}.Z1x(:)))/Norm;

                        % y-planes
                        [~,Sy] = Poynting(obj,obj.OUT.data{n_obj}.X1y,obj.OUT.data{n_obj}.Y1y,obj.OUT.data{n_obj}.Z1y);
                        obj.OUT.data{n_obj}.S1y(:,:,time)	= squeeze(Sy);
                        [~,Sy] = Poynting(obj,obj.OUT.data{n_obj}.X2y,obj.OUT.data{n_obj}.Y2y,obj.OUT.data{n_obj}.Z2y);
                        obj.OUT.data{n_obj}.S2y(:,:,time)	= squeeze(Sy);

                        Norm = (max(obj.OUT.data{n_obj}.X1y(:))-min(obj.OUT.data{n_obj}.X1y(:)))*(max(obj.OUT.data{n_obj}.Z1y(:))-min(obj.OUT.data{n_obj}.Z1y(:)));
                        Y = obj.OUT.data{n_obj}.Y1y(1);
                        S_aver = S_aver + integral2(@(X,Z)(interpn(obj.xEy,obj.yEy,obj.zEy,obj.Ey,X,Y*ones(size(X)),Z,'cubic').*interpn(obj.xHz,obj.yHz,obj.zHz,obj.Hz,X,Y*ones(size(X)),Z,'cubic') ...
                                                          -interpn(obj.xEz,obj.yEz,obj.zEz,obj.Ez,X,Y*ones(size(X)),Z,'cubic').*interpn(obj.xHy,obj.yHy,obj.zHy,obj.Hy,X,Y*ones(size(X)),Z,'cubic')), ...
                            min(obj.OUT.data{n_obj}.X1y(:)),max(obj.OUT.data{n_obj}.X1y(:)), min(obj.OUT.data{n_obj}.Z1y(:)),max(obj.OUT.data{n_obj}.Z1y(:)))/Norm;
                        Y = obj.OUT.data{n_obj}.Y2y(1);
                        S_aver = S_aver - integral2(@(X,Z)(interpn(obj.xEy,obj.yEy,obj.zEy,obj.Ey,X,Y*ones(size(X)),Z,'cubic').*interpn(obj.xHz,obj.yHz,obj.zHz,obj.Hz,X,Y*ones(size(X)),Z,'cubic') ...
                                                          -interpn(obj.xEz,obj.yEz,obj.zEz,obj.Ez,X,Y*ones(size(X)),Z,'cubic').*interpn(obj.xHy,obj.yHy,obj.zHy,obj.Hy,X,Y*ones(size(X)),Z,'cubic')), ...
                            min(obj.OUT.data{n_obj}.X1y(:)),max(obj.OUT.data{n_obj}.X1y(:)), min(obj.OUT.data{n_obj}.Z1y(:)),max(obj.OUT.data{n_obj}.Z1y(:)))/Norm;

                        % z-planes
                        [~,~,Sz] = Poynting(obj,obj.OUT.data{n_obj}.X1z,obj.OUT.data{n_obj}.Y1z,obj.OUT.data{n_obj}.Z1z);
                        obj.OUT.data{n_obj}.S1z(:,:,time)	= Sz;
                        [~,~,Sz] = Poynting(obj,obj.OUT.data{n_obj}.X2z,obj.OUT.data{n_obj}.Y2z,obj.OUT.data{n_obj}.Z2z);
                        obj.OUT.data{n_obj}.S2z(:,:,time)	= Sz;
                        
                        Norm = (max(obj.OUT.data{n_obj}.X1z(:))-min(obj.OUT.data{n_obj}.X1z(:)))*(max(obj.OUT.data{n_obj}.Y1z(:))-min(obj.OUT.data{n_obj}.Y1z(:)));
                        Z = obj.OUT.data{n_obj}.Z1z(1);
                        S_aver = S_aver + integral2(@(X,Y)(interpn(obj.xEy,obj.yEy,obj.zEy,obj.Ey,X,Y,Z*ones(size(X)),'cubic').*interpn(obj.xHz,obj.yHz,obj.zHz,obj.Hz,X,Y,Z*ones(size(X)),'cubic') ...
                                                          -interpn(obj.xEz,obj.yEz,obj.zEz,obj.Ez,X,Y,Z*ones(size(X)),'cubic').*interpn(obj.xHy,obj.yHy,obj.zHy,obj.Hy,X,Y,Z*ones(size(X)),'cubic')), ...
                            min(obj.OUT.data{n_obj}.X1z(:)),max(obj.OUT.data{n_obj}.X1z(:)), min(obj.OUT.data{n_obj}.Y1z(:)),max(obj.OUT.data{n_obj}.Y1z(:)))/Norm;
                        Z = obj.OUT.data{n_obj}.Z2z(1);
                        S_aver = S_aver - integral2(@(X,Y)(interpn(obj.xEy,obj.yEy,obj.zEy,obj.Ey,X,Y,Z*ones(size(X)),'cubic').*interpn(obj.xHz,obj.yHz,obj.zHz,obj.Hz,X,Y,Z*ones(size(X)),'cubic') ...
                                                          -interpn(obj.xEz,obj.yEz,obj.zEz,obj.Ez,X,Y,Z*ones(size(X)),'cubic').*interpn(obj.xHy,obj.yHy,obj.zHy,obj.Hy,X,Y,Z*ones(size(X)),'cubic')), ...
                            min(obj.OUT.data{n_obj}.X1z(:)),max(obj.OUT.data{n_obj}.X1z(:)), min(obj.OUT.data{n_obj}.Y1z(:)),max(obj.OUT.data{n_obj}.Y1z(:)))/Norm;

                        obj.OUT.data{n_obj}.S_aver(time) = S_aver;
                end
            end
        end
        
        %% calculate source function at given time and location
        function obj = add_Ez_inc(obj,Position,time,sourceStr)
%             if (sourceStr.ppw <= 0)
%                 fprintf(stderr,'E_z(inc): ezIncInit() must be called before ezInc.\n Points per wavelength must be positive.\n');
%                 exit(-1);
%             end
            
            source = SourceFunction(time,sourceStr,obj.CourantNumber(1));
            switch obj.Dimensionality
                case 1
                    if sourceStr.isAdd
                        obj.Ez(Position(1)) = obj.Ez(Position(1)) + source;
                    else
                        obj.Ez(Position(1)) = source;
                    end
                    fprintf('%d, %g\n',time,obj.Ez(Position(1)))
                case 2
                    if sourceStr.isAdd
                        obj.Ez(Position(1),Position(2)) = obj.Ez(Position(1),Position(2)) + source(1);
                    else
                        obj.Ez(Position(1),Position(2)) = source(1);
                    end
                    fprintf('%d, %g\n',time,obj.Ez(Position(1),Position(2)))
                case 3
                    if sourceStr.isAdd
                        obj.Ez(Position(1),Position(2),Position(3)) = obj.Ez(Position(1),Position(2),Position(3)) + source(1);
                    else
                        obj.Ez(Position(1),Position(2),Position(3)) = source(1);
                    end
                    fprintf('%d, %g\n',time,obj.Ez(Position(1),Position(2),Position(3)))
            end
        end
        
        function obj = add_Hy_inc(obj,Position,time,sourceStr)
            global imp0
            source = SourceFunction(time,sourceStr,obj.CourantNumber(1)) / imp0;
            switch obj.Dimensionality
                case 1
                    if sourceStr.isAdd
                        obj.Hy(Position(1)) = obj.Hy(Position(1)) - source;
                    else
                        obj.Hy(Position(1)) = -source;
                    end
                    fprintf('%d, %g\n',time,obj.Hy(Position(1)))
                case 2
                    if sourceStr.isAdd
                        obj.Hy(Position(1),Position(2)) = obj.Hy(Position(1),Position(2)) - source;
                    else
                        obj.Hy(Position(1),Position(2)) = -source;
                    end
                    fprintf('%d, %g\n',time,obj.Hy(Position(1),Position(2)))
                case 3
                    if sourceStr.isAdd
                        obj.Hy(Position(1),Position(2),Position(3)) = obj.Hy(Position(1),Position(2),Position(3)) - source;
                    else
                        obj.Hy(Position(1),Position(2),Position(3)) = -source;
                    end
                    fprintf('%d, %g\n',time,obj.Hy(Position(1),Position(2),Position(3)))
            end
        end
        
        %% field components outputs
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
        function out = snapshotH(obj,type)
            if ~exist('type','var')
                type = '';
            end
            if strcmpi(type,'intensity')
                switch obj.Dimensionality
                    case 1
                        out = obj.Hy.^2;
                    case 2
                        if obj.isTE
                            out = obj.snapshotH('x0').^2 + obj.snapshotH('y0').^2;
                        else
                            out = obj.Hz.^2;
                        end
                    case 3
                        out = obj.snapshotH('x0').^2 + obj.snapshotH('y0').^2 + obj.snapshotH('z0').^2;
                end
                return
            end
            
            switch obj.Dimensionality
                case 1
                    out = obj.Hy;
                case 2
                    if obj.isTE
                        switch type
                            % return as it
                            case 'x'
                                out = obj.Hx;
                            case 'y'
                                out = obj.Hy;
                            % fit to space grid
                            case 'x0'
                                out = (obj.Hx(1:end-1,:)+obj.Hx(2:end,:))/2;
                            case 'y0'
                                out = (obj.Hy(:,1:end-1)+obj.Hy(:,2:end))/2;
                            otherwise
                                out = obj.Hx;
                        end
                    else
                        out = obj.Hz;
                    end
                case 3
                    switch type
                        % return as it
                        case 'x'
                            out = obj.Hx;
                        case 'y'
                            out = obj.Hy;
                        case 'z'
                            out = obj.Hz;
                        % fit to space grid
                        case 'x0'
                            out = (obj.Hx(1:end-1,:,:)+obj.Hx(2:end,:,:))/2;
                        case 'y0'
                            out = (obj.Hy(:,1:end-1,:)+obj.Hy(:,2:end,:))/2;
                        case 'z0'
                            out = (obj.Hz(:,:,1:end-1)+obj.Hz(:,:,2:end))/2;
                        otherwise
                            out = obj.Ez;
                    end
            end
        end
        
        function [Sx,Sy,Sz] = Poynting(obj,X,Y,Z)
            if nargin==1
                Sx = real(obj.snapshotE('y0').*obj.snapshotH('z0')-obj.snapshotE('z0').*obj.snapshotH('y0'));
                Sy = real(obj.snapshotE('z0').*obj.snapshotH('x0')-obj.snapshotE('x0').*obj.snapshotH('z0'));
                Sz = real(obj.snapshotE('x0').*obj.snapshotH('y0')-obj.snapshotE('y0').*obj.snapshotH('x0'));
            else
                tEx = interpn(obj.xEx,obj.yEx,obj.zEx,obj.Ex,X,Y,Z,'cubic');
                tEy = interpn(obj.xEy,obj.yEy,obj.zEy,obj.Ey,X,Y,Z,'cubic');
                tEz = interpn(obj.xEz,obj.yEz,obj.zEz,obj.Ez,X,Y,Z,'cubic');
                tHx = interpn(obj.xHx,obj.yHx,obj.zHx,obj.Hx,X,Y,Z,'cubic');
                tHy = interpn(obj.xHy,obj.yHy,obj.zHy,obj.Hy,X,Y,Z,'cubic');
                tHz = interpn(obj.xHz,obj.yHz,obj.zHz,obj.Hz,X,Y,Z,'cubic');
                Sx = real(tEy.*tHz-tEz.*tHy);
                Sy = real(tEz.*tHx-tEx.*tHz);
                Sz = real(tEx.*tHy-tEy.*tHx);

                
            end
        end

        %% get grid coordinates
        function [X,Y,Z] = getGridXYZ(obj,type)
            if ~exist('type','var')
                type = '';
            end
            
            switch obj.Dimensionality
                case 1
                    switch type
                        case 'Ez'
                            X = linspace(obj.FDTDspace(1,1),obj.FDTDspace(1,2),obj.Nx+1);
                        case 'Hy'
                            X = linspace(obj.FDTDspace(1,1)+obj.dx/2,obj.FDTDspace(1,2)-obj.dx/2,obj.Nx);
                        otherwise
                            % averaged grid
                            X = linspace(obj.FDTDspace(1,1)+obj.dx/2,obj.FDTDspace(1,2)-obj.dx/2,obj.Nx);
                    end
                case 2
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
                case 3
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
        end
        
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
    end

    methods(Access = protected)
        function coef = coefABC(obj,coefEz_y,coefHy_z,order)
            % tmp = Courant/sqrt(mu_r*eps_r)
            switch obj.Dimensionality
                case 1
                    tmp1 = sqrt(mean(obj.(coefEz_y)(1:2))*obj.(coefHy_z)(1));
                    tmp2 = sqrt(mean(obj.(coefEz_y)(end-(0:1)))*obj.(coefHy_z)(end));
                case 2
                    tmp1 = sqrt(mean(obj.(coefEz_y)(1:2,1))*obj.(coefHy_z)(1,1));
                    tmp2 = sqrt(mean(obj.(coefEz_y)(end-(0:1),1))*obj.(coefHy_z)(end,1));
                case 3
                    tmp1 = sqrt(mean(obj.(coefEz_y)(1:2,1))*obj.(coefHy_z)(1,1));
                    tmp2 = sqrt(mean(obj.(coefEz_y)(end-(0:1),1))*obj.(coefHy_z)(end,1));
%                     tmp1 = obj.CourantNumber;
%                     tmp2 = obj.CourantNumber;
            end
            if order==1
                coef = [(tmp1-1)./(tmp1+1),(tmp2-1)./(tmp2+1)];
            else
                coef = [[-(1/tmp1-2+tmp1);-2*(tmp1-1/tmp1);4*(tmp1+1/tmp1)]/(1/tmp1+2+tmp1), ...
                        [-(1/tmp2-2+tmp2);-2*(tmp2-1/tmp2);4*(tmp2+1/tmp2)]/(1/tmp2+2+tmp2)];
            end
        end
        
    end
end

%%
function coef = coefABC2(S,order)
    if order==1
        coef(1,2) = (S-1)/(S+1);
        coef(1,1) = coef(2,1); % (start:end)
    else
        coef(1:3,2) = [-(1./S-2+S);
                         -2*(S-1./S);
                          4*(S+1./S)]./(1./S+2+S);
        coef(1:3,1) = coef(1:3,2); % (start:end),(0,+/-1,+/-2)
    end
end

function source = SourceFunction(time,sourceStr,CourantNumber)
    switch sourceStr.Type
    case 'gauss'
        source = exp(-(CourantNumber/sourceStr.ppw*(time-sourceStr.timeDelay)).^2);
    case 'Ricker'
        % The Ricker Wavelet
        arg = (pi*((CourantNumber*time - sourceStr.timeDelay)/sourceStr.ppw - 1)).^2;
        source = (1 - 2*arg).*exp(-arg);%.*sin(sourceStr.omega*time); % omega = (omega*dt)
    end
%     source = exp(-((time - delay - location / cdtds) / width)^2);
%     source = sin(2*pi/ppw * (cdtds * time - location));
%  arg = (pi * ((cdtds * time - location) / ppw - 1.0))^2;
%     source = (1 - 2*arg)*exp(-arg);
end
