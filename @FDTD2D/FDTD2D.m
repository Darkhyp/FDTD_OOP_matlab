classdef FDTD2D < handle
%#codegen

    properties (Access = public)
    
        % electro-magnetic field
        Ez, Hx,Hy
        Ex,Ey, Hz

        % grid coordinates
        x, y
        xEx,yEx, xEy,yEy, xEz,yEz
        xHx,yHx, xHy,yHy, xHz,yHz
        
        % structure for out-field 
        OUT
        
    end
    properties (Access = private)
        
        isTE        = true;
        
        % Perfect-matched layers
        PMLtype	= 0; % no PML
        PML
        
        StabilityFactor = 0.99;
        Courant

        Nx, Ny
        N0x, N0y
        FDTDspace
        dx,dy, dt
        
        coefHx,coefHx_y,coefHx_z, coefHy,coefHy_x,coefHy_z, coefHz,coefHz_x,coefHz_y % H-field coefficients
        coefEx,coefEx_y,coefEx_z, coefEy,coefEy_x,coefEy_z, coefEz,coefEz_x,coefEz_y % E-field coefficients
        
        mu_xx, mu_yy, mu_zz
        epsilon_xx, epsilon_yy, epsilon_zz
        sigma_xx, sigma_yy, sigma_zz

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
%         function set.Courant(obj,Courant)
%             obj.Courant = Courant;
%         end
        function Courant = get_Courant(obj)
            Courant = obj.Courant;
        end
        
        function obj = FDTD2D(DomainSize, Space, PML, CourantNumber, isTE)
            global c
            
            obj.StabilityFactor = 0.99;
            
            obj.FDTDspace = Space;
            
            if ~isempty(PML)
                obj.PML     = PML;
				switch PML.type
					case 'split PML'
						obj.PMLtype	= 1;
					case 'UPML' % unsplit PML 
						obj.PMLtype	= 2;
				end
            end
            
			obj.N0x = DomainSize(1);
			obj.N0y = DomainSize(2);
            obj.Nx = obj.N0x;
            obj.Ny = obj.N0y;

			obj.dx = diff(Space(1,:))/obj.N0x;
			obj.dy = diff(Space(2,:))/obj.N0y;

			if nargin==4 % check polarization
				obj.isTE = isTE;
			end
			
			if obj.PMLtype>0
				obj.Nx = obj.N0x + 2*obj.PML.layers;
				obj.Ny = obj.N0y + 2*obj.PML.layers;
				
				obj.PML.Wx      = obj.PML.layers*obj.dx;
				obj.PML.x1      = Space(1,1);
				obj.PML.x2      = Space(1,2);

				obj.PML.Wy      = obj.PML.layers*obj.dy;
				obj.PML.y1      = Space(2,1);
				obj.PML.y2      = Space(2,2);

				obj.FDTDspace = Space(1:2,:) + [-obj.PML.Wx,obj.PML.Wx; -obj.PML.Wy,obj.PML.Wy];
			end                    
			
			if obj.isTE
				obj.Hx          = zeros(obj.Nx+1, obj.Ny);
				obj.Hy          = zeros(obj.Nx,   obj.Ny+1);
				obj.Ez          = zeros(obj.Nx+1, obj.Ny+1);

                switch obj.PMLtype
					case 1
						% Ez grid (2:end-1,2:end-1,:)
						obj.PML.QyHx	= zeros(obj.Nx-1, obj.Ny-1);
						obj.PML.QxHy	= zeros(obj.Nx-1, obj.Ny-1);

						% Hx grid
						obj.PML.QyEz	= zeros(obj.Nx+1, obj.Ny);
						obj.PML.QzEy	= zeros(obj.Nx+1, obj.Ny);
						% Hy grid
						obj.PML.QxEz	= zeros(obj.Nx, obj.Ny+1);
						obj.PML.QzEx	= zeros(obj.Nx, obj.Ny+1);
					case 2
                        obj.PML.Mx      = zeros( obj.Nx+1, obj.Ny);
                        obj.PML.My      = zeros( obj.Nx,   obj.Ny+1);

                        obj.PML.Fz      = zeros(obj.Nx-1, obj.Ny-1);
                        obj.PML.Gz      = zeros(obj.Nx-1, obj.Ny-1);
                end                    
                obj.mu_xx       = ones( obj.Nx+1, obj.Ny);
                obj.mu_yy       = ones( obj.Nx,   obj.Ny+1);
                
                obj.epsilon_zz	= ones( obj.Nx-1, obj.Ny-1);
                obj.sigma_zz	= zeros(obj.Nx-1, obj.Ny-1);

                obj.coefHx		= ones( obj.Nx+1, obj.Ny);
                obj.coefHx_z	= ones( obj.Nx+1, obj.Ny);

                obj.coefHy		= ones( obj.Nx, obj.Ny+1);
                obj.coefHy_z	= ones( obj.Nx, obj.Ny+1);

                obj.coefEz		= ones( obj.Nx-1, obj.Ny-1);
                obj.coefEz_x	= ones( obj.Nx-1, obj.Ny-1);
                obj.coefEz_y	= ones( obj.Nx-1, obj.Ny-1);
            else
                obj.Hz          = zeros(obj.Nx,   obj.Ny);

                obj.Ex          = zeros(obj.Nx,   obj.Ny+1);
                obj.Ey          = zeros(obj.Nx+1, obj.Ny);

                switch obj.PMLtype
					case 1
						% Ex grid (:,2:end-1)
						obj.PML.QyHz	= zeros(obj.Nx, obj.Ny-1);
						obj.PML.QzHy	= zeros(obj.Nx, obj.Ny-1);
						% Ey grid (2:end-1,:,2:end-1)
						obj.PML.QxHz	= zeros(obj.Nx-1, obj.Ny);
						obj.PML.QzHx	= zeros(obj.Nx-1, obj.Ny);
						
						% Hz grid
						obj.PML.QxEy	= zeros(obj.Nx, obj.Ny);
						obj.PML.QyEx	= zeros(obj.Nx, obj.Ny);
					case 2
                        obj.PML.Mz      = zeros( obj.Nx,   obj.Ny);

                        obj.PML.Fx      = zeros( obj.Nx, obj.Ny-1);
                        obj.PML.Gx      = zeros( obj.Nx, obj.Ny-1);

                        obj.PML.Fy      = zeros( obj.Nx-1, obj.Ny);
                        obj.PML.Gy      = zeros( obj.Nx-1, obj.Ny);
                end                    
                obj.mu_zz       = ones( obj.Nx,   obj.Ny);

                obj.epsilon_xx	= ones( obj.Nx,   obj.Ny-1);
                obj.sigma_xx	= zeros(obj.Nx,   obj.Ny-1);

                obj.epsilon_yy	= ones( obj.Nx-1, obj.Ny);
                obj.sigma_yy	= zeros(obj.Nx-1, obj.Ny);

                obj.coefHz		= ones( obj.Nx, obj.Ny);
                obj.coefHz_x	= ones( obj.Nx, obj.Ny);
                obj.coefHz_y	= ones( obj.Nx, obj.Ny);

                obj.coefEx		= ones( obj.Nx, obj.Ny-1);
                obj.coefEx_z	= ones( obj.Nx, obj.Ny-1);

                obj.coefEy		= ones( obj.Nx-1, obj.Ny);
                obj.coefEy_z	= ones( obj.Nx-1, obj.Ny);
			end

            obj.dt = obj.StabilityFactor/(sqrt(1/obj.dx^2+1/obj.dy^2)*c);
            obj.Courant = obj.dt*c./[obj.dx,obj.dy]; % Courant number
%                     obj.Courant = 1/sqrt(2); % Courant number

        end
        
    end
    
    methods (Access = public)

        % update magnetic field
        obj = initGrid(obj, Objects, Courant)
        
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
        
        [Sx,Sy,Ex,Ey,Ez,Hx,Hy,Hz] = Poynting(obj,X,Y)

        %% get grid coordinates
        [X,Y] = getGridXYZ(obj,type)

    end

    methods(Access = protected)
        
        coef = coefABC(obj,coefEz_y,coefHy_z,order)
        
    end
end




