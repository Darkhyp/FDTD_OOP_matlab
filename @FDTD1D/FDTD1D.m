classdef FDTD1D < handle
%#codegen

    properties (Access = public)
    
        % electro-magnetic field
        Ez, Hy

        % grid coordinates
        x
        xEz
        xHy
        
        % structure for out-field 
        OUT
        
    end
    properties (Access = private)
        
        % Perfect-matched layers
        PMLtype	= 0; % no PML
        PML
        
        StabilityFactor = 0.99;
        Courant

        Nx
        N0x
        FDTDspace
        dx
        dt
        
        coefHy,coefHy_z % H-field coefficients
        coefEz,coefEz_y % E-field coefficients
        
        mu_yy
        epsilon_zz
        sigma_zz
        
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
        
        function obj = FDTD1D(DomainSize, Space, PML)
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
            obj.Nx = obj.N0x;

			obj.dx = diff(Space(1,:))/obj.N0x;

			if obj.PMLtype>0
				obj.Nx = obj.Nx + 2*obj.PML.layers;

				obj.PML.Wx      = obj.PML.layers*obj.dx;
				obj.PML.x1      = Space(1,1);
				obj.PML.x2      = Space(1,2);

				obj.FDTDspace = Space(1,:) + [-obj.PML.Wx,obj.PML.Wx];

                switch obj.PMLtype
					case 1
						% Ez grid (2:end-1,2:end-1,:)
						obj.PML.QyHx	= zeros(obj.Nx-1,1);
						obj.PML.QxHy	= zeros(obj.Nx-1,1);

						% Hy grid
						obj.PML.QxEz	= zeros(obj.Nx,1);
						obj.PML.QzEx	= zeros(obj.Nx,1);
					case 2
                        obj.PML.My      = zeros( obj.Nx,   1);

                        obj.PML.Fz      = zeros( obj.Nx-1, 1);
                        obj.PML.Gz      = zeros( obj.Nx-1, 1);
                end
			end                    
		
			obj.Hy          = zeros(obj.Nx, 1);
			
			obj.Ez          = zeros(obj.Nx+1, 1);
            
            obj.mu_yy       = ones( obj.Nx,  1);

            obj.epsilon_zz	= ones( obj.Nx-1, 1);
            obj.sigma_zz	= zeros(obj.Nx-1, 1);

            obj.coefHy      = ones( obj.Nx, 1);
            obj.coefHy_z	= ones( obj.Nx, 1);

            obj.coefEz      = ones( obj.Nx-1, 1);
            obj.coefEz_y	= ones( obj.Nx-1, 1);

%           obj.Ey = @(obj)  obj.Ez;
%           obj.Hz = @(obj) -obj.Hy;

            obj.dt = obj.StabilityFactor/(sqrt(1/obj.dx^2)*c);
            obj.Courant = obj.dt*c./obj.dx; % Courant number
%             obj.Courant = 1; % Courant number

        end
        
    end
    
    methods (Access = public)

        % init grid
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
        
        [Sx,Ez,Hy] = Poynting(obj,X)

        %% get grid coordinates
        X = getGridXYZ(obj,type)

    end

    methods(Access = protected)
        
        coef = coefABC(obj,coefEz_y,coefHy_z,order)
        
    end
end




