%% initialize TFSF boundary
 function obj = initTFSF(obj,TFSFregion,TFSFsource)
    obj.TFSF.region	= TFSFregion;
    obj.TFSF.source	= TFSFsource;
    obj.TFSF.source.omega = TFSFsource.omega*obj.dt;
    switch obj.Dimensionality
        case {2,3}
            % for TE Ey -> -Ez, Hz -> Hy
            obj.TFSF.fields1D = FDTD(obj.Nx+1,obj.FDTDspace(1,:),obj.PML);	% create the 1D FDTD grid
            obj.TFSF.fields1D.initGrid([],obj.CourantNumber);               % initialize the 1D FDTD grid with Courant number for 2d/3D model
            obj.TFSF.fields1D.initTFSF(1,obj.TFSF.source);
            if isempty(obj.PML)
                obj.TFSF.fields1D.initABC('2abc');                          % initialize ABC
			else
				obj.TFSF.region	= TFSFregion + obj.PML.layers;
            end
    end
end
