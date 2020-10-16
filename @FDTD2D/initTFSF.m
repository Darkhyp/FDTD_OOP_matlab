%% initialize TFSF boundary
 function obj = initTFSF(obj,TFSFregion,TFSFsource)
%#codegen

obj.TFSF.region	= TFSFregion;
obj.TFSF.source	= TFSFsource;
if isfield(TFSFsource,'omega')
    obj.TFSF.source.omega = TFSFsource.omega*obj.dt;
end


% for TE Ey -> -Ez, Hz -> Hy
obj.TFSF.fields1D = FDTD1D(obj.N0x,obj.FDTDspace(1,:),obj.PML);	% create the 1D FDTD grid
obj.TFSF.fields1D.initGrid([],obj.Courant);                     % initialize the 1D FDTD grid with Courant number for 2d/3D model

if isempty(obj.PML)
    obj.TFSF.fields1D.initTFSF(1,TFSFsource);
	obj.TFSF.fields1D.initABC('2abc');                          % initialize ABC
else
    obj.TFSF.fields1D.initTFSF(1+obj.PML.layers,TFSFsource);
	obj.TFSF.region	= TFSFregion(1:2,:) + obj.PML.layers;
end

end
