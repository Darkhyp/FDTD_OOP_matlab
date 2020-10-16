%% initialize TFSF boundary
 function obj = initTFSF(obj,TFSFregion,TFSFsource)
%#codegen

obj.TFSF.region	= TFSFregion(1,1);
obj.TFSF.source	= TFSFsource;
if isfield(TFSFsource,'omega')
    obj.TFSF.source.omega = TFSFsource.omega*obj.dt;
end


end
