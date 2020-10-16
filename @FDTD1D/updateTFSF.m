%% Apply TFSF
function obj = updateTFSF(obj,time)
%#codegen
if isempty(obj.TFSF) % check if initTFSF() has been called
	fprintf('\nUse initTFSF() method before updateTFSF()!...')
%                 exit(-1);
end
global imp0

if isempty(obj.PML)
    % correct Hy adjacent to TFSF boundary
    obj.Hy(obj.TFSF.region) = obj.Hy(obj.TFSF.region) ...
        - obj.coefHy_z(obj.TFSF.region)*SourceFunction(time, obj.TFSF.source,obj.Courant(1));
    %{
    figure
    T = 1:500;
    plot(T*obj.dt/1e15,SourceFunction(T, obj.TFSF.source,obj.Courant(1)),'.-')
        
    %}
    
    % correct Ez adjacent to TFSF boundary (time+0.5-(-0.5))
    obj.Ez(obj.TFSF.region+1) = obj.Ez(obj.TFSF.region+1) ...
        + SourceFunction(time+1,obj.TFSF.source,obj.Courant(1));
else
    % correct Hy adjacent to TFSF boundary
    obj.Hy(obj.TFSF.region) = obj.Hy(obj.TFSF.region) ...
        - SourceFunction(time, obj.TFSF.source,obj.Courant(1))/imp0;
    % correct Ez adjacent to TFSF boundary (time+0.5-(-0.5))
    obj.Ez(obj.TFSF.region+1) = obj.Ez(obj.TFSF.region+1) ...
        + SourceFunction(time+1,obj.TFSF.source,obj.Courant(1));
end

end
