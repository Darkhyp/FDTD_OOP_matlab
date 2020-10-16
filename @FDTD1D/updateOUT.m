%% Store field data
function obj = updateOUT(obj,n_time)
%#codegen

if isempty(obj.OUT)
    fprintf('Objects are not defined ...\n')
    return
end

% calculate Poynting
for n_obj=length(obj.OUT.Objects):-1:1
    switch obj.OUT.data{n_obj}.type
        case 1
            [Sx,Hy,Ez] = Poynting(obj,obj.OUT.data{n_obj}.X);
            % plane perpendicular to x-axis (Sx)
            obj.OUT.data{n_obj}.Hy(n_time)      = Hy;
            obj.OUT.data{n_obj}.Ez(n_time)      = Ez;
            obj.OUT.data{n_obj}.Sx(n_time)      = Sx;
            obj.OUT.data{n_obj}.Sx_aver(n_time)	= Sx;
    end
end

end
