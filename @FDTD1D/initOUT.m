%% initialize out field data
function obj = initOUT(obj,Objects, MaxTime)
%#codegen

if ~isempty(Objects)
    % check objects
    obj.OUT.Objects = Objects;
    for n_obj=length(obj.OUT.Objects):-1:1
        if obj.PMLtype>0
            ind = obj.OUT.Objects{n_obj}(1,:) + obj.PML.layers;
        else
            ind = obj.OUT.Objects{n_obj}(1,:);
        end

        tmp = diff(ind(1,:),[],2);

        obj.OUT.data{n_obj}.X = obj.x(ind(1,1):ind(1,2));
                
        obj.OUT.data{n_obj}.type	= 1;
        obj.OUT.data{n_obj}.Sx      = zeros([tmp(1)+1;MaxTime]');
        obj.OUT.data{n_obj}.Sx_aver	= obj.OUT.data{n_obj}.Sx;
        obj.OUT.data{n_obj}.Hy      = obj.OUT.data{n_obj}.Sx;
        obj.OUT.data{n_obj}.Ez      = obj.OUT.data{n_obj}.Sx;
    end
end

end

