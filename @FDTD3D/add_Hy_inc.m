function obj = add_Hy_inc(obj,Position,time,sourceStr)
%#codegen

global imp0
source = SourceFunction(time,sourceStr,obj.CourantNumber(1)) / imp0;

if sourceStr.isAdd
	obj.Hy(Position(1),Position(2),Position(3)) = obj.Hy(Position(1),Position(2),Position(3)) - source;
else
	obj.Hy(Position(1),Position(2),Position(3)) = -source;
end
fprintf('%d, %g\n',time,obj.Hy(Position(1),Position(2),Position(3)))

end
