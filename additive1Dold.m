% 1D FDTD simulation with an additive source.

SIZE = 200;
ez = zeros(1,SIZE);
hy = zeros(1,SIZE);
imp0 = 377.0; % impedance
maxTime = 400;

basename = "sim";
additiveSource = 50;


frame = 0;

%FILE *snapshot;


figure
hold on
view(3)

% do time stepping
mmH = 1:SIZE-1;
mmE = 2:SIZE;
for qTime =0:maxTime-1
  % update magnetic field
  hy(1,mmH) += (ez(1,mmH+1) - ez(1,mmH)) / imp0;
  
  % update electric field
  ez(1,mmE) += (hy(1,mmE) - hy(1,mmE-1)) * imp0;
  
  % use additive source at node $additiveSource
  ez(additiveSource) += exp(-(qTime-30)^2 / 100);
  
  % show/write snapshot if time a multiple of 10
  if (mod(qTime,10) == 0)
  %{
    sprintf(filename, "%s.%d", basename, frame++);
    snapshot=fopen(filename, "w");
    for mm=1:SIZE
      fprintf(snapshot, "%g\n", ez[mm]);
      end
    fclose(snapshot);
    %}
    waterfall(1:SIZE,qTime,ez);
  end
end %end of time-stepping