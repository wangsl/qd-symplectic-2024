
function [] = PlotCRP()

global H2eV 
global HO2Data

E = HO2Data.CRP.energies*H2eV;
CRP = -HO2Data.CRP.CRP;

persistent CRPPlot
persistent fig_index

if isempty(CRPPlot) 
  figure(2)
  CRPPlot = plot(E, CRP, 'b-', 'LineWidth', 2, 'YDataSource', ...
		 'CRP');
  
  set(gca, 'xtick', [0.8:0.2:1.8]);
  
  xlabel('Energy (eV)')
  ylabel('Reaction probabilities')
  
  axis([0.8 1.8 0 0.70])
  
  grid on
  grid minor
end

refreshdata(CRPPlot, 'caller');

drawnow

if isempty(fig_index) 
  fig_index = 0;
else
  fig_index = fig_index + 1;
end

index = strcat('000000000', int2str(fig_index));
fig_name = strcat('crp/', index(end-8:end), '.png');

fprintf(' Figure: %s\n', fig_name)

figure(2);
print(fig_name, '-dpng');


