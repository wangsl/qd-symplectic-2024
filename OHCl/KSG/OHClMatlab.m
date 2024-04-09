
function [] = OHClMatlab()

global OHClData

fprintf(' From OHClMatlab\n')


if mod(OHClData.time.steps, OHClData.options.steps_to_copy_psi_from_device_to_host) == 0
  PlotCRP();
  if OHClData.options.calculate_reaction_probabilities == 1
    fprintf(' Save reaction probabilities to %s\n', OHClData.CRP.mat_file);
    CRP = OHClData.CRP;
    save(OHClData.CRP.mat_file, 'CRP');
  end
end

if mod(OHClData.time.steps, OHClData.options.steps_to_copy_psi_from_device_to_host) == 0
  PlotPotWave()
end

return
