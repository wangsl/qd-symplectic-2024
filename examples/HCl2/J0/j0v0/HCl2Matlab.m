
function [] = HCl2Matlab()

global HCl2Data

fprintf(' From HCl2Matlab\n')

if mod(HCl2Data.time.steps, HCl2Data.options.steps_to_copy_psi_from_device_to_host) == 0
  PlotCRP();
  if HCl2Data.options.calculate_reaction_probabilities == 1
    fprintf(' Save reaction probabilities to %s\n', HCl2Data.CRP.mat_file);
    CRP = HCl2Data.CRP;
    save(HCl2Data.CRP.mat_file, 'CRP');
  end
end

if mod(HCl2Data.time.steps, HCl2Data.options.steps_to_copy_psi_from_device_to_host) == 0
  PlotPotWave()
end

return

