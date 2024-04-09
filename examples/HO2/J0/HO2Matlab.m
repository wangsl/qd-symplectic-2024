
function [] = HO2Matlab()

global HO2Data

fprintf(' From HO2Matlab\n')

if mod(HO2Data.time.steps, HO2Data.options.steps_to_copy_psi_from_device_to_host) == 0
  PlotCRP();
  if HO2Data.options.calculate_reaction_probabilities == 1
    fprintf(' Save reaction probabilities to %s\n', HO2Data.CRP.mat_file);
    CRP = HO2Data.CRP;
    save(HO2Data.CRP.mat_file, 'CRP');
  end
end

if mod(HO2Data.time.steps, HO2Data.options.steps_to_copy_psi_from_device_to_host) == 0
  PlotPotWave()
end

return

