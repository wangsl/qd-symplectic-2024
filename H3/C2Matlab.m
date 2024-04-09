
% $Id$

function [] = C2Matlab() 

global H3Data

if mod(H3Data.time.steps, H3Data.options.steps_to_copy_psi_from_device_to_host) == 0
  PlotCRP();
  if H3Data.options.calculate_reaction_probabilities == 1
    fprintf(' Save reaction probabilities to %s\n', H3Data.CRP.mat_file);
    CRP = H3Data.CRP;
    save(H3Data.CRP.mat_file, 'CRP');
  end
end

if mod(H3Data.time.steps, H3Data.options.steps_to_copy_psi_from_device_to_host) == 0
  PlotPotWave()
end


return

fprintf(' From C2Matlab\n')

if mod(H3Data.time.steps, 50) == 0
  if H3Data.CRP.calculate_CRP == 1
    CRP = H3Data.CRP;
    save(H3Data.options.CRPMatFile, 'CRP');
  end
end

if mod(H3Data.time.steps, 50) == 0
  PlotPotWave()
end

return

