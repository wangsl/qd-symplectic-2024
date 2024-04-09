
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

