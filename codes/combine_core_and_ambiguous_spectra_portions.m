
function [core_spect,ambiguous_spect] = combine_core_and_ambiguous_spectra_portions(in_spect,in_direction,in_scale,roi_directional_extrema_core,roi_directional_extrema_ambiguous)

directional_res = median(diff(in_direction));
inds_per_180 = floor(pi/directional_res);

copy_spect_core = in_spect;
copy_spect_ambiguous = 0*in_spect;

inds_nan = isnan(roi_directional_extrema_core(:,1)) | isnan(roi_directional_extrema_core(:,2));
roi_directional_extrema_core(inds_nan,:) = NaN;

ind_max_core = find(~isnan(roi_directional_extrema_core(:,1)),1,'last');
ind_max_ambiguous = find(~isnan(roi_directional_extrema_ambiguous(:,1)),1,'last');

for i = 1:ind_max_core

    ind_left_core = find(in_direction<roi_directional_extrema_core(i,1),1,'last');
    ind_right_core = find(in_direction>roi_directional_extrema_core(i,2),1,'first');

    inds_core = ind_left_core:ind_right_core;

    ind_left_ambiguous = find(in_direction<roi_directional_extrema_ambiguous(i,1),1,'last');
    ind_right_ambiguous = find(in_direction>roi_directional_extrema_ambiguous(i,2),1,'first');

    inds_ambiguous = ind_left_ambiguous:ind_right_ambiguous;

    core_ambig_diff = median(inds_core) - median(inds_ambiguous);

    if core_ambig_diff == 0

        offset_sign = 1;

    else

        offset_sign = core_ambig_diff/abs(core_ambig_diff);

    end

    copy_spect_ambiguous(i,inds_ambiguous + offset_sign*inds_per_180) = in_spect(i,inds_ambiguous);

    copy_spect_core(i,1:ind_left_core-1) = 0;
    copy_spect_core(i,ind_right_core+1:end) = 0;

end

copy_spect_core(ind_max_core:end,:) = 0;
copy_spect_ambiguous(ind_max_ambiguous:end,:) = 0;

core_spect = copy_spect_core;
ambiguous_spect = copy_spect_ambiguous;