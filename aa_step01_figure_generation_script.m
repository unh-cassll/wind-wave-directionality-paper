% Figure generation command center
%
% "The directionality of wind waves
% from equilibrium to gravity-capillary scales"
%
% N. Laxague and co-authors, 2026
%
% First argument: figure number(s) to generate [default: 1:num_figs]
% Second argument: print flag [true or false]
% Third argument: wind source override, for the four-way source comparison in
%                 aa_step02_wind_source_comparison [default: the preamble value]
% Fourth argument: wave-age bin spacing override ('quantile', 'log' or
%                  'linear'), for sensitivity reruns [default: the preamble value]
%
function aa_step01_figure_generation_script(varargin)

addpath codes/
addpath codes/m_map/
addpath data/

fsize = 24;
figure_style(fsize)

full_pos = [50 100 1100 1000];

raster_dpi = 300;
raster_dpi_string = ['-r' num2str(raster_dpi)];

num_figs = 19;

if length(varargin) < 1
    fig_list = 1:num_figs;
else
    fig_list = varargin{1};
end
if length(varargin) < 2
    print_flag = false;
else
    print_flag = varargin{2};
end
if length(varargin) < 3
    wind_source_override = '';
else
    wind_source_override = varargin{3};
end
if length(varargin) < 4
    wave_age_spacing_override = '';
else
    wave_age_spacing_override = varargin{4};
end
fig_folder = 'figs/';

% Global figure settings
example_run_ind = 109;            % spectra/spreading example (figs 5, 8, B1; initial-submission case)
breaking_example_run_ind = 173;   % breaking example (figs 9, 10): run 109 has no Lambda data; 114 = v2 breaking case (theta_br +11.71)
dk = 2.16;
k_low = 2*dk;
k_high = 100;
f_high = 7.5;
nu_high = 4;
U_low = 1;
U_high = 13;
dU = 2;
wave_age_lims = [10 90];
wave_age_bin_spacing = 'log';  % 'quantile' (equal count), 'log' (wave age ~ 1/U10, skewed), or 'linear'
n_wave_age_bins = 5;                % quantile and log spacing
d_wave_age = 15;                    % linear spacing only
wind_source_choice = 'EC_COARE_ASIT';  % 'EC_ASIT', 'COARE_ASIT', 'EC_COARE_ASIT' (EC with ASIT
                                       % COARE backfill)

diving_board_heading_deg = 255;   % instrument boom heading [deg true]

ASIT_EC_wind_gate_deg = 90;

ASIT_COARE_wind_gate_deg = 135;

equilibrium_slope_cut = -4.3;     % reject runs whose fitted slope below the break is already at the
                                  % saturation value: no equilibrium range, so f_eq_end is the swell flank
ustar_cut = 0;                    

recompute_derived_products = false;
if ~isempty(wind_source_override)
    wind_source_choice = wind_source_override;
end
if ~isempty(wave_age_spacing_override)
    wave_age_bin_spacing = wave_age_spacing_override;
end
save('data/global_figure_settings.mat','example_run_ind','breaking_example_run_ind','k_high','k_low','f_high','nu_high','U_low','U_high','dU','wave_age_lims','d_wave_age','wave_age_bin_spacing','n_wave_age_bins','wind_source_choice','equilibrium_slope_cut','diving_board_heading_deg','ASIT_EC_wind_gate_deg','ASIT_COARE_wind_gate_deg','ustar_cut','recompute_derived_products')

% Example-case values consumed by the manuscript via \input fragments
pk = load('data/ASIT2019_peak_wave_phase_speed.mat');
wind = wind_forcing();
U10_coare = wind.U10;
us_ec = wind.ustar;
write_tex_macros('figs/tex/values_example_cases.tex', ...
    {'ExRunID','ExUten','ExWaveAge','BrExRunID','BrExUten','BrExWaveAge'}, ...
    {sprintf('%d',example_run_ind), ...
     sprintf('%.2f',U10_coare(example_run_ind)), ...
     sprintf('%.1f',pk.C_p(example_run_ind)/us_ec(example_run_ind)), ...
     sprintf('%d',breaking_example_run_ind), ...
     sprintf('%.2f',U10_coare(breaking_example_run_ind)), ...
     sprintf('%.1f',pk.C_p(breaking_example_run_ind)/us_ec(breaking_example_run_ind))})

%% Handle input figure list

for fignum = fig_list

    switch fignum

        %% Figure 01 - map of Cape Cod, zoom into ASIT location

        case 1

            figure(fignum);clf
            set(fignum,'Position',full_pos.*[1 1 1 0.55])

            generate_MV_ASIT_map(fignum,fsize)

            if print_flag

                figure(fignum)
                pause(0.5)
                print([fig_folder 'CapeCod_MV_ASIT_bathy.svg'],'-dsvg')
                pause(0.5)
                print([fig_folder 'CapeCod_MV_ASIT_bathy.png'],'-dpng',raster_dpi_string)

            end


            %% Figure 02 - windrose (combine with ASIT schematic)

        case 2

            figure(fignum);clf
            set(fignum,'Position',full_pos.*[1 1 0.5 0.5])
            wind_speed_direction_rose(fignum,fsize);

            if print_flag

                figure(fignum)
                pause(0.5)
                print([fig_folder 'ASIT_Winter2019-2020_windrose.svg'],'-dsvg')
                pause(0.5)
                print([fig_folder 'ASIT_Winter2019-2020_windrose.png'],'-dpng',raster_dpi_string)

            end

            %% Figure 03 - stacked timeseries/histograms of key variables

        case 3

            figure(fignum);clf
            set(fignum,'Position',full_pos.*[1 1 1 1])

            stacked_timeseries(fignum,fsize)

            if print_flag

                figure(fignum)
                pause(0.5)
                print([fig_folder 'stacked_timeseries_histograms.svg'],'-dsvg')
                pause(0.5)
                print([fig_folder 'stacked_timeseries_histograms.png'],'-dpng',raster_dpi_string)

            end

            %% Figure 04 - directional spectrum and direction-integrated outputs

        case 4

            figure(fignum);clf
            set(fignum,'Position',full_pos.*[1 1 0.5 1])

            wavenumber_directional_and_omni_spect(fignum,fsize)

            if print_flag

                figure(fignum)
                pause(0.5)
                print([fig_folder 'wavenumber_directional_and_omni_spect.svg'],'-dsvg')
                pause(0.5)
                print([fig_folder 'wavenumber_directional_and_omni_spect.png'],'-dpng',raster_dpi_string)

            end

            %% Figure 05 - frequency and wavenumber omnispect

        case 5

            figure(fignum);clf
            set(fignum,'Position',full_pos)

            binned_omnispect(fignum,fsize)

            if print_flag

                figure(fignum)
                pause(0.5)
                print([fig_folder 'frequency_wavenumber_omnispect.svg'],'-dsvg')
                pause(0.5)
                print([fig_folder 'frequency_wavenumber_omnispect.png'],'-dpng',raster_dpi_string)

            end

            %% Figure 06 - frequency and wavenumber spectral subranges

        case 6

            figure(fignum);clf
            set(fignum,'Position',full_pos.*[1 1 0.5 1.33])

            frequency_wavenumber_spectral_subranges(fignum,fsize)

            if print_flag

                figure(fignum)
                pause(0.5)
                print([fig_folder 'frequency_wavenumber_spectral_subranges.svg'],'-dsvg')
                pause(0.5)
                print([fig_folder 'frequency_wavenumber_spectral_subranges.png'],'-dpng',raster_dpi_string)

            end

            %% Figure 07 - directional spectra comparison: f, k

        case 7

            figure(fignum);clf
            set(fignum,'Position',full_pos.*[1 1 1 0.75])

            directional_spectra_spreading_delta(fignum,fsize)

            if print_flag

                figure(fignum)
                pause(0.5)
                print([fig_folder 'directional_spectra_spreading_delta.svg'],'-dsvg')
                pause(0.5)
                print([fig_folder 'directional_spectra_spreading_delta.png'],'-dpng',raster_dpi_string)

            end

            %% Figure 08 - breaking crest length distribution and S_ds

        case 8

            figure(fignum);clf
            set(fignum,'Position',full_pos.*[1 1 0.5 1])

            lambda_S_ds_example(fignum,fsize)

            if print_flag

                figure(fignum)
                pause(0.5)
                print([fig_folder 'lambda_S_ds_example.svg'],'-dsvg')
                pause(0.5)
                print([fig_folder 'lambda_S_ds_example.png'],'-dpng',raster_dpi_string)

            end

            %% Figure 09 - short wave spectrum with S_ds contours and wavenumber-integrated quantities (e.g., saturation-range spreading function and S_ds)

        case 9

            figure(fignum);clf
            set(fignum,'Position',full_pos.*[1 1 0.5 1.33])

            spectra_S_ds_contours_example(fignum,fsize)

            if print_flag

                figure(fignum)
                pause(0.5)
                print([fig_folder 'spectra_S_ds_contours_example.svg'],'-dsvg')
                pause(0.5)
                print([fig_folder 'spectra_S_ds_contours_example.png'],'-dpng',raster_dpi_string)

            end

            %% Figure 10 - S_ds(theta), binned by wave age

        case 10

            figure(fignum);clf
            set(fignum,'Position',full_pos.*[1 1 0.5 1])

            S_ds_theta_all(fignum,fsize)

            if print_flag

                figure(fignum)
                pause(0.5)
                print([fig_folder 'S_ds_theta_all.svg'],'-dsvg')
                pause(0.5)
                print([fig_folder 'S_ds_theta_all.png'],'-dpng',raster_dpi_string)

            end

            %% Figure 11 - wave breaking momentum flux and breaking direction

        case 11

            figure(fignum);clf
            set(fignum,'Position',full_pos.*[1 1 0.5 1])

            breaking_wave_momentum_and_energy_flux(fignum,fsize)

            if print_flag

                figure(fignum)
                pause(0.5)
                print([fig_folder 'breaking_direction_rel_wind_by_wave_age.svg'],'-dsvg')
                pause(0.5)
                print([fig_folder 'breaking_direction_rel_wind_by_wave_age.png'],'-dpng',raster_dpi_string)

            end

            %% Figure 12 - wave-relative wind direction with wind speed

        case 12

            figure(fignum);clf
            set(fignum,'Position',full_pos.*[1 1 1 1.1])

            wind_wave_subrange_directions(fignum,fsize)

            if print_flag

                figure(fignum)
                pause(0.5)
                print([fig_folder 'wind_wave_subrange_directions.svg'],'-dsvg')
                pause(0.5)
                print([fig_folder 'wind_wave_subrange_directions.png'],'-dpng',raster_dpi_string)

            end

            %% Figure 13 - scaled wavenumber spectra and transition wavenumbers

        case 13

            figure(fignum);clf
            set(fignum,'Position',full_pos.*[1 1 0.55 1])

            scaled_spectra_and_transition_wavenumbers(fignum,fsize)

            if print_flag

                figure(fignum)
                pause(0.5)
                print([fig_folder 'scaled_spectra_and_transition_wavenumbers.svg'],'-dsvg')
                pause(0.5)
                print([fig_folder 'scaled_spectra_and_transition_wavenumbers.png'],'-dpng',raster_dpi_string)

            end

            %% Figure 14 - normalized transition wavenumber

        case 14

            figure(fignum);clf
            set(fignum,'Position',full_pos.*[1 1 1 0.55])

            normalized_transition_wavenumber(fignum,fsize)

            if print_flag

                figure(fignum)
                pause(0.5)
                print([fig_folder 'normalized_transition_wavenumber.svg'],'-dsvg')
                pause(0.5)
                print([fig_folder 'normalized_transition_wavenumber.png'],'-dpng',raster_dpi_string)

            end

            %% Figure 15-16 - directional spreading, binned by wave age and plotted against normalized wavenumber

        case 15

            figure(fignum);clf
            set(fignum,'Position',full_pos.*[1 1 1 1])

            wave_age_binned_directional_spreading(fignum,fsize)

            if print_flag

                % explicit handle: the generator also creates figure fignum+1,
                % and pause() lets the newer figure steal gcf in batch mode
                h_fig = figure(fignum);
                pause(0.5)
                print(h_fig,[fig_folder 'frequency_wavenumber_directional_spreading_binned.svg'],'-dsvg')
                pause(0.5)
                print(h_fig,[fig_folder 'frequency_wavenumber_directional_spreading_binned.png'],'-dpng',raster_dpi_string)

            end

        case 16

            figure(fignum);clf
            set(fignum,'Position',full_pos.*[1 1 1 1])

            wave_age_binned_directional_spreading(fignum-1,fsize)

            if print_flag

                h_fig = figure(fignum);
                pause(0.5)
                print(h_fig,[fig_folder 'normalized_directional_spreading_binned.svg'],'-dsvg')
                pause(0.5)
                print(h_fig,[fig_folder 'normalized_directional_spreading_binned.png'],'-dpng',raster_dpi_string)

            end

            %% Figure A1

        case 17

            figure(fignum);clf
            set(fignum,'Position',full_pos.*[1 1 1 0.5])

            all_wind_speed_stress_profiles(fignum,fsize)

            if print_flag

                figure(fignum)
                pause(0.5)
                print([fig_folder 'all_wind_speed_stress_profiles.svg'],'-dsvg')
                pause(0.5)
                print([fig_folder 'all_wind_speed_stress_profiles.png'],'-dpng',raster_dpi_string)

            end

            %% Figure B1

        case 18

            figure(fignum);clf
            set(fignum,'Position',full_pos.*[1 0.5 0.5 1.25])

            directional_spectra_spreading_revisit(fignum,fsize)

            if print_flag

                figure(fignum)
                pause(0.5)
                print([fig_folder 'directional_spectra_spreading_revisit.svg'],'-dsvg')
                pause(0.5)
                print([fig_folder 'directional_spectra_spreading_revisit.png'],'-dpng',raster_dpi_string)

            end

            %% Figure C1

        case 19

            figure(fignum);clf
            set(fignum,'Position',full_pos.*[1 0.5 0.5 1])

            spectral_subrange_slope_extraction_example(fignum,fsize)

            if print_flag

                figure(fignum)
                pause(0.5)
                print([fig_folder 'spectral_subrange_slope_extraction_example.svg'],'-dsvg')
                pause(0.5)
                print([fig_folder 'spectral_subrange_slope_extraction_example.png'],'-dpng',raster_dpi_string)

            end

    end

end

