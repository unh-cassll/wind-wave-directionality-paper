
close all;clear;clc

%% grab data files from Zenodo repository

base_url = 'https://zenodo.org/records/18881301/files/';

filenames = {'ASIT2019_EPSS_directional_spectra.nc','ASIT2019_supporting_environmental_observations.nc','ASIT2019_supporting_data_compilation.nc','ASIT2019_wave_spectra_stats_timeseries_empirical_gain.nc'};

for n = 1:length(filenames)

    fn = filenames{n};

    if exist(['data/' fn],'file')

        disp([fn ' already in data/'])

    else

        disp(['Downloading ' fn ' from Zenodo...'])
        websave(fn,[base_url fn]);
        movefile(fn,'data/');

    end

end

%% grab M_MAP toolbox

remote_name = 'https://www.eos.ubc.ca/%7Erich/m_map1.4.zip';

if exist('codes/m_map/','dir')

    disp('M_MAP already retrieved');

else

    websave('m_map.zip',remote_name);
    unzip('m_map.zip','codes/');
    delete('m_map.zip');

end
