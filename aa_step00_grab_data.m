
close all;clear;clc

%% grab data files from Zenodo repository
% Evergreen Zenodo concept record: the API redirects to the newest published
% version, so this never needs editing when a new dataset version is uploaded
% (concept DOI: 10.5281/zenodo.17388349).

concept_recid = '17388349';

filenames = {'ASIT2019_EPSS_directional_spectra.nc','ASIT2019_supporting_environmental_observations.nc','ASIT2019_supporting_data_compilation.nc','ASIT2019_wave_spectra_stats_timeseries_empirical_gain.nc'};

% Resolve the concept record to the latest version, mapping filename -> URL
record = webread(['https://zenodo.org/api/records/' concept_recid],weboptions('Timeout',30));
record_files = record.files;
file_keys = strings(numel(record_files),1);
file_urls = strings(numel(record_files),1);
for n = 1:numel(record_files)
    if iscell(record_files); f = record_files{n}; else; f = record_files(n); end
    file_keys(n) = string(f.key);
    if isfield(f.links,'content')
        file_urls(n) = string(f.links.content);
    else
        file_urls(n) = string(f.links.self);
    end
end
disp(['Latest Zenodo record: ' num2str(record.id) ' (' num2str(numel(record_files)) ' files available)'])

% Download without timeout (files range up to several hundred MB)
download_opts = weboptions('Timeout',Inf);

for n = 1:length(filenames)

    fn = filenames{n};

    if exist(['data/' fn],'file')

        disp([fn ' already in data/'])

    else

        ind = find(file_keys == fn,1);

        if isempty(ind)

            warning([fn ' not in latest Zenodo record; skipping'])
            continue

        end

        disp(['Downloading ' fn ' from Zenodo...'])
        websave(fn,file_urls(ind),download_opts);
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
    copyfile('data/gshhs_f.b','codes/m_map/data/gshhs_f.b')

end
