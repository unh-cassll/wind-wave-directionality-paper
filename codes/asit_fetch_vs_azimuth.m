% Fetch from the ASIT position along a given wind coming-from azimuth [deg
% true], by ray-casting against the GSHHS full-resolution coastline of the
% surrounding region (Cape Cod, the Islands, Buzzards Bay, the RI coast).
% Open-ocean azimuths with no land inside the region are capped at 200 km.
%
% The per-degree fetch table is computed once and cached in data/; queries
% interpolate on the cached table. Requires codes/m_map on the path
%
% Nathan Laxague 2026
%
function fetch_km = asit_fetch_vs_azimuth(az_query_deg)

fetch_cap_km = 200;
asit_position = [41.333150 -70.573135];   % [lat lon], as in generate_MV_ASIT_map

cache_name = repo_data_path('ASIT_fetch_vs_azimuth.mat');
coast_name = repo_data_path('ASIT_region_coast.mat');

if exist(cache_name,'file')

    c = load(cache_name);

else

    % Extract the regional coastline once via m_map
    if ~exist(coast_name,'file')
        m_proj('miller','long',[-72.5 -69.2],'lat',[40.6 42.6]);
        m_gshhs_f('save',coast_name);
    end

    coast = load(coast_name);
    ncst = coast.ncst;

    % Local flat-earth coordinates about ASIT [km]
    lat0 = asit_position(1);
    lon0 = asit_position(2);
    xc = (ncst(:,1)-lon0)*111.320*cosd(lat0);
    yc = (ncst(:,2)-lat0)*110.574;

    % Coastline segments, dropping the NaN-separated polygon breaks
    p1 = [xc(1:end-1) yc(1:end-1)];
    p2 = [xc(2:end) yc(2:end)];

    inds_ok = all(isfinite([p1 p2]),2);
    p1 = p1(inds_ok,:);
    p2 = p2(inds_ok,:);

    e = p2 - p1;

    az = 0:359;
    fetch_table = fetch_cap_km*ones(size(az));

    for n = 1:length(az)

        dx = sind(az(n));
        dy = cosd(az(n));

        % Intersect the ray t*[dx dy] with each segment p1 + s*e
        Det = dy*e(:,1) - dx*e(:,2);
        t = (p1(:,2).*e(:,1) - p1(:,1).*e(:,2))./Det;   % distance along ray [km]
        s = (dx*p1(:,2) - dy*p1(:,1))./Det;             % position along segment

        inds_hit = t > 0.01 & s >= 0 & s <= 1;

        if any(inds_hit)
            fetch_table(n) = min([t(inds_hit); fetch_cap_km]);
        end

    end

    save(cache_name,'az','fetch_table','fetch_cap_km','asit_position')
    c = load(cache_name);

end

% Circular interpolation onto the query azimuths
az_wrap = [c.az c.az(1)+360];
fetch_wrap = [c.fetch_table c.fetch_table(1)];

fetch_km = interp1(az_wrap,fetch_wrap,mod(az_query_deg,360));
