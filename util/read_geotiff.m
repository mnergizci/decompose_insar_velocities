function [lon,lat,data,dx,dy] = read_geotiff(tif_file,out_type)
%% read_geotiff.m
% Read geotif file and output data, coord vectors, and grid spacing.
%
% INPUT:                                                           
%   tif_file: input string of tif file
%   out_type: force data type of output (currently just single or double)
% OUTPUT:    
%   lon, lat: vectors of coordinates, lat counts downwards
%   data: loaded tif data
%   dx, dy: posting spacing of lon and lat
% Andrew Watson     24-08-2021

% open geotiff
try
    [data,georef] = readgeoraster(tif_file);
catch
    error(['Cannot open ', tif_file])
    return
end
    
% convert to double
data = double(data);

try
if strcmp(georef.RasterInterpretation,'cells')
    disp('Tif is in cell format, converting to postings.')
    if strcmpi(georef.CoordinateSystemType,'geographic')
        dx = georef.CellExtentInLongitude;
        dy = georef.CellExtentInLatitude;
        latlim = [georef.LatitudeLimits(1)+dy/2 georef.LatitudeLimits(2)-dy/2];
        lonlim = [georef.LongitudeLimits(1)+dx/2 georef.LongitudeLimits(2)-dx/2];
    else  % Sometimes GDAL writes the tifs slightly differently.....
        dx = georef.CellExtentInWorldX;
        dy = georef.CellExtentInWorldY;
        latlim = [georef.YWorldLimits(1)+dy/2 georef.YWorldLimits(2)-dy/2];
        lonlim = [georef.XWorldLimits(1)+dx/2 georef.XWorldLimits(2)-dx/2];
    end
    rasterSize = georef.RasterSize;
    georef = georefpostings(latlim,lonlim,rasterSize);
end
catch err
    throw(err)
end

% grid spacing
dx = georef.SampleSpacingInLongitude;
dy = georef.SampleSpacingInLatitude;

% get coords and grid spacing
lon = georef.LongitudeLimits(1) : dx : georef.LongitudeLimits(2);
lat = georef.LatitudeLimits(2) : -dy : georef.LatitudeLimits(1);

if nargin == 2 & out_type == 'single'
    data = single(data);
end

end

