// (c) 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <iostream>
#include <armadillo>

#include "../include/inputs.h"
#include "../include/report.h"
#include "../include/grid.h"
#include "../include/planets.h"
#include "../include/sizes.h"
#include "../include/fill_grid.h"

using namespace arma;

void Grid::init_geo_grid(Planets planet, Inputs input, Report &report) {

  std::string function="Grid::init_geo_grid";
  static int iFunction = -1;
  report.enter(function, iFunction);  

  // This is just an example:

  Inputs::grid_input_struct grid_input = input.get_grid_inputs();

  long iLon, iLat, iAlt, index;
  float longitude, latitude, altitude;

  // Longitudes:
  // - Make a 1d vector
  // - copy it into the 3d cube
  fvec lon1d(nLons);
  float dlon = (grid_input.lon_max - grid_input.lon_min) / nGeoLons;
  for (iLon=0; iLon < nLons; iLon++)
    lon1d(iLon) = (iLon-nGCs+0.5) * dlon;

  for (iLat=0; iLat < nLats; iLat++) {
    for (iAlt=0; iAlt < nAlts; iAlt++) {
      geoLon_scgc.subcube(0,iLat,iAlt,nLons-1,iLat,iAlt) = lon1d;
    }
  }

  // Latitudes:
  // - Make a 1d vector
  // - copy it into the 3d cube
  fvec lat1d(nLats);
  float dlat = (grid_input.lat_max - grid_input.lat_min) / nGeoLats;
  for (iLat=0; iLat < nLats; iLat++)
    lat1d(iLat) = grid_input.lat_min + (iLat-nGCs+0.5) * dlat;
  
  for (iLon=0; iLon < nLons; iLon++) {
    for (iAlt=0; iAlt < nAlts; iAlt++) {
      geoLat_scgc.subcube(iLon,0,iAlt,iLon,nLats-1,iAlt) = lat1d;
    }
  }
  
  
  fvec alt1d(nAlts);

  if (grid_input.IsUniformAlt) {
    for (iAlt = 0; iAlt < nAlts; iAlt++)
      alt1d(iAlt) =
	grid_input.alt_min +
	(iAlt-nGeoGhosts) * grid_input.dalt;
  }
  for (iLon = 0; iLon < nLons; iLon++) {
    for (iLat = 0; iLat < nLats; iLat++) {
      geoAlt_scgc.tube(iLon,iLat) = alt1d;
    }
  }

  float altitudes[nGeoAltsG];

  IsGeoGrid = 1;
  
  // Make a uniform grid in altitude:
  if (grid_input.IsUniformAlt) {
    for (iAlt=0; iAlt < nGeoAltsG; iAlt++) {
      altitudes[iAlt] =
	grid_input.alt_min + float(iAlt-nGeoGhosts)*grid_input.dalt;
    }
  } else {
    // calc_stretched_altitudes(input, planet, altitudes);
  }

  for (iLon = 0; iLon < nGeoLonsG; iLon++) {
    longitude = grid_input.lon_min + (float(iLon-nGeoGhosts)+0.5) * dlon;
    for (iLat = 0; iLat < nGeoLatsG; iLat++) {
      latitude = grid_input.lat_min + (float(iLat-nGeoGhosts)+0.5) * dlat;
      for (iAlt = 0; iAlt < nGeoAltsG; iAlt++) {
	index = ijk_geo_s3gc(iLon,iLat,iAlt);
	geoLon_s3gc[index] = longitude;
	geoLat_s3gc[index] = latitude;
	geoAlt_s3gc[index] = altitudes[iAlt];
      }
    }
  }

  // Calculate the radius, etc:

  fill_grid_radius(planet, report);

  fill_grid_bfield(planet, input, report);
  
  report.exit(function);  

}
