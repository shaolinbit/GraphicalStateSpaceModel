#ifndef RADARVAR_H_INCLUDED
#define RADARVAR_H_INCLUDED

#include <gtsam/nonlinear/NonlinearFactor.h>
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <gtsam/nonlinear/Values.h>

namespace gtsam
{

extern int RadarWindowIndex;



#define RADAR_CONSTANTV_INTEGER 500000
#define RADAR_CONSTANTH_INTEGER 1000000
#define CONSTANT_RADAR_TS  0.05
#define SIN_WINDOW_LENGTH 10
extern NonlinearFactorGraph RADARWindowgraph;
extern NonlinearFactorGraph RADARPriorigraph;
extern Values RADARvalues;

extern std::vector<int> RADARKeys;
extern std::vector<Matrix> RADARmnresults;

extern int RADARHKey;
extern int RADARXKey;
extern int RADARXVKey;
};


#endif // GSSMVAR_H_INCLUDED
