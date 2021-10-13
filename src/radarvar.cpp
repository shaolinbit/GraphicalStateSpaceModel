#include "radarvar.h"

namespace gtsam
{
int RadarWindowIndex;
NonlinearFactorGraph RADARWindowgraph;
NonlinearFactorGraph RADARPriorigraph;
 Values RADARvalues;
std::vector<int> RADARKeys;
std::vector<Matrix> RADARmnresults;
int RADARHKey;
int RADARXKey;
int RADARXVKey;
};
