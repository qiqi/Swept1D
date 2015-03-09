#include "pde_common.h"

const double DT = 0.00001, DX = 0.1;
const int nStepsPerPixel = 1000, nPixel = 100;

inline void updateStep(SpatialPoint<1,1>& point)
{
    const double CFL = DT / DX / DX;
    double u  = point.inputs(0),
           uL = point.nbr(0).inputs(0),
           uR = point.nbr(1).inputs(0);
    point.outputs(0) = u + CFL * (uL + uR - 2 * u);
}

void init(SpatialPoint<0,1>& point) {
    const double PI = atan(1.0) * 4;
    point.outputs(0) = cos(point.x / 128. * 19 * PI) * 2.;
}
