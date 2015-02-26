#include "pde_common.h"

const double DT = 0.005;
const int nStepsPerPixel = 1000, nPixel = 100;

void uxxStep0(SpatialPoint<1,2>& point)
{
    double dx = point.x - point.nbr(0).x;
    double u = point.inputs(0),
           uL = point.nbr(0).inputs(0),
           uR = point.nbr(1).inputs(0);
    point.outputs(0) = u;
    point.outputs(1) = (uL + uR - 2 * u) / (dx * dx);
}

void updateStep0(SpatialPoint<2,2>& point)
{
    double dx = point.x - point.nbr(0).x;
    double u = point.inputs(0),
           uL = point.nbr(0).inputs(0),
           uR = point.nbr(1).inputs(0);
    double uxx = point.inputs(1),
           uxxL = point.nbr(0).inputs(1),
           uxxR = point.nbr(1).inputs(1);
    double conv = (uR*uR - uL*uL) / (4 * dx);
    double diff = ((uL + uxxL) + (uR + uxxR) - 2 * (u + uxx)) / (dx * dx);
    double dudt = -conv - diff;
    point.outputs(0) = u;
    point.outputs(1) = u + 0.5 * DT * dudt;
}

void uxxStep1(SpatialPoint<2,3>& point)
{
    double dx = point.x - point.nbr(0).x;
    double u0 = point.inputs(0),
           u  = point.inputs(1),
           uL = point.nbr(0).inputs(1),
           uR = point.nbr(1).inputs(1);
    point.outputs(0) = u0;
    point.outputs(1) = u;
    point.outputs(2) = (uL + uR - 2 * u) / (dx * dx);
}

void updateStep1(SpatialPoint<3,1>& point)
{
    double dx = point.x - point.nbr(0).x;
    double u0 = point.inputs(0);
    double u = point.inputs(1),
           uL = point.nbr(0).inputs(1),
           uR = point.nbr(1).inputs(1);
    double uxx = point.inputs(2),
           uxxL = point.nbr(0).inputs(2),
           uxxR = point.nbr(1).inputs(2);
    double conv = (uR*uR - uL*uL) / (4 * dx);
    double diff = ((uL + uxxL) + (uR + uxxR) - 2 * (u + uxx)) / (dx * dx);
    double dudt = -conv - diff;
    point.outputs(0) = u0 + DT * dudt;
}

void init(SpatialPoint<0,1>& point) {
    const double PI = atan(1.0) * 4;
    point.outputs(0) = cos(point.x / 128. * 19 * PI) * 2.;
}

