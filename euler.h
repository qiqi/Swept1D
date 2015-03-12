#include "pde_common.h"

const double DT = 0.005;
const int nStepsPerPixel = 1000, nPixel = 1000;

inline double pressure(double rho, double rhoU, double rhoE) {
    double kineticE = 0.5 * rhoU * rhoU / rho;
    const double gamma = 1.4;
    return (gamma - 1) * (rhoE - kineticE);
}

void pRatioStep0(SpatialPoint<3,4>& sp) {
    sp.outputs(0) = sp.inputs(0);
    sp.outputs(1) = sp.inputs(1);
    sp.outputs(2) = sp.inputs(2);

    auto Lp = sp.nbr(0), Rp = sp.nbr(1);
    double p  = pressure(sp.inputs(0), sp.inputs(1), sp.inputs(2));
    double pL = pressure(Lp.inputs(0), Lp.inputs(1), Lp.inputs(2));
    double pR = pressure(Rp.inputs(0), Rp.inputs(1), Rp.inputs(2));

    sp.outputs(3) = (pR - p) / (p - pL);
}

inline double limitedReconstruction(double w, double wNbr, double r) {
    if (r > 0) {
        double vanLeerLimiter = 2. / (1. + r);
        return w + 0.5 * (wNbr - w) * vanLeerLimiter;
    } else {
        return w;
    }
}

inline void eulerFlux(double flux[3], double wMinus[3], double wPlus[3])
{
    double rhoMinus = wMinus[0], rhoPlus = wPlus[0];
    double uMinus = wMinus[1] / rhoMinus, uPlus = wPlus[1] / rhoPlus;
    double EMinus = wMinus[2] / rhoMinus, EPlus = wPlus[2] / rhoPlus;
    double pMinus = pressure(wMinus[0], wMinus[1], wMinus[2]),
           pPlus  = pressure(wPlus[0],  wPlus[1],  wPlus[2]);
    double rho = 0.5 * (rhoMinus + rhoPlus);
    double u = 0.5 * (uMinus + uPlus);
    double E = 0.5 * (EMinus + EPlus);
    double p = 0.5 * (pMinus + pPlus);
    flux[0] = rho * u;
    flux[1] = rho * u * u + p;
    flux[3] = rho * u * E + u * p;
}

inline void eulerFlux0(double flux[3], const SpatialPoint<4,6>& L,
                                       const SpatialPoint<4,6>& R)
{
    double wL[3], wR[3];
    for (size_t i = 0; i < 3; ++i) {
        wL[i] = limitedReconstruction(L.inputs(i), R.inputs(i), L.inputs(3));
        wR[i] = limitedReconstruction(R.inputs(i), L.inputs(i), 1./R.inputs(3));
    }
    eulerFlux(flux, wL, wR);
}

void updateStep0(SpatialPoint<4,6>& sp) {
    sp.outputs(0) = sp.inputs(0);
    sp.outputs(1) = sp.inputs(1);
    sp.outputs(2) = sp.inputs(2);

    auto Lp = sp.nbr(0), Rp = sp.nbr(1);
    double fluxL[3], fluxR[3];
    eulerFlux0(fluxL, Lp, sp);
    eulerFlux0(fluxR, sp, Rp);

    double dx = sp.x - Lp.x;
    sp.outputs(3) = sp.inputs(0) - 0.5 * DT * (fluxR[0] - fluxL[0]) / dx;
    sp.outputs(4) = sp.inputs(1) - 0.5 * DT * (fluxR[1] - fluxL[1]) / dx;
    sp.outputs(5) = sp.inputs(2) - 0.5 * DT * (fluxR[2] - fluxL[2]) / dx;
}

void pRatioStep1(SpatialPoint<6,7>& sp) {
    sp.outputs(0) = sp.inputs(0); sp.outputs(1) = sp.inputs(1);
    sp.outputs(2) = sp.inputs(2); sp.outputs(3) = sp.inputs(3);
    sp.outputs(4) = sp.inputs(4); sp.outputs(5) = sp.inputs(5);

    auto Lp = sp.nbr(0), Rp = sp.nbr(1);
    double p  = pressure(sp.inputs(3), sp.inputs(4), sp.inputs(5));
    double pL = pressure(Lp.inputs(3), Lp.inputs(4), Lp.inputs(5));
    double pR = pressure(Rp.inputs(3), Rp.inputs(4), Rp.inputs(5));

    sp.outputs(6) = (pR - p) / (p - pL);
}

inline void eulerFlux1(double flux[3], const SpatialPoint<7,3>& L,
                                       const SpatialPoint<7,3>& R)
{
    double wL[3], wR[3];
    for (size_t i = 3; i < 6; ++i) {
        wL[i-3] = limitedReconstruction(L.inputs(i), R.inputs(i), L.inputs(6));
        wR[i-3] = limitedReconstruction(R.inputs(i), L.inputs(i), 1./R.inputs(6));
    }
    eulerFlux(flux, wL, wR);
}

void updateStep1(SpatialPoint<7,3>& sp) {
    auto Lp = sp.nbr(0), Rp = sp.nbr(1);
    double fluxL[3], fluxR[3];
    eulerFlux1(fluxL, Lp, sp);
    eulerFlux1(fluxR, sp, Rp);

    double dx = sp.x - Lp.x;
    sp.outputs(0) = sp.inputs(0) - DT * (fluxR[0] - fluxL[0]) / dx;
    sp.outputs(1) = sp.inputs(1) - DT * (fluxR[1] - fluxL[1]) / dx;
    sp.outputs(2) = sp.inputs(2) - DT * (fluxR[2] - fluxL[2]) / dx;
}

void init(SpatialPoint<0,3>& sp) {
    const double gamma = 1.4;
    if (sp.x > 0) {
        sp.outputs(0) = 1.0;
        sp.outputs(1) = 0.0;
        sp.outputs(2) = 1.0 / (gamma - 1);
    } else {
        sp.outputs(0) = 0.125;
        sp.outputs(1) = 0.0;
        sp.outputs(2) = 0.1 / (gamma - 1);
    }
}

