#include "pde_common.h"

const double DT = 0.005;
const int nStepsPerPixel = 1000, nPixel = 10000;

inline double pressure(double rho, double rhoU, double rhoE) {
    double kineticE = 0.5 * rhoU * rhoU / rho;
    const double gamma = 1.4;
    return (gamma - 1) * (rhoE - kineticE);
}

void pRatioStep0(const LocalInputs1D<3>& inputs,
                         LocalOutputs1D<5>& outputs,
                         const LocalMesh&) {
    outputs[0] = inputs[0];
    outputs[1] = inputs[1];
    outputs[2] = inputs[2];

    double p  = pressure(inputs[0],        inputs[1],        inputs[2]);
    double pL = pressure(inputs[0].nbr(0), inputs[1].nbr(0), inputs[2].nbr(0));
    double pR = pressure(inputs[0].nbr(1), inputs[1].nbr(1), inputs[2].nbr(1));

    outputs[3] = (pR - p) / (p - pL);
}

inline double limitedReconstruction(double w, double wNbr, double r) {
    if (r > 0) {
        double vanLeerLimiter = 2. / (1. + r);
        return w + 0.5 * (wNbr - w) * vanLeerLimiter;
    } else {
        return w;
    }
}

inline void eulerFlux(double flux[3], double wMinus[3], double wPlus[3]) {
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

void updateStep0(const LocalInputs1D<2>& inputs,
                 LocalOutputs1D<2>& outputs,
                 const LocalMesh& mesh)
{
    outputs[0] = inputs[0];
    outputs[1] = inputs[1];
    outputs[2] = inputs[2];

    double wLminus[3], wLplus[3];
    for (size_t i = 0; i < 3; ++i) {
        wLminus[i] = limitedReconstruction(
                     inputs[i].nbr(0), inputs[i], inputs[3].nbr(0));
        wLplus[i]  = limitedReconstruction(
                     inputs[i], inputs[i].nbr(0), 1./ inputs[3]);
    }
    double fluxL[3];
    eulerFlux(fluxL, wLminus, wLplus);

    double wRminus[3], wRplus[3];
    for (size_t i = 0; i < 3; ++i) {
        wRminus[i] = limitedReconstruction(
         inputs[i], inputs[i].nbr(1), inputs[3]);
        wRplus[i]  = limitedReconstruction(
         inputs[i].nbr(1), inputs[i], 1./ inputs[3].nbr(1));
    }
    double fluxR[3];
    eulerFlux(fluxR, wRminus, wRplus);

    outputs[3] = inputs[0] - 0.5 * DT * (fluxR[0] - fluxL[0]) / mesh.dx;
    outputs[4] = inputs[1] - 0.5 * DT * (fluxR[1] - fluxL[1]) / mesh.dx;
    outputs[5] = inputs[2] - 0.5 * DT * (fluxR[2] - fluxL[2]) / mesh.dx;
}

void uxxStep1(const LocalInputs1D<2>& inputs,
              LocalOutputs1D<3>& outputs,
              const LocalMesh& mesh) {
    double u0 = inputs[0],
           u  = inputs[1],
           uL = inputs[1].nbr(0),
           uR = inputs[1].nbr(1);
    outputs[0] = u0;
    outputs[1] = u;
    outputs[2] = (uL + uR - 2 * u) / (mesh.dx * mesh.dx);
}

void updateStep1(const LocalInputs1D<3>& inputs,
                 LocalOutputs1D<1>& outputs,
                 const LocalMesh& mesh) {
    double u0 = inputs[0];
    double u  = inputs[1],
           uL = inputs[1].nbr(0),
           uR = inputs[1].nbr(1);
    double uxx  = inputs[2],
           uxxL = inputs[2].nbr(0),
           uxxR = inputs[2].nbr(1);
    double conv = (uR*uR - uL*uL) / (4 * mesh.dx);
    double diff = ((uL + uxxL) + (uR + uxxR) - 2 * (u + uxx))
                / (mesh.dx * mesh.dx);
    double dudt = -conv - diff;
    outputs[0] = u0 + DT * dudt;
}

void init(LocalOutputs1D<1>& u, const LocalMesh& mesh) {
    const double PI = atan(1.0) * 4;
    u[0] = cos(mesh.x / 128. * 19 * PI) * 2.;
}
