#include "pde_common.h"

const double DT = 0.005;
const int nStepsPerPixel = 1000, nPixel = 1000;

void uxxStep0(const LocalInputs1D<1>& inputs,
              LocalOutputs1D<2>& outputs,
              const LocalMesh& mesh) {
    double u = inputs[0],
           uL = inputs[0].nbr(0),
           uR = inputs[0].nbr(1);
    outputs[0] = u;
    outputs[1] = (uL + uR - 2 * u) / (mesh.dx * mesh.dx);
}

void updateStep0(const LocalInputs1D<2>& inputs,
                 LocalOutputs1D<2>& outputs,
                 const LocalMesh& mesh) {
    double u = inputs[0],
           uL = inputs[0].nbr(0),
           uR = inputs[0].nbr(1);
    double uxx = inputs[1],
           uxxL = inputs[1].nbr(0),
           uxxR = inputs[1].nbr(1);
    double conv = (uR*uR - uL*uL) / (4 * mesh.dx);
    double diff = ((uL + uxxL) + (uR + uxxR) - 2 * (u + uxx))
                / (mesh.dx * mesh.dx);
    double dudt = -conv - diff;
    outputs[0] = u;
    outputs[1] = u + 0.5 * DT * dudt;
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
