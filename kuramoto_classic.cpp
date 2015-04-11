#include <cmath>
#include <cassert>
#include <ctime>
#include "kuramoto.h"
#include "pde_classic.h"

int main(int argc, char*argv[])
{
    if (argc != 2) {
        exit(-1);
    }
    size_t numGridPerProc = atoi(argv[1]);
    ClassicDiscretization1D disc(numGridPerProc, 0.5, init);
    disc.colorMap.red.set(0, -2., 2.);
    disc.colorMap.green.set(0, -2., 2.);
    disc.colorMap.blue.set(0, -2., 2.);

    std::clock_t startTime = std::clock();

    for (int iPixel = 0; iPixel < nPixel; ++iPixel) {
        for (int iStep = 0; iStep < nStepsPerPixel; ++iStep) {
            disc.applyOp(uxxStep0);
            disc.applyOp(updateStep0);
            disc.applyOp(uxxStep1);
            disc.applyOp(updateStep1);
        }
        disc.variablesToColor(iPixel);
        disc.writePng("kura");
    }

    if (iProc() == 0) {
        std::clock_t endTime = std::clock();
        double totalTime = (endTime - startTime) / (double)CLOCKS_PER_SEC;
        std::cout << totalTime * 1000000 / nStepsPerPixel / nPixel / 4
                  << " microseconds per SubStep" << std::endl;
    }
}
