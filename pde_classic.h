#ifndef PDE_CLASSIC_H
#define PDE_CLASSIC_H

#include<vector>
#include<cstring>
#include<cassert>
#include<cstdlib>
#include<cmath>
#include<mpi.h>
#include"pde_common.h"
#include"PngWriter.hpp"

// ====== implementation of ClassicDiscretization1D ======= //

class ClassicDiscretization1D {
    private:
    size_t numGrids_;
    double dx_;
    double x0_;

    int numVariables_, numLastOutput_;
    double * variablesData_, * inputs_, * outputs_;

    SpatialPoint<0,0> * spatialPoints_;

    char* pngFilename_;
    PngWriter png_;

    void commonInit_() {
        MPI_Init(0, 0);

        int iProc;
        MPI_Comm_rank(MPI_COMM_WORLD, &iProc);
        x0_ = numGrids_ * dx_ * iProc;
    }

    public:
    virtual ~ClassicDiscretization1D() {
        if (variablesData_) {
            delete variablesData_;
        }
        for (size_t i = 0; i < numGrids_ + 2; ++i) {
            spatialPoints_[i].~SpatialPoint();
        }
        free(spatialPoints_);
        MPI_Finalize();
    }

    public:
    template<size_t numVar>
    ClassicDiscretization1D(int numGrids, double dx,
         void (*localOperator)(SpatialPoint<0, numVar>&))
    :
        numGrids_(numGrids), dx_(dx),
        numVariables_(numVar), numLastOutput_(numVar),
        pngFilename_(0), png_(0,0)
    {
        commonInit_();

        variablesData_ = new double[2 * numVar * (numGrids_ + 2)];
        inputs_ = variablesData_;
        outputs_ = variablesData_ + numVar;

        spatialPoints_ = (SpatialPoint<0,0>*)
            malloc((numGrids_ + 2) * sizeof(SpatialPoint<0,0>));
        for (size_t i = 0; i < numGrids_ + 2; ++i) {
            auto p = spatialPoints_ + i;
            auto pL = (i == 0) ? 0 : spatialPoints_ + i - 1;
            auto pR = (i == numGrids_ + 1) ? 0 : spatialPoints_ + i + 1;
            size_t iShift = i * numVariables_ * 2;
            new(p) SpatialPoint<0,0>(
                    x0_ + (i - 1) * dx_, iShift, inputs_, outputs_, pL, pR);
            localOperator(*(SpatialPoint<0, numVar>*)p);
        }

        ClassicSyncer1D sync(
                pOutput_(1), pOutput_(numGrids_),
                pOutput_(0), pOutput_(numGrids_+1), numVar);
        sync.waitTillDone();
    }

    private:

    void resizeWorkspace_(size_t newNumVariables) {
        double * newVariablesData
            = new double[2 * newNumVariables * (numGrids_ + 2)];
        double * newInputs = newVariablesData;
        double * newOutputs = newVariablesData + newNumVariables;
        for (size_t i = 0; i < numGrids_ + 2; ++ i) {
            void * dest = newOutputs + i * 2 * newNumVariables;
            void * src = outputs_ + i * 2 * numVariables_;
            memcpy(dest, src, sizeof(double) * numVariables_);
        }
        delete[] variablesData_;
        numVariables_ = newNumVariables;
        variablesData_ = newVariablesData;
        inputs_ = newInputs;
        outputs_ = newOutputs;

        for (size_t i = 0; i < numGrids_ + 2; ++i) {
            auto p = spatialPoints_ + i;
            auto pL = (i == 0) ? 0 : p - 1;
            auto pR = (i == numGrids_ + 1) ? 0 : p + 1;
            size_t iShift = i * numVariables_ * 2;
            p->~SpatialPoint();
            new(p) SpatialPoint<0,0>(
                    x0_ + (i - 1) * dx_, iShift, inputs_, outputs_, pL, pR);
        }
    }

    double * pInput_(size_t iGrid) {
        assert(iGrid <= numGrids_ + 1);
        return inputs_ + iGrid * 2 * numVariables_;
    }
    double * pOutput_(size_t iGrid) {
        assert(iGrid <= numGrids_ + 1);
        return outputs_ + iGrid * 2 * numVariables_;
    }

    public:
    template<size_t numInput, size_t numOutput>
    void applyOp(void (*localOperator)(SpatialPoint<numInput,numOutput>& point))
    {
        assert(numInput == numLastOutput_);
        numLastOutput_ = numOutput;
        if (numOutput > numVariables_) {
            resizeWorkspace_(numOutput);
        }

        double * tmp = inputs_;
        inputs_ = outputs_;
        outputs_ = tmp;

        localOperator((SpatialPoint<numInput, numOutput>&)spatialPoints_[1]);
        localOperator((SpatialPoint<numInput, numOutput>&)spatialPoints_[numGrids_]);

        ClassicSyncer1D sync(
                pOutput_(1), pOutput_(numGrids_),
                pOutput_(0), pOutput_(numGrids_+1), numOutput);

        for (size_t iGrid = 2; iGrid < numGrids_; ++iGrid) {
            localOperator((SpatialPoint<numInput, numOutput>&)spatialPoints_[iGrid]);
        }
        sync.waitTillDone();
    }

    // ------------ write to png file ------------- //

    public:

    struct {
        class VarColor_ {
            private:
            size_t iVar_;
            double low_, high_;

            public:
            VarColor_() : iVar_(0), low_(0.0), high_(INFINITY) {}
            void set(int iVar, double low, double high) {
                iVar_ = iVar;
                low_ = low;
                high_ = high;
            }

            void assertIVarLessThan(size_t numVars) {
                assert(iVar_ < numVars);
            }

            double map(double* pVal) {
                return (pVal[iVar_] - low_) / (high_ - low_);
            }
        } red, green, blue;

        void assertIVarLessThan(int numVars) {
            red.assertIVarLessThan(numVars);
            green.assertIVarLessThan(numVars);
            blue.assertIVarLessThan(numVars);
        }
    } colorMap;

    void variablesToColor(int iStep)
    {
        for (size_t iGrid = 0; iGrid < numGrids_; ++iGrid) {
            colorMap.assertIVarLessThan(numVariables_);
            double* pGrid = outputs_ + 2 * iGrid * numVariables_;
            double r = colorMap.red.map(pGrid),
                   g = colorMap.green.map(pGrid),
                   b = colorMap.blue.map(pGrid);
            png_.set(iGrid, iStep, r, g, b);
        }
    }

    void writePng() {
        if (pngFilename_) {
            png_.write(pngFilename_);
        }
    }

    void writePng(const char* filename) {
        if (pngFilename_) {
            delete[] pngFilename_;
        }

        pngFilename_ = new char[strlen(filename) + strlen(mpiRankString()) + 5];
        strcpy(pngFilename_, filename);
        strcat(pngFilename_, mpiRankString());
        strcat(pngFilename_, ".png");
        writePng();
    }
};

#endif
