// This file is part of Diamond1D
// Copyright (C) 2015 Qiqi Wang, qiqi@mit.edu
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as
// published by the Free Software Foundation, either version 3 of the
// License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Affero General Public License for more details.
//
// You should have received a copy of the GNU Affero General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef PDE_COMMON_H
#define PDE_COMMON_H

#include<cstring>
#include<cmath>
#include<cstdlib>
#include<cassert>
#include<vector>
#include<mpi.h>

// A PDE integration scheme can be decomposed into
// multiple "substeps".  Each "substep" has a number of inputs and
// a number of outputs.  The outputs at each grid are computed from
// the inputs at the same grid and its nearest neighbors.
// 
// For example, gradient computation on a compact stencil is such a "substep".
// The inputs are the value of field variables; the outputs are the gradient
// of these variables computed with a compact stencil.  Another example of a
// "substep" is flux computation and accumulation in finite volume.  
// The inputs are the conservative variables and their gradients;
// the outputs are the divergence of fluxes, or time derivatives
// of the conservative variables.
//
// Every "substep" must use the same inputs as the outputs of
// the previous "substep".  This means some "substeps" must "forward"
// some variables not involved in computation, by setting some of
// its outputs to the values of its corresponding input.
// In the gradient computation "substep", for example, the output
// should include not only the gradients of the conservative variables,
// but also the conservative variables themselves.  This way, the outputs
// include all the necessary inputs of the next "substep", flux accumulation.
//

// ====== common interface for local schemes ======= //
// 
// Each "substep" is implemented as a function:
// 
// void substep(SpatialPoint<numIn, numOut>& point);
// 
// where numInput and numOutput are concrete numbers.
// For example, scalar gradient computation in 1D can be
// implemented as
//
// void scalarGradient(SpatialPoint<numIn, numOut>& point)
// {
//     point.outputs(0) = point.inputs(0); // "forwarding" the field value
//     double leftVal = point.nbr(0).inputs(0);
//     double rightVal = point.nbr(1).inputs(0);
//     point.outputs(1) = (rightVal - leftVal) / (2 * dx);
// }

template<size_t numInput, size_t numOutput>
class SpatialPoint {
    public:
    const double x;

    private:
    size_t iShift_;
    double * const & inputs_;
    double * const & outputs_;
    SpatialPoint const * const nbr_[2];

    public:
    SpatialPoint(double x, size_t iShift,
                 double *& inputs, double *& outputs,
                 const SpatialPoint* nbrL, const SpatialPoint* nbrR)
    : x(x), iShift_(iShift), inputs_(inputs), outputs_(outputs),
      nbr_{nbrL, nbrR} 
    {
    }

    const SpatialPoint& nbr(size_t i) {
        assert(i < 2);
        assert(nbr_[i] != 0);
        return *(nbr_[i]);
    }

    double inputs(size_t i) const {
        assert(i < numInput);
        return inputs_[iShift_ + i];
    }

    double & outputs(size_t i) {
        assert(i < numOutput);
        return outputs_[iShift_ + i];
    }
};

// ====== common interface of Discretization ======= //
// This code is here to illustrate the proper interface.
// It does not have to be inherited.

class Discretization {
    public:
    template<size_t numVar>
    Discretization(int numGrids, double dx,
         void (&localOperator)(SpatialPoint<0, numVar> local));

    template<size_t numInput, size_t numOutput>
    void applyOp(void (&localOperator)(SpatialPoint<numInput,numOutput>& point));
};

// ====== common utilities for pde_diamond and pde_classic ====== //

class DiscretizationBase {
    public:
    DiscretizationBase() {
        MPI_Init(0,0);
    }
    virtual ~DiscretizationBase() {
        MPI_Finalize();
    }
};

int iProc() {
    int iPr;
    MPI_Comm_rank(MPI_COMM_WORLD, &iPr);
    return iPr;
}

int iProcLeft() {
    int iProc, numProc;
    MPI_Comm_rank(MPI_COMM_WORLD, &iProc);
    MPI_Comm_size(MPI_COMM_WORLD, &numProc);
    return (iProc + numProc - 1) % numProc;
}

int iProcRight() {
    int iProc, numProc;
    MPI_Comm_rank(MPI_COMM_WORLD, &iProc);
    MPI_Comm_size(MPI_COMM_WORLD, &numProc);
    return (iProc + 1) % numProc;
}

const char* mpiRankString()
{
    static char mpiRankStr[256] = "";
    if (mpiRankStr[0] == 0) {
        int iProc, numProc;
        MPI_Comm_size(MPI_COMM_WORLD, &numProc);
        MPI_Comm_rank(MPI_COMM_WORLD, &iProc);

        sprintf(mpiRankStr, "%d", numProc);
        int maxStrLen = strlen(mpiRankStr);

        sprintf(mpiRankStr, "%d", iProc);
        int strLen = strlen(mpiRankStr);
        int padLen = maxStrLen - strLen;
        memmove(mpiRankStr + padLen, mpiRankStr, strLen + 1);
        memset(mpiRankStr, '0', padLen);
    }
    return mpiRankStr;
}

class ClassicSyncer1D {
    // Only used for syncing initial data
    private:
    MPI_Request reqs_[4];

    public:
    ClassicSyncer1D(double * pLeft, double * pRight,
                    double * pLeftGhost, double * pRightGhost, int numVar)
    {
        MPI_Isend(pLeft, numVar, MPI_DOUBLE,
                  iProcLeft(), 0, MPI_COMM_WORLD, reqs_);
        MPI_Isend(pRight, numVar, MPI_DOUBLE,
                  iProcRight(), 1, MPI_COMM_WORLD, reqs_ + 1);

        MPI_Irecv(pRightGhost, numVar, MPI_DOUBLE,
                  iProcRight(), 0, MPI_COMM_WORLD, reqs_ + 2);
        MPI_Irecv(pLeftGhost, numVar, MPI_DOUBLE,
                  iProcLeft(), 1, MPI_COMM_WORLD, reqs_ + 3);
    }

    void waitTillDone() {
        MPI_Status stats[4];
        MPI_Waitall(4, reqs_, stats);
    }

    ~ClassicSyncer1D() {
        waitTillDone();
    }
};

#endif
