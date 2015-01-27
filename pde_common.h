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

#ifndef DIAMOND_H
#define DIAMOND_H

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
// void substep(LocalInputs1D<numInput>& inputs,
//              LocalOutputs1D<numOutput>& inputs,
//              LocalMesh& mesh);
// 
// where numInput and numOutput are concrete numbers.
// For example, scalar gradient computation in 1D can be
// implemented as
//
// void scalarGradient(LocalInputs1D<1>& inputs,
//                     LocalOutputs1D<2>& outputs,
//                     LocalMesh& mesh)
// {
//     outputs[0] = inputs[0]); // "forwarding" the field value
//     double leftVal = inputs[0].nbr(0);
//     double rightVal = inputs[0].nbr(1);
//     outputs[1] = (rightVal - leftVal) / (2 * mesh.dx));
// }

template<size_t numVar>
class LocalInputs1D {
    private:
    const double *pVal_, *pNbr_[2];

    public:
    class SingleVar_ {
        private:
            size_t iVar_;
            LocalInputs1D<numVar>* pInputs_;
        public:
        void assign(LocalInputs1D<numVar>* inputs, size_t iVar) {
            pInputs_ = inputs;
            iVar_ = iVar;
        }

        operator double() const {
            return pInputs_->pVal_[iVar_];
        }

        double nbr(size_t iNbr) const {
            assert(iNbr < 2);
            return pInputs_->pNbr_[iNbr][iVar_];
        }
    };

    private:
    SingleVar_ singleVars_[numVar];

    public:
    LocalInputs1D() : pVal_(0) {}
    LocalInputs1D(const double* pVal, const double* pValL, const double* pValR)
    {
        pVal_ = pVal;
        pNbr_[0] = pValL;
        pNbr_[1] = pValR;
        for (int iVar = 0; iVar < numVar; ++iVar) {
            singleVars_[iVar].assign(this, iVar);
        }
    }

    const SingleVar_& operator[] (size_t iVar) const {
        assert(iVar < numVar);
        return singleVars_[iVar];
    }
};

template<size_t numVar>
class LocalOutputs1D {
    private:
    double *pVal_;

    public:
    LocalOutputs1D(double* pVal) : pVal_(pVal) {}

    double& operator[](size_t iVar) {
        return pVal_[iVar];
    }
};

struct LocalMesh {
    double x, dx;
    LocalMesh(double x, double dx) : x(x), dx(dx) {}
    LocalMesh() : x(nan("")), dx(nan("")) {}
};

// ====== common interface of Discretization ======= //
// This code is here to illustrate the proper interface.
// It does not have to be inherited.

class Discretization {
    public:
    template<size_t numVar>
    Discretization(int numGrids, double dx,
         void (&localOperator)(
               LocalOutputs1D<numVar>&, const LocalMesh&));

    template<size_t numInput, size_t numOutput>
    void applyOp(void (&localOperator)(
                 const LocalInputs1D<numInput>& inputs,
                 LocalOutputs1D<numOutput>& outputs,
                 const LocalMesh& mesh));
};

// ====== common utilities for pde_diamond and pde_classic ====== //

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
    double *data_;
    size_t numGrids_, numVar_;
    MPI_Request reqs_[4];

    double * pGrid_(int iGrid) {
        return data_ + iGrid * numVar_;
    }

    public:
    ClassicSyncer1D(double * data, int numGrids, int numVar)
    : data_(data), numGrids_(numGrids), numVar_(numVar)
    {
        const int iGridLeft = 1, iGridRight = numGrids_;
        MPI_Isend(pGrid_(iGridLeft), numVar, MPI_DOUBLE,
                  iProcLeft(), 0, MPI_COMM_WORLD, reqs_);
        MPI_Isend(pGrid_(iGridRight), numVar, MPI_DOUBLE,
                  iProcRight(), 1, MPI_COMM_WORLD, reqs_ + 1);

        const int iGridLeftGhost = 0, iGridRightGhost = numGrids_ + 1;
        MPI_Irecv(pGrid_(iGridRightGhost), numVar, MPI_DOUBLE,
                  iProcRight(), 0, MPI_COMM_WORLD, reqs_ + 2);
        MPI_Irecv(pGrid_(iGridLeftGhost), numVar, MPI_DOUBLE,
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
