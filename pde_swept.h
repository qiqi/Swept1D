#ifndef PDE_SWEPT_H
#define PDE_SWEPT_H

#include<vector>
#include<cstring>
#include<cassert>
#include<cstdlib>
#include<algorithm>
#include<cmath>
#include<mpi.h>
#include"pde_common.h"
#include"PngWriter.hpp"

//       . .
//     . . . .
//   . . . . . .  
// + + . . . . # #
// + + + . . # # # #   
// + + + + # # # # # # 
// + + + ` ` # # # # . .
// + + ` ` ` ` # # . . . .
// . ` ` ` ` ` ` . . . . . .   initial data

// ====== implementation of DiamondDiscretization1D ======= //
//
class LocalOperatorBase {
    // I define this base class, so that a list of LocalOperatrs
    // with different numbers of inputs and outputs can be
    // stored in a queue, and called sequentially using the same
    // virtual function "apply".
    // The first two arguments of "apply" should be pointers to
    // LocalInputs1D<numInput> and LocalOutputs1D<numOutput>
    // instances.
    public:
    virtual void applyToArray(size_t nGrid, const double* pInputsBegin,
            double* pOutputsBegin, double x0, double dx) const = 0;
    virtual void applyToIterator(
            std::vector<SpatialPoint<0,0>>::iterator pBegin,
            std::vector<SpatialPoint<0,0>>::iterator pEnd) const = 0;
    virtual size_t numInputs() const = 0;
    virtual size_t numOutputs() const = 0;
    virtual ~LocalOperatorBase() {}
};

template<size_t numInput, size_t numOutput>
class LocalOperator : public LocalOperatorBase {
    public:
    virtual size_t numInputs() const { return numInput; }
    virtual size_t numOutputs() const { return numOutput; }
    private:
    void (&operator_)(SpatialPoint<numInput,numOutput>& point);

    public:
    LocalOperator(void (&localOperator)(SpatialPoint<numInput,numOutput>& point))
    : operator_(localOperator) {}

    virtual ~LocalOperator() {}

    // virtual void apply(const double* pInputs, const double* pInputsL,
    //                    const double* pInputsR, double* pOutputs,
    //                    const LocalMesh& mesh) const
    // {
    //     LocalInputs1D<numInput> inputs(pInputs, pInputsL, pInputsR);
    //     LocalOutputs1D<numOutput> outputs(pOutputs);
    //     operator_(inputs, outputs, mesh);
    // }

    virtual void applyToArray(size_t nGrid, const double* pInputsBegin,
                              double* pOutputsBegin, double x0, double dx) const
    {
        double (*pOutput)[numOutput] = (double (*)[numOutput])pOutputsBegin;
        const double (*pInput)[numInput] = (const double (*)[numInput])pInputsBegin;

        for (size_t iGrid = 0; iGrid < nGrid; ++iGrid) {
            double x = x0 + iGrid * dx;
            SpatialPoint<numInput, numOutput> p(x, pInput[iGrid], pOutput[iGrid]);
            SpatialPoint<numInput, numOutput> pL(x - dx, pInput[iGrid-1], 0);
            SpatialPoint<numInput, numOutput> pR(x + dx, pInput[iGrid+1], 0);
            p.addNeighbors(&pL, &pR);

            operator_(p);
        }
    }

    virtual void applyToIterator(
            std::vector<SpatialPoint<0,0>>::iterator pBegin,
            std::vector<SpatialPoint<0,0>>::iterator pEnd) const
    {
        for (auto it = pBegin; it < pEnd; ++it) {
            auto p = (SpatialPoint<numInput, numOutput>&) (*it);
            operator_(p);
        }
    }
};

class LocalVariablesQueue {
    private:
    size_t maxBytes_;
    char* pData_;
    size_t enqueueByte_, dequeueByte_;
    MPI_Request req_;
    int hasSendOrRecvBeenCalled_;

    public:
    LocalVariablesQueue(size_t maxBytes)
    : maxBytes_(maxBytes), enqueueByte_(0), dequeueByte_(0),
      hasSendOrRecvBeenCalled_(0)
    {
        pData_ = (char*)malloc(maxBytes);
    }

    virtual ~LocalVariablesQueue() {
        free(pData_);
    }

    void Irecv(int iProc, int tag) {
        assert(!hasSendOrRecvBeenCalled_);
        hasSendOrRecvBeenCalled_ = 1;
        MPI_Irecv(pData_, maxBytes_, MPI_BYTE, iProc, tag, MPI_COMM_WORLD, &req_);
        // std::cout << "Irecv()..." << req_ << " iProc:" << iProc << ",Tag:" << tag << std::endl;
    }

    void Isend(int iProc, int tag) {
        assert(!hasSendOrRecvBeenCalled_);
        hasSendOrRecvBeenCalled_ = 1;
        MPI_Isend(pData_, maxBytes_, MPI_BYTE, iProc, tag, MPI_COMM_WORLD, &req_);
        // std::cout << "Isend()..." << req_ << " iProc:" << iProc << ",Tag:" << tag << std::endl;
    }

    void waitForSendOrRecv() {
        assert(hasSendOrRecvBeenCalled_);
        // std::cout << "waitForSendOrRecv()..." << req_ << std::endl;
        MPI_Wait(&req_, MPI_STATUS_IGNORE);
        // std::cout << "waitForSendOrRecv() done." << std::endl;
    }

    void enqueue(size_t numVar, const double* pData) {
        *(size_t*)(pData_ + enqueueByte_) = numVar;
        memcpy(pData_ + enqueueByte_ + sizeof(size_t), pData,
               sizeof(double) * numVar);
        enqueueByte_ += sizeof(size_t) + sizeof(double) * numVar;
        assert(enqueueByte_ <= maxBytes_);
    }

    size_t dequeue(double* pData) {
        size_t numVar = *(size_t*)(pData_ + dequeueByte_);
        memcpy(pData, pData_ + dequeueByte_ + sizeof(size_t),
               sizeof(double) * numVar);
        dequeueByte_ += sizeof(size_t) + sizeof(double) * numVar;
        assert(dequeueByte_ <= maxBytes_);
        return numVar;
    }
};

class DiamondTop {
    // TODO
    // An instance of this class computes one space-time diamond
    // When constructing the instance, one feeds in a series of LocalMesh
    // instances, spaning the spatial domain of this diamond, as well
    // as a series of LocalOperator instances, spanning the time domain
    // of this diamond. From these inputs, we know how much input data to
    // expect from either the initial condition or the "foundation" of
    // the diamonds.  Computation of the whole diamond is performed
    // upon its construction. Afterwards, call "getRoof" function to obtain
    // the outputs.
    private:
    double x0_, dx_;
    std::vector<const LocalOperatorBase*> localOperators_;

    size_t numVariables_;
    double * variablesData_;

    double * varAtGrid_(size_t iGrid) {
        return variablesData_ + iGrid * numVariables_;
    }

    LocalVariablesQueue leftRoof_, rightRoof_;

    void enqueueTwoGrids_(LocalVariablesQueue& q, double * pData)
    {
        q.enqueue(numVariables_, pData);
        q.enqueue(numVariables_, pData + numVariables_);
    }

    size_t roofBytes_()
    {
        size_t queueBytes = 0;
        std::for_each(localOperators_.begin(), localOperators_.end(),
                      [&] (const LocalOperatorBase * op) {
            if (op == 0) return;
            size_t localVarBytes = sizeof(size_t)
                + sizeof(double) * op->numInputs();
            if (op == localOperators_.back()) {
                localVarBytes += sizeof(size_t)
                    + sizeof(double) * op->numOutputs();
            }
            queueBytes += localVarBytes * 2;
        });
        return queueBytes;
    }

    void compute_() {
        size_t nStep = localOperators_.size();
        size_t nGrid = 2 * (nStep + 1);

        for (size_t iStep = 0; iStep < nStep; ++ iStep)
        {
            auto op = localOperators_[iStep];
            if (op == 0) continue;

            enqueueTwoGrids_(leftRoof_, varAtGrid_(0));
            enqueueTwoGrids_(rightRoof_, varAtGrid_(nGrid - 2));

            size_t numInput = op->numInputs();
            size_t numOutput = op->numOutputs();
            assert(numVariables_ == numInput);

            nGrid -= 2;
            const double * pInput = variablesData_;
            double * pOutput = new double[numOutput * nGrid];

            double x0 = x0_ + (iStep + 1) * dx_;
            op->applyToArray(nGrid, pInput + numInput, pOutput, x0, dx_);

            delete[] variablesData_;
            variablesData_ = pOutput;
            numVariables_ = numOutput;
        }

        if (localOperators_.back()) {
            enqueueTwoGrids_(leftRoof_, varAtGrid_(0));
            enqueueTwoGrids_(rightRoof_, varAtGrid_(nGrid - 2));
            assert(nGrid == 2);
        }
    }

    public:
    DiamondTop(std::vector<const LocalOperatorBase*>::const_iterator opBegin,
               std::vector<const LocalOperatorBase*>::const_iterator opEnd,
               double x0, double dx, const double * data,
               int iProcLeftRoofGoesTo, int tagLeftRoofGoesTo,
               int iProcRightRoofGoesTo, int tagRightRoofGoesTo)
    :
        x0_(x0), dx_(dx), localOperators_(opBegin, opEnd),
        leftRoof_(roofBytes_()),
        rightRoof_(roofBytes_())
    {
        size_t nStep = localOperators_.size();
        size_t nGrid = 2 * (nStep + 1);
        numVariables_ = (*opBegin)->numInputs();
        size_t nDouble = numVariables_ * nGrid;

        variablesData_ = new double[nDouble];
        memcpy(variablesData_, data, sizeof(double) * nDouble);
        compute_();

        leftRoof_.Isend(iProcLeftRoofGoesTo, tagLeftRoofGoesTo);
        rightRoof_.Isend(iProcRightRoofGoesTo, tagRightRoofGoesTo);
        delete[] variablesData_;
    }

    virtual ~DiamondTop() {
        leftRoof_.waitForSendOrRecv();
        rightRoof_.waitForSendOrRecv();
    }
};

class DiamondBottom {
    // TODO
    // An instance of this class computes one space-time diamond
    // When constructing the instance, one feeds in a series of LocalMesh
    // instances, spaning the spatial domain of this diamond, as well
    // as a series of LocalOperator instances, spanning the time domain
    // of this diamond. From these inputs, we know how much input data to
    // expect from either the initial condition or the "foundation" of
    // the diamonds.  Computation of the whole diamond is performed
    // upon its construction. Afterwards, call "getRoof" function to obtain
    // the outputs.
    private:
    double x0_, dx_;
    std::vector<const LocalOperatorBase*> localOperators_;

    size_t numVariables_;
    double * variablesData_;

    double * varAtGrid_(size_t iGrid) {
        return variablesData_ + iGrid * numVariables_;
    }

    LocalVariablesQueue leftFoundation_, rightFoundation_;

    void dequeueTwoGrids_(LocalVariablesQueue& q, double * pData)
    {
        size_t numVar = q.dequeue(pData);
        assert(numVariables_ == numVar);
        numVar = q.dequeue(pData + numVariables_);
        assert(numVariables_ == numVar);
    }

    size_t foundationBytes_()
    {
        size_t queueBytes = 0;
        std::for_each(localOperators_.begin(), localOperators_.end(),
                      [&] (const LocalOperatorBase * op) {
            if (op == 0) return;
            size_t localVarBytes = sizeof(size_t)
                + sizeof(double) * op->numInputs();
            queueBytes += localVarBytes * 2;
        });
        return queueBytes;
    }

    void compute_() {
        size_t nStep = localOperators_.size();
        size_t nGrid = 4;

        for (size_t iStep = 0; iStep < nStep; ++ iStep)
        {
            auto op = localOperators_[iStep];
            if (op == 0) continue;

            dequeueTwoGrids_(leftFoundation_, varAtGrid_(0));
            dequeueTwoGrids_(rightFoundation_, varAtGrid_(nGrid - 2));

            size_t numInput = op->numInputs();
            size_t numOutput = op->numOutputs();
            assert(numVariables_ == numInput);

            if (iStep < nStep - 1) {
                nGrid += 2;
            } else {
                nGrid -= 2;
            }
            const double * pInput = variablesData_;
            double * pOutput = new double[numOutput * nGrid];

            double x0 = x0_ + (nStep - iStep - 1) * dx_;
            if (iStep < nStep - 1) {
                op->applyToArray(nGrid - 4, pInput, pOutput + numOutput * 2, x0, dx_);
            } else {
                op->applyToArray(nGrid - 4, pInput, pOutput, x0, dx_);
            }

            delete[] variablesData_;
            variablesData_ = pOutput;
            numVariables_ = numOutput;
        }
        assert(nGrid == nStep * 2 || localOperators_.back() == 0);
    }

    public:
    DiamondBottom(std::vector<const LocalOperatorBase*>::const_iterator opBegin,
                  std::vector<const LocalOperatorBase*>::const_iterator opEnd,
                  double x0, double dx,
                  int iProcLeftFoundationIsFrom, int tagLeftFoundationIsFrom,
                  int iProcRightFoundationIsFrom, int tagRightFoundationIsFrom)
    :
        x0_(x0), dx_(dx), localOperators_(opBegin, opEnd),
        leftFoundation_(foundationBytes_()),
        rightFoundation_(foundationBytes_())
    {
        leftFoundation_.Irecv(iProcLeftFoundationIsFrom,
                              tagLeftFoundationIsFrom);
        rightFoundation_.Irecv(iProcRightFoundationIsFrom,
                               tagRightFoundationIsFrom);
        leftFoundation_.waitForSendOrRecv();
        rightFoundation_.waitForSendOrRecv();

        numVariables_ = (*opBegin)->numInputs();
        variablesData_ = new double[4 * numVariables_];
        compute_();
    }

    const double * top() {
        return variablesData_;
    }

    ~DiamondBottom() {
        delete variablesData_;
    }
};

class SweptDiscretization1D {
    private:
    size_t nGrid_;
    double x0_, dx_;

    std::vector<const LocalOperatorBase*> localOperators_;

    DiamondTop * pDiamondTop_;
    DiamondBottom * pDiamondBottom_;

    size_t initialNumVar_;
    double * initialData_;

    typedef enum {
        LEFT_DIAMOND,
        RIGHT_DIAMOND
    } WorkingSide_;
    WorkingSide_ isWorkingOn_;

    void clearLocalOperators_() {
        std::for_each(localOperators_.begin(), localOperators_.end(),
                    [] (const LocalOperatorBase * op) { if (op) delete op; });
        localOperators_.clear();
    }

    public:
    ~SweptDiscretization1D() {
        while (localOperators_.size() > 0) {
            applyOpPointer_(0);
        }
        MPI_Finalize();
    }

    template<size_t numVar>
    SweptDiscretization1D(size_t numGrids, double dx,
         void (&localOperator)(SpatialPoint<0, numVar>&))
    :
        nGrid_(numGrids), dx_(dx),
        pDiamondTop_(0), pDiamondBottom_(0),
        initialNumVar_(numVar),
        isWorkingOn_(LEFT_DIAMOND)
    {
        MPI_Init(0, 0);
        x0_ = numGrids * dx * iProc();

        initialData_ = new double[numGrids * numVar];
        for (size_t iGrid = 0; iGrid < numGrids; ++iGrid) {
            SpatialPoint<0,numVar> p(x0_ + iGrid * dx, 0, initialData_ + iGrid * numVar);
            localOperator(p);
        }
    }

    private:
    DiamondTop * computeDiamondTop_() {
        // std::cout << iProc() << ": computing top\n";
        const double * pData;
        if (initialData_) {
            assert(pDiamondBottom_ == 0);
            pData = initialData_;
            assert(initialNumVar_ == localOperators_[0]->numInputs());
        } else {
            assert(initialData_ == 0);
            pData = pDiamondBottom_->top();
        }

        const int tagLeftwards = 1, tagRightwards = 2, tagToSelf = 0;
        DiamondTop * top;
        if (isWorkingOn_ == LEFT_DIAMOND) {
            top = new DiamondTop(
                    localOperators_.begin(), localOperators_.end(),
                    x0_, dx_, pData,
                    iProcLeft(), tagLeftwards, iProc(), tagToSelf);
            isWorkingOn_ = RIGHT_DIAMOND;
        } else {
            top = new DiamondTop(
                    localOperators_.begin(), localOperators_.end(),
                    x0_, dx_, pData,
                    iProc(), tagToSelf, iProcRight(), tagRightwards);
            isWorkingOn_ = LEFT_DIAMOND;
        }

        if (initialData_) {
            delete[] initialData_;
            initialData_ = 0;
        }
        return top;
    }

    DiamondBottom * computeDiamondBottom_() {
        // std::cout << iProc() << ": computing bottom\n";
        const int tagLeftwards = 1, tagRightwards = 2, tagToSelf = 0;
        if (isWorkingOn_ == LEFT_DIAMOND) {
            return new DiamondBottom(
                    localOperators_.begin(), localOperators_.end(),
                    x0_, dx_,
                    iProcLeft(), tagRightwards, iProc(), tagToSelf);
        } else {
            return new DiamondBottom(
                    localOperators_.begin(), localOperators_.end(),
                    x0_, dx_,
                    iProc(), tagToSelf, iProcRight(), tagLeftwards);
        }
    }

    void applyOpPointer_(LocalOperatorBase * p) {
        localOperators_.push_back(p);
        if (pDiamondTop_) {  // compute diamond bottom next
            if (localOperators_.size() == nGrid_ / 2) {
                pDiamondBottom_ = computeDiamondBottom_();
                delete pDiamondTop_;
                pDiamondTop_ = 0;
                clearLocalOperators_();
            }
        } else {
            if (localOperators_.size() == nGrid_ / 2 - 1) {
                pDiamondTop_ = computeDiamondTop_();
                if (pDiamondBottom_) {
                    delete pDiamondBottom_;
                }
            }
        }
    }

    public:
    template<size_t numInput, size_t numOutput>
    void applyOp(void (&localOperator)(SpatialPoint<numInput, numOutput>&))
    {
        auto p = new LocalOperator<numInput, numOutput>(localOperator);
        applyOpPointer_(p);
    }
};

#endif
