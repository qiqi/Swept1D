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
    virtual void applyToArray(void * pointsBegin, void * pointsEnd) const = 0;
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
    void (* const operator_)(SpatialPoint<numInput,numOutput>& point);

    public:
    LocalOperator(void (*localOperator)(SpatialPoint<numInput,numOutput>& point))
    : operator_(localOperator) {}

    virtual ~LocalOperator() {}

    virtual void applyToArray(void * pBegin, void * pEnd) const
    {
        SpatialPoint<numInput, numOutput>
            * pointBegin = pBegin, * pointEnd = pEnd;
        for (auto point = pointBegin; point < pointEnd; ++point) {
            operator_(*point);
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

class DiamondHalf {
    protected:
    size_t numPoints_;
    SpatialPoint<0,0> * spatialPoints_;

    size_t numVariables_;
    double * variablesData_, * inputs_, * outputs_;

    double * pInput_(size_t iPoint) {
        assert(iPoint <= numPoints_ - 1);
        return inputs_ + iPoint * 2 * numVariables_;
    }
    double * pOutput_(size_t iPoint) {
        assert(iPoint <= numPoints_ - 1);
        return outputs_ + iPoint * 2 * numVariables_;
    }

    void resizeWorkspace_(size_t newNumVariables)
    {
        assert (newNumVariables > numVariables_);
        double * newVariablesData
            = new double[2 * newNumVariables * numPoints_];
        double * newInputs = newVariablesData;
        double * newOutputs = newVariablesData + newNumVariables;

        if (numVariables_) {
            for (size_t i = 0; i < numPoints_; ++ i) {
                void * dest = newInputs + i * 2 * newNumVariables;
                void * src = inputs_ + i * 2 * numVariables_;
                memcpy(dest, src, sizeof(double) * numVariables_);

                dest = newOutputs + i * 2 * newNumVariables;
                src = outputs_ + i * 2 * numVariables_;
                memcpy(dest, src, sizeof(double) * numVariables_);
            }
            delete[] variablesData_;
        }

        numVariables_ = newNumVariables;
        variablesData_ = newVariablesData;
        inputs_ = newInputs;
        outputs_ = newOutputs;

        for (size_t i = 0; i < numPoints_; ++i) {
            auto p = spatialPoints_ + i;
            auto pL = (i == 0) ? 0 : p - 1;
            auto pR = (i == numPoints_ - 1) ? 0 : p + 1;
            size_t iShift = i * numVariables_ * 2;
            double x = p->x;
            p->~SpatialPoint();
            new(p) SpatialPoint<0,0>(x, iShift, inputs_, outputs_, pL, pR);
        }
    }

    template<typename LocalOpIter>
    void ensureSufficientMemory_(LocalOpIter opBegin, LocalOpIter opEnd)
    {
        size_t maxNumVar = 0;
        std::for_each(opBegin, opEnd, [&] (const LocalOperatorBase * op) {
            if (op == 0) return;
            if (op->numInputs() > maxNumVar) maxNumVar = op->numInputs();
            if (op->numOutputs() > maxNumVar) maxNumVar = op->numOutputs();
        });

        if (maxNumVar > numVariables_) {
            resizeWorkspace_(maxNumVar);
        }
    }

    void swapInputOutput_() {
        double * tmp = inputs_;
        inputs_ = outputs_;
        outputs_ = tmp;
    }

    public:

    DiamondHalf(size_t numPoints, double x0, double dx)
    :
        numPoints_(numPoints), numVariables_(0)
    {
        spatialPoints_ = (SpatialPoint<0,0>*)
            malloc(numPoints_ * sizeof(SpatialPoint<0,0>));
        for (size_t i = 0; i < numPoints_; ++i) {
            auto p = spatialPoints_ + i;
            new(p) SpatialPoint<0,0>(
                    x0 + i * dx, 0, inputs_, outputs_, nullptr, nullptr);
        }
    }

    virtual ~DiamondHalf() {
        for (size_t i = 0; i < numPoints_; ++i) {
            auto p = spatialPoints_ + i;
            p->~SpatialPoint<0,0>();
        }
        free(spatialPoints_);

        if (numVariables_) {
            delete[] variablesData_;
        }
    }
};

class DiamondTop : public DiamondHalf {
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
    int iProcLeftRoofGoesTo_, tagLeftRoofGoesTo_;
    int iProcRightRoofGoesTo_, tagRightRoofGoesTo_;

    public:
    DiamondTop(size_t numPoints, double x0, double dx,
               int iProcLeftRoofGoesTo, int tagLeftRoofGoesTo,
               int iProcRightRoofGoesTo, int tagRightRoofGoesTo)
    :
        DiamondHalf(numPoints, x0, dx),
        iProcLeftRoofGoesTo_(iProcLeftRoofGoesTo),
        tagLeftRoofGoesTo_(tagLeftRoofGoesTo),
        iProcRightRoofGoesTo_(iProcRightRoofGoesTo),
        tagRightRoofGoesTo_(tagRightRoofGoesTo)
    { }

    private:
    template<typename LocalOpIter>
    size_t roofBytes_(LocalOpIter opBegin, LocalOpIter opEnd)
    {
        size_t queueBytes = 0;
        std::for_each(opBegin, opEnd, [&] (const LocalOperatorBase * op) {
            if (op == 0) return;
            size_t localVarBytes = sizeof(size_t)
                + sizeof(double) * op->numInputs();
            queueBytes += localVarBytes * 2;
        });
        const LocalOperatorBase * op = *(opEnd - 1);
        if (op) {
            size_t localVarBytes = sizeof(size_t)
                + sizeof(double) * op->numOutputs();
            queueBytes += localVarBytes * 2;
        }
        return queueBytes;
    }

    void copyInitialData_(const double* pInitialData,
                          size_t dataSpacing, size_t numVar) {
        assert(numVar <= dataSpacing);
        assert(numVar <= numVariables_);
        for (size_t i = 0; i < numPoints_; ++i) {
            const double * src = pInitialData + i * dataSpacing;
            memcpy(pOutput_(i), src, sizeof(double) * numVar);
        }
    }

    public:
    template<typename LocalOpIter>
    void computeOps(LocalOpIter opBegin, LocalOpIter opEnd,
                    const double* pInitialData, size_t dataSpacing)
    {
        ensureSufficientMemory_(opBegin, opEnd);
        size_t roofBytes = roofBytes_(opBegin, opEnd);
        LocalVariablesQueue leftRoof(roofBytes), rightRoof(roofBytes);

        const LocalOperatorBase* opFirst = *opBegin;
        size_t numPrevOutput = opFirst->numInputs();
        copyInitialData_(pInitialData, dataSpacing, numPrevOutput);

        assert(numPoints_ == 2 * (opEnd - opBegin + 1));
        size_t iActivePointBegin = 0, iActivePointEnd = numPoints_;

        std::for_each(opBegin, opEnd, [&](const LocalOperatorBase* op)
        {
            if (op == 0) return;
            size_t numInput = op->numInputs();
            size_t numOutput = op->numOutputs();
            assert(numInput == numPrevOutput);

            leftRoof.enqueue(numInput, pInput_(iActivePointBegin));
            leftRoof.enqueue(numInput, pInput_(iActivePointBegin + 1));
            rightRoof.enqueue(numInput, pInput_(iActivePointEnd - 2));
            rightRoof.enqueue(numInput, pInput_(iActivePointEnd - 1));

            ++iActivePointBegin;
            --iActivePointEnd;

            swapInputOutput_();

            op->applyToArray(spatialPoints_ + iActivePointBegin,
                             spatialPoints_ + iActivePointEnd);
            numPrevOutput = numOutput;
        });

        const LocalOperatorBase * op = *(opEnd - 1);
        if (op) {
            assert(iActivePointEnd - iActivePointBegin == 2);
            leftRoof.enqueue(numPrevOutput, pInput_(iActivePointBegin));
            leftRoof.enqueue(numPrevOutput, pInput_(iActivePointBegin + 1));
            rightRoof.enqueue(numPrevOutput, pInput_(iActivePointEnd - 2));
            rightRoof.enqueue(numPrevOutput, pInput_(iActivePointEnd - 1));
        }

        leftRoof.Isend(iProcLeftRoofGoesTo_, tagLeftRoofGoesTo_);
        rightRoof.Isend(iProcRightRoofGoesTo_, tagRightRoofGoesTo_);
    }
};

class DiamondBottom : public DiamondHalf {
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
    int iProcLeftFoundationIsFrom_, tagLeftFoundationIsFrom_;
    int iProcRightFoundationIsFrom_, tagRightFoundationIsFrom_;

    public:
    DiamondBottom(size_t numPoints, double x0, double dx,
                  int iProcLeftFoundationIsFrom, int tagLeftFoundationIsFrom,
                  int iProcRightFoundationIsFrom, int tagRightFoundationIsFrom)
    :   // two more points than corresponding DiamondTop
        DiamondHalf(numPoints + 2, x0, dx),
        iProcLeftFoundationIsFrom_(iProcLeftFoundationIsFrom),
        tagLeftFoundationIsFrom_(tagLeftFoundationIsFrom),
        iProcRightFoundationIsFrom_(iProcRightFoundationIsFrom),
        tagRightFoundationIsFrom_(tagRightFoundationIsFrom)
    { } 

    private:
    template<typename LocalOpIter>
    size_t foundationBytes_(LocalOpIter opBegin, LocalOpIter opEnd)
    {
        size_t queueBytes = 0;
        std::for_each(opBegin, opEnd, [&] (const LocalOperatorBase * op) {
            if (op == 0) return;
            size_t localVarBytes = sizeof(size_t)
                + sizeof(double) * op->numInputs();
            queueBytes += localVarBytes * 2;
        });
        return queueBytes;
    }

    public:
    template<typename LocalOpIter>
    void computeOps(LocalOpIter opBegin, LocalOpIter opEnd)
    {
        ensureSufficientMemory_(opBegin, opEnd);
        size_t foundationBytes = foundationBytes_(opBegin, opEnd);
        LocalVariablesQueue leftFoundation(foundationBytes),
                            rightFoundation(foundationBytes);

        leftFoundation.Irecv(iProcLeftFoundationIsFrom_,
                             tagLeftFoundationIsFrom_);
        rightFoundation.Irecv(iProcRightFoundationIsFrom_,
                              tagRightFoundationIsFrom_);
        leftFoundation.waitForSendOrRecv();
        rightFoundation.waitForSendOrRecv();

        const LocalOperatorBase* opFirst = *opBegin;
        size_t numPrevOutput = opFirst->numInputs();

        assert(numPoints_ == 2 * (opEnd - opBegin + 1));
        size_t iActivePointBegin = numPoints_ / 2 - 1,
               iActivePointEnd = numPoints_ / 2 + 1;

        std::for_each(opBegin, opEnd, [&](const LocalOperatorBase* op)
        {
            if (op == 0) return;
            assert(op->numInputs() == numPrevOutput);
            size_t numOutput = op->numOutputs();

            size_t numDequeue;
            numDequeue = leftFoundation.dequeue(pInput_(iActivePointBegin - 1));
            numDequeue = leftFoundation.dequeue(pInput_(iActivePointBegin));
            numDequeue = rightFoundation.dequeue(pInput_(iActivePointEnd - 1));
            numDequeue = rightFoundation.dequeue(pInput_(iActivePointEnd));

            op->applyToArray(spatialPoints_ + iActivePointBegin,
                             spatialPoints_ + iActivePointEnd);

            --iActivePointBegin;
            ++iActivePointEnd;

            swapInputOutput_();
            numPrevOutput = numOutput;
        });
        const LocalOperatorBase * op = *(opEnd - 1);
        if (op) {
            assert(iActivePointBegin == 0);
            assert(iActivePointEnd == numPoints_);
        }
    }

    const double * finalData() {
        return inputs_;
    }

    size_t dataSpacing() {
        return numVariables_ * 2;
    }
};

const int tagLeftwards = 1, tagRightwards = 2, tagToSelf = 0;

class SweptDiscretization1D : public DiscretizationBase {
    private:
    size_t numPoints_;
    std::vector<const LocalOperatorBase*> localOperators_;

    DiamondTop diamondTopLeft_, diamondTopRight_;
    DiamondBottom diamondBottomLeft_, diamondBottomRight_;

    size_t initialNumVar_;
    double * initialData_;

    typedef enum {
        LEFT_DIAMOND,
        RIGHT_DIAMOND
        LEFT_DIAMOND,
        RIGHT_DIAMOND
    } WorkingDiamond_;
    WorkingDiamond_ isWorkingOn_;

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
    SweptDiscretization1D(size_t numPoints, double dx,
         void (&localOperator)(SpatialPoint<0, numVar>&))
    :
        DiscretizationBase(),
        numPoints_(numPoints),
        diamondTopLeft_(numPoints, numPoints * iProc() * dx, dx,
                        iProcLeft(), tagLeftwards, iProc(), tagToSelf),
        diamondTopRight_(numPoints, (numPoints / 2 + numPoints * iProc()) * dx, dx,
                         iProc(), tagToSelf, iProcRight(), tagRightwards),
        diamondBottomLeft_(numPoints, numPoints * dx * iProc(), dx,
                           iProcLeft(), tagRightwards, iProc(), tagToSelf),
        diamondBottomRight_(numPoints, (numPoints / 2 + numPoints * iProc()) * dx, dx,
                            iProc(), tagToSelf, iProcRight(), tagLeftwards),
        initialNumVar_(numVar),
        isWorkingOn_(LEFT_DIAMOND)
    {
        initialData_ = new double[numPoints * numVar];
        double x0 = numPoints * dx * iProc();
        for (size_t iPoint = 0; iPoint < numPoints; ++iPoint) {
            SpatialPoint<0,numVar> p(x0 + iPoint * dx, 0, initialData_ + iPoint * numVar);
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
