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
            * pointBegin = (SpatialPoint<numInput, numOutput>*) pBegin,
            * pointEnd = (SpatialPoint<numInput, numOutput>*) pEnd;
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
    bool hasSendOrRecvBeenCalled_, hasSendOrRecvFinished_;

    public:
    LocalVariablesQueue(size_t maxBytes)
    : maxBytes_(maxBytes), enqueueByte_(0), dequeueByte_(0),
      hasSendOrRecvBeenCalled_(false), hasSendOrRecvFinished_(false)
    {
        pData_ = (char*)malloc(maxBytes);
    }

    virtual ~LocalVariablesQueue() {
        if (!hasSendOrRecvFinished_) {
            waitForSendOrRecv();
        }
        free(pData_);
    }

    void Irecv(int iProc, int tag) {
        assert(!hasSendOrRecvBeenCalled_);
        hasSendOrRecvBeenCalled_ = true;
        MPI_Irecv(pData_, maxBytes_, MPI_BYTE, iProc, tag, MPI_COMM_WORLD, &req_);
        // std::cout << "Irecv()..." << req_ << " iProc:" << iProc << ",Tag:" << tag << std::endl;
    }

    void Isend(int iProc, int tag) {
        assert(!hasSendOrRecvBeenCalled_);
        hasSendOrRecvBeenCalled_ = true;
        MPI_Isend(pData_, maxBytes_, MPI_BYTE, iProc, tag, MPI_COMM_WORLD, &req_);
        // std::cout << "Isend()..." << req_ << " iProc:" << iProc << ",Tag:" << tag << std::endl;
    }

    void waitForSendOrRecv() {
        assert(hasSendOrRecvBeenCalled_);
        assert(!hasSendOrRecvFinished_);
        // std::cout << "waitForSendOrRecv()..." << req_ << std::endl;
        MPI_Wait(&req_, MPI_STATUS_IGNORE);
        // std::cout << "waitForSendOrRecv() done." << std::endl;
        hasSendOrRecvFinished_ = true;
    }

    void enqueue(size_t numVar, const double* pData)
    {
        assert(!hasSendOrRecvBeenCalled_);

        *(size_t*)(pData_ + enqueueByte_) = numVar;
        memcpy(pData_ + enqueueByte_ + sizeof(size_t), pData,
               sizeof(double) * numVar);
        enqueueByte_ += sizeof(size_t) + sizeof(double) * numVar;
        assert(enqueueByte_ <= maxBytes_);
    }

    size_t dequeue(double* pData)
    {
        assert(hasSendOrRecvFinished_);
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

    LocalVariablesQueue * pLeftRoof, * pRightRoof;

    public:
    DiamondTop(size_t numPoints, double x0, double dx,
               int iProcLeftRoofGoesTo, int tagLeftRoofGoesTo,
               int iProcRightRoofGoesTo, int tagRightRoofGoesTo)
    :
        DiamondHalf(numPoints, x0, dx),
        iProcLeftRoofGoesTo_(iProcLeftRoofGoesTo),
        tagLeftRoofGoesTo_(tagLeftRoofGoesTo),
        iProcRightRoofGoesTo_(iProcRightRoofGoesTo),
        tagRightRoofGoesTo_(tagRightRoofGoesTo),
        pLeftRoof(nullptr), pRightRoof(nullptr)
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
        if (pLeftRoof) delete pLeftRoof;
        if (pRightRoof) delete pRightRoof;

        ensureSufficientMemory_(opBegin, opEnd);
        size_t roofBytes = roofBytes_(opBegin, opEnd);
        pLeftRoof = new LocalVariablesQueue(roofBytes);
        pRightRoof = new LocalVariablesQueue(roofBytes);

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

            pLeftRoof->enqueue(numInput, pInput_(iActivePointBegin));
            pLeftRoof->enqueue(numInput, pInput_(iActivePointBegin + 1));
            pRightRoof->enqueue(numInput, pInput_(iActivePointEnd - 2));
            pRightRoof->enqueue(numInput, pInput_(iActivePointEnd - 1));

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
            pLeftRoof->enqueue(numPrevOutput, pInput_(iActivePointBegin));
            pLeftRoof->enqueue(numPrevOutput, pInput_(iActivePointBegin + 1));
            pRightRoof->enqueue(numPrevOutput, pInput_(iActivePointEnd - 2));
            pRightRoof->enqueue(numPrevOutput, pInput_(iActivePointEnd - 1));
        }

        pLeftRoof->Isend(iProcLeftRoofGoesTo_, tagLeftRoofGoesTo_);
        pRightRoof->Isend(iProcRightRoofGoesTo_, tagRightRoofGoesTo_);
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
            size_t numOutput = op->numOutputs();
            assert(op->numInputs() == numPrevOutput);

            size_t numDequeue;
            numDequeue = leftFoundation.dequeue(pInput_(iActivePointBegin - 1));
            assert(numDequeue == op->numInputs());
            numDequeue = leftFoundation.dequeue(pInput_(iActivePointBegin));
            assert(numDequeue == op->numInputs());
            numDequeue = rightFoundation.dequeue(pInput_(iActivePointEnd - 1));
            assert(numDequeue == op->numInputs());
            numDequeue = rightFoundation.dequeue(pInput_(iActivePointEnd));
            assert(numDequeue == op->numInputs());

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

class SweptDiscretization1D : public DiscretizationBase {
    private:
    enum {
        TAG_Leftwards = 101,
        TAG_Rightwards = 102,
        TAG_ToSelf = 100
    };

    size_t numPoints_;
    std::vector<const LocalOperatorBase*> localOperators_;

    DiamondTop diamondTopLeft_, diamondTopRight_;
    DiamondBottom diamondBottomLeft_, diamondBottomRight_;

    size_t initialNumVar_;
    double * initialData_;

    typedef enum {
        LEFT_DIAMOND_TOP,
        RIGHT_DIAMOND_TOP,
        LEFT_DIAMOND_BOTTOM,
        RIGHT_DIAMOND_BOTTOM
    } WorkingDiamond_;
    WorkingDiamond_ isWorkingOn_;

    public:
    ~SweptDiscretization1D() {
        while (localOperators_.size() > 0) {
            applyOpPointer_(0);
        }
    }

    template<size_t numVar>
    SweptDiscretization1D(size_t numPoints, double dx,
         void (&localOperator)(SpatialPoint<0, numVar>&))
    :
        DiscretizationBase(),
        numPoints_(numPoints),
        diamondTopLeft_(numPoints, numPoints * iProc() * dx, dx,
                        iProcLeft(), TAG_Leftwards, iProc(), TAG_ToSelf),
        diamondTopRight_(numPoints, (numPoints / 2 + numPoints * iProc()) * dx, dx,
                         iProc(), TAG_ToSelf, iProcRight(), TAG_Rightwards),
        diamondBottomLeft_(numPoints, numPoints * dx * iProc(), dx,
                           iProcLeft(), TAG_Rightwards, iProc(), TAG_ToSelf),
        diamondBottomRight_(numPoints, (numPoints / 2 + numPoints * iProc()) * dx, dx,
                            iProc(), TAG_ToSelf, iProcRight(), TAG_Leftwards),
        initialNumVar_(numVar),
        isWorkingOn_(LEFT_DIAMOND_TOP)
    {
        initialData_ = new double[numPoints * numVar];
        double x0 = numPoints * dx * iProc();
        for (size_t iPoint = 0; iPoint < numPoints; ++iPoint) {
            double x = x0 + iPoint * dx;
            size_t iShift = iPoint * numVar;
            SpatialPoint<0,numVar> p(x, iShift, initialData_, initialData_,
                                     nullptr, nullptr);
            localOperator(p);
        }
    }

    private:
    void computeFirstDiamondTop_() {
        assert(isWorkingOn_ == LEFT_DIAMOND_TOP);
        assert(initialData_);
        assert(initialNumVar_ == localOperators_[0]->numInputs());

        DiamondTop & top = diamondTopLeft_;

        top.computeOps(localOperators_.begin(), localOperators_.end(),
                       initialData_, initialNumVar_);

        delete[] initialData_;
        initialData_ = 0;
    }

    void computeDiamondTop_() {
        // std::cout << iProc() << ": computing top\n";
        assert(isWorkingOn_ == LEFT_DIAMOND_TOP ||
               isWorkingOn_ == RIGHT_DIAMOND_TOP);
        assert(initialData_ == 0);

        DiamondTop & top = (isWorkingOn_ == LEFT_DIAMOND_TOP) ?
                            diamondTopLeft_ : diamondTopRight_;
        DiamondBottom & bottom = (isWorkingOn_ == LEFT_DIAMOND_TOP) ?
                            diamondBottomLeft_ : diamondBottomRight_;

        top.computeOps(localOperators_.begin(), localOperators_.end(),
                bottom.finalData(), bottom.dataSpacing());
    }

    void computeDiamondBottom_() {
        // std::cout << iProc() << ": computing bottom\n";
        assert(isWorkingOn_ == LEFT_DIAMOND_BOTTOM ||
               isWorkingOn_ == RIGHT_DIAMOND_BOTTOM);

        DiamondBottom & bottom = (isWorkingOn_ == LEFT_DIAMOND_BOTTOM) ?
                            diamondBottomLeft_ : diamondBottomRight_;
        bottom.computeOps(localOperators_.begin(), localOperators_.end());
    }

    void clearLocalOperators_() {
        std::for_each(localOperators_.begin(), localOperators_.end(),
                    [] (const LocalOperatorBase * op) { if (op) delete op; });
        localOperators_.clear();
    }

    void applyOpPointer_(LocalOperatorBase * p) {
        localOperators_.push_back(p);
        if (initialData_ && localOperators_.size() == numPoints_ / 2 - 1) {
            computeFirstDiamondTop_();
            isWorkingOn_ = RIGHT_DIAMOND_BOTTOM;
            return;
        }
        switch (isWorkingOn_) {
            case LEFT_DIAMOND_TOP:
                if (localOperators_.size() == numPoints_ / 2 - 1) {
                    computeDiamondTop_();
                    isWorkingOn_ = RIGHT_DIAMOND_BOTTOM;
                }
                break;
            case RIGHT_DIAMOND_TOP:
                if (localOperators_.size() == numPoints_ / 2 - 1) {
                    computeDiamondTop_();
                    isWorkingOn_ = LEFT_DIAMOND_BOTTOM;
                }
                break;
            case LEFT_DIAMOND_BOTTOM:
                if (localOperators_.size() == numPoints_ / 2) {
                    computeDiamondBottom_();
                    isWorkingOn_ = LEFT_DIAMOND_TOP;
                    clearLocalOperators_();
                }
                break;
            case RIGHT_DIAMOND_BOTTOM:
                if (localOperators_.size() == numPoints_ / 2) {
                    computeDiamondBottom_();
                    isWorkingOn_ = RIGHT_DIAMOND_TOP;
                    clearLocalOperators_();
                }
                break;
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
