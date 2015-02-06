#ifndef DIAMOND_SCHEME_H
#define DIAMOND_SCHEME_H

#include<vector>
#include<deque>
#include<cstring>
#include<cassert>
#include<cstdlib>
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
// . ` ` ` ` ` ` . . . . . .
// _ _ _ _ _ _ _ _ _ _ _ _ _   initial data

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
    // virtual void apply(const double* pInputs, const double* pInputsL,
    //                    const double* pInputsR, double* pOutputs,
    //                    const LocalMesh& mesh) const = 0;
    virtual void applyIterator(const double* pInputsBegin,
            const double* pInputsLBegin,
            const double* pInputsRBegin, double* pOutputsBegin,
            std::vector<LocalMesh>::const_iterator meshBegin,
            std::vector<LocalMesh>::const_iterator meshEnd) const = 0;
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
    void (&operator_)(
         const LocalInputs1D<numInput>& inputs,
         LocalOutputs1D<numOutput>& outputs,
         const LocalMesh& mesh);

    public:
    LocalOperator(void (&localOperator)(
                  const LocalInputs1D<numInput>& inputs,
                  LocalOutputs1D<numOutput>& outputs,
                  const LocalMesh& mesh))
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

    virtual void applyIterator(const double* pInputsBegin,
            const double* pInputsLBegin,
            const double* pInputsRBegin, double* pOutputsBegin,
            std::vector<LocalMesh>::const_iterator meshBegin,
            std::vector<LocalMesh>::const_iterator meshEnd) const
    {
        auto pInputs  = pInputsBegin;
        auto pInputsL = pInputsLBegin;
        auto pInputsR = pInputsRBegin;
        auto pOutputs = pOutputsBegin;
        for (auto mesh = meshBegin; mesh < meshEnd; ++mesh) {
            LocalInputs1D<numInput> inputs(pInputs, pInputsL, pInputsR);
            LocalOutputs1D<numOutput> outputs(pOutputs);
            operator_(inputs, outputs, *mesh);
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

    int isSendOrRecvComplete() {
        if (!hasSendOrRecvBeenCalled_) {
            return false;
        }
        else {
            int testResult;
            MPI_Test(&req_, &testResult, MPI_STATUS_IGNORE);
            return testResult;
        }
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

class Diamond1D {
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
    std::vector<LocalMesh> localMeshes_;
    std::vector<const LocalOperatorBase*> localOperators_;

    // This is where the actual computation is done
    int numVariables_;
    double * variablesData_;

    inline double * varAtGrid_(size_t iGrid) {
        return variablesData_ + iGrid * numVariables_;
    }

    LocalVariablesQueue * pLeftFoundation_, * pRightFoundation_;
    LocalVariablesQueue * pLeftRoof_, * pRightRoof_;

    void step_(const double* pInput, double* pOutput,
               size_t iFirst, size_t iLast,
               const LocalOperatorBase* pOp)
    {
        size_t numInput = pOp->numInputs();
        size_t numOutput = pOp->numOutputs();
        pOp->applyIterator(pInput + iFirst * numInput,
                pInput + (iFirst - 1) * numInput,
                pInput + (iFirst + 1) * numInput,
                pOutput + iFirst * numOutput,
                localMeshes_.begin() + iFirst - 1,
                localMeshes_.begin() + iLast);
    }

    void dequeueTwoGrids_(LocalVariablesQueue* pFoundation, double * pData)
    {
        size_t numVar = pFoundation->dequeue(pData);
        assert(numVariables_ == numVar);
        numVar = pFoundation->dequeue(pData + numVar);
        assert(numVariables_ == numVar);
    }

    bool isHalfDiamond_() {
        if (localMeshes_.size() == localOperators_.size() + 1) {
            return false;
        } else {
            assert(localMeshes_.size() == localOperators_.size() * 2);
            return true;
        }
    }

    void computeFirstHalf_() 
    {
        assert(!isHalfDiamond_());

        const size_t nGrid = localMeshes_.size();
        numVariables_ = localOperators_[0]->numInputs();
        variablesData_ = new double[numVariables_ * (nGrid + 2)];
        dequeueTwoGrids_(pLeftFoundation_, varAtGrid_(nGrid / 2 - 1));
        dequeueTwoGrids_(pRightFoundation_, varAtGrid_(nGrid / 2 + 1));

        size_t lastStep = nGrid/2 - 2;
        for (size_t iStep = 0; iStep <= lastStep; ++ iStep)
        {
            assert(numVariables_ == localOperators_[iStep]->numInputs());
            size_t numOutput = localOperators_[iStep]->numOutputs();
            double * pOutput = new double[numOutput * (nGrid + 2)];

            const size_t iFirst = nGrid/2 - iStep, iLast = nGrid/2 + iStep + 1;
            assert(iFirst > 0 && iLast <= nGrid);
            step_(variablesData_, pOutput, iFirst, iLast, localOperators_[iStep]);

            delete[] variablesData_;
            variablesData_ = pOutput;
            numVariables_ = numOutput;

            dequeueTwoGrids_(pLeftFoundation_, varAtGrid_(iFirst - 2));
            dequeueTwoGrids_(pRightFoundation_, varAtGrid_(iLast + 1));
        }
    }

    void enqueueTwoGrids(LocalVariablesQueue* pRoof, double * pData)
    {
        pRoof->enqueue(numVariables_, pData);
        pRoof->enqueue(numVariables_, pData + numVariables_);
    }

    void computeSecondHalf_()
    {
        const size_t nGrid = localMeshes_.size();

        size_t firstStep = nGrid/2 - 1, lastStep = nGrid - 2;
        if (isHalfDiamond_()) {
            lastStep -= firstStep;
            firstStep = 0;
        }

        for (size_t iStep = firstStep; iStep <= lastStep; ++ iStep)
        {
            assert(numVariables_ == localOperators_[iStep]->numInputs());
            size_t numOutput = localOperators_[iStep]->numOutputs();
            double * pOutput = new double[numOutput * (nGrid + 2)];

            const size_t iFirst = nGrid/2 - (lastStep - iStep),
                         iLast  = nGrid/2 + (lastStep - iStep) + 1;
            assert(iFirst > 0 && iLast <= nGrid);
            step_(variablesData_, pOutput, iFirst, iLast,
                  localOperators_[iStep]);

            delete[] variablesData_;
            variablesData_ = pOutput;
            numVariables_ = numOutput;

            enqueueTwoGrids(pLeftRoof_, varAtGrid_(iFirst));
            enqueueTwoGrids(pRightRoof_, varAtGrid_(iLast - 1));
        }
        delete[] variablesData_;
    }

    public:
    // constructor from initial condition
    Diamond1D(std::vector<LocalMesh>::const_iterator firstGrid,
              std::vector<LocalMesh>::const_iterator lastGrid,
              std::deque<const LocalOperatorBase*>::const_iterator firstOp,
              std::deque<const LocalOperatorBase*>::const_iterator lastOp,
              double* initialData, int numVar)
    :
        localMeshes_(firstGrid, lastGrid),
        localOperators_(firstOp, lastOp),
        pLeftRoof_(0), pRightRoof_(0)
    {
        // std::cout << "Creating Diamond with initial data" << std::endl;

        assert(localMeshes_.size() == localOperators_.size() * 2);
        numVariables_ = numVar;
        variablesData_ = new double[numVar * (localMeshes_.size() + 2)];

        size_t bytesToCopy = sizeof(double) * numVar * localMeshes_.size();
        memcpy(variablesData_ + numVar, initialData, bytesToCopy);

        ClassicSyncer1D sync(variablesData_, localMeshes_.size(), numVar);
        sync.waitTillDone();

        // This tells the "compute" function to only process the top half
        // of the diamond
        pLeftFoundation_ = pRightFoundation_ = 0;
    }

    private:
    size_t foundationBytes_()
    {
        if (isHalfDiamond_()) {
            return 0;
        } else {
            size_t foundationBytes = 0;
            for (size_t iOp = 0; iOp < localMeshes_.size() / 2; ++iOp) {
                size_t localVarBytes = sizeof(size_t)
                    + sizeof(double) * localOperators_[iOp]->numInputs();
                foundationBytes += localVarBytes * 2;
            }
            return foundationBytes;
        }
    }

    size_t roofBytes_()
    {
        size_t firstOp = localMeshes_.size() / 2 - 1;
        size_t lastOp = localMeshes_.size() - 2;

        if (isHalfDiamond_()) {
            lastOp -= firstOp;
            firstOp = 0;
        }

        size_t roofBytes = 0;
        for (size_t iOp = firstOp; iOp <= lastOp; ++iOp) {
            size_t localVarBytes = sizeof(size_t)
                + sizeof(double) * localOperators_[iOp]->numOutputs();
            roofBytes += localVarBytes * 2;
        }
        return roofBytes;
    }

    public:
    // constructor from other diamonds
    Diamond1D(std::vector<LocalMesh>::const_iterator firstGrid,
              std::vector<LocalMesh>::const_iterator lastGrid,
              std::deque<const LocalOperatorBase*>::const_iterator firstOp,
              std::deque<const LocalOperatorBase*>::const_iterator lastOp,
              int iProcLeftFoundationIsFrom, int tagLeftFoundationIsFrom,
              int iProcRightFoundationIsFrom, int tagRightFoundationIsFrom)
    :
        localMeshes_(firstGrid, lastGrid),
        localOperators_(firstOp, lastOp),
        pLeftRoof_(0), pRightRoof_(0)
    {
        // std::cout << "Creating full Diamond" << std::endl;

        assert(localMeshes_.size() % 2 == 0);
        assert(localMeshes_.size() == localOperators_.size() + 1);

        pLeftFoundation_ = new LocalVariablesQueue(foundationBytes_());
        pLeftFoundation_->Irecv(iProcLeftFoundationIsFrom,
                                tagLeftFoundationIsFrom);
        pRightFoundation_ = new LocalVariablesQueue(foundationBytes_());
        pRightFoundation_->Irecv(iProcRightFoundationIsFrom,
                                 tagRightFoundationIsFrom);
    }

    void compute(int iProcLeftRoofGoesTo, int tagLeftRoofGoesTo,
                 int iProcRightRoofGoesTo, int tagRightRoofGoesTo)
    {
        // std::cout << "Computing Diamond" << std::endl;

        if (pLeftFoundation_ && pRightFoundation_) {
            pLeftFoundation_->waitForSendOrRecv();
            pRightFoundation_->waitForSendOrRecv();
            computeFirstHalf_();
            delete pLeftFoundation_;
            delete pRightFoundation_;
        }

        pLeftRoof_ = new LocalVariablesQueue(roofBytes_());
        pRightRoof_ = new LocalVariablesQueue(roofBytes_());
        computeSecondHalf_();
        pLeftRoof_->Isend(iProcLeftRoofGoesTo, tagLeftRoofGoesTo);
        pRightRoof_->Isend(iProcRightRoofGoesTo, tagRightRoofGoesTo);

        // std::cout << "Finish computing Diamond" << std::endl;
    }

    bool isSendComplete() {
        if (pLeftRoof_ == 0 || pRightRoof_ == 0) {
            return false;
        }
        return pLeftRoof_->isSendOrRecvComplete()
            && pRightRoof_->isSendOrRecvComplete();
    }

    ~Diamond1D() {
        delete pLeftRoof_;
        delete pRightRoof_;
    }
};

class DiamondDiscretization1D {
    private:
    std::deque<const LocalOperatorBase*> localOperators_;
    std::vector<LocalMesh> localMeshes_;
    std::deque<Diamond1D> diamonds_;
    size_t initialNumVar_;
    double* initialData_;

    public:
    ~DiamondDiscretization1D() {
        while (localOperators_.size() > 0) {
            delete localOperators_.front();
            localOperators_.pop_front();
        }
        MPI_Finalize();
    }

    template<size_t numVar>
    DiamondDiscretization1D(int numGrids, double dx,
         void (&localOperator)(
               LocalOutputs1D<numVar>&, const LocalMesh&))
    {
        MPI_Init(0, 0);

        double x0 = numGrids * dx * iProc();

        assert(numGrids % 2 == 0);
        for (int iGrid = 0; iGrid < numGrids * 3 / 2; ++iGrid) {
            localMeshes_.emplace_back(x0 + iGrid * dx, dx);
        }

        initialNumVar_ = numVar;
        initialData_ = new double[numGrids * numVar];
        for (int iGrid = 0; iGrid < numGrids; ++iGrid) {
            LocalOutputs1D<numVar> localVar(initialData_ + iGrid * numVar);
            localOperator(localVar, localMeshes_[iGrid]);
        }
    }

    private:
    size_t numGrids_() {
        assert(localMeshes_.size() % 3 == 0);
        return localMeshes_.size() / 3 * 2;
    }

    bool lastDiamondOnTheLeft_;

    void initialDiamond_() {
        assert(diamonds_.empty());
        diamonds_.emplace_back(localMeshes_.begin(),
                               localMeshes_.begin() + numGrids_(),
                               localOperators_.begin(),
                               localOperators_.end(),
                               initialData_, initialNumVar_);

        const int tagLeftwards = 1, tagToSelf = 0;
        diamonds_.front().compute(iProcLeft(), tagLeftwards,
                                   iProc(), tagToSelf);

        lastDiamondOnTheLeft_ = true;
    }

    void newDiamond_() {
        const int tagLeftwards = 1, tagRightwards = 2, tagToSelf = 0;
        if (lastDiamondOnTheLeft_) {
            diamonds_.emplace_back(localMeshes_.begin(),
                                   localMeshes_.begin() + numGrids_(),
                                   localOperators_.begin(),
                                   localOperators_.end(),
                                   iProc(), tagToSelf,
                                   iProcRight(), tagLeftwards);
            lastDiamondOnTheLeft_ = false;
            diamonds_.back().compute(iProc(), tagToSelf,
                                   iProcRight(), tagRightwards);
        }
        else {
            diamonds_.emplace_back(localMeshes_.begin(),
                                   localMeshes_.begin() + numGrids_(),
                                   localOperators_.begin(),
                                   localOperators_.end(),
                                   iProcLeft(), tagRightwards,
                                   iProc(), tagToSelf);
            lastDiamondOnTheLeft_ = true;
            diamonds_.back().compute(iProcLeft(), tagLeftwards,
                                     iProc(), tagToSelf);
        }
    }

    public:
    template<size_t numInput, size_t numOutput>
    void applyOp(void (&localOperator)(
                 const LocalInputs1D<numInput>& inputs,
                 LocalOutputs1D<numOutput>& outputs,
                 const LocalMesh& mesh))
    {
        localOperators_.push_back(new LocalOperator<numInput, numOutput>(localOperator));
        
        if (initialData_ && localOperators_.size() == numGrids_() / 2)
        {
            initialDiamond_();

            localOperators_.pop_front();
            delete[] initialData_;
            initialData_ = 0;
        }

        if (localOperators_.size() == numGrids_() - 1)
        {
            newDiamond_();

            for (size_t iStep = 0; iStep < numGrids_() / 2; ++iStep) {
                delete localOperators_.front();
                localOperators_.pop_front();
            }

            while (diamonds_.size() && diamonds_.front().isSendComplete()) {
                diamonds_.pop_front();
            }
        }
    }
};

#endif
