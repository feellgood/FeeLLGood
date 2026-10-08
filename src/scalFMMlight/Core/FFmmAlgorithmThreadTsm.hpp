// See LICENCE file at project root
#ifndef FFMMALGORITHMTHREADTSM_HPP
#define FFMMALGORITHMTHREADTSM_HPP

#include "../Utils/FAssert.hpp"
#include "../Utils/FGlobal.hpp"
#include "../Containers/FOctree.hpp"

#include <omp.h>

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FFmmAlgorithmThreadTsm
* @brief
* Please read the license
*/

enum FFmmOperations
    {
    FFmmP2P  = (1 << 0),
    FFmmP2M  = (1 << 1),
    FFmmM2M  = (1 << 2),
    FFmmM2L  = (1 << 3),
    FFmmL2L  = (1 << 4),
    FFmmL2P  = (1 << 5),
    FFmmNearField = FFmmP2P,
    FFmmFarField  = (FFmmP2M|FFmmM2M|FFmmM2L|FFmmL2L|FFmmL2P),
    FFmmNearAndFarFields = (FFmmNearField|FFmmFarField)
    };


/** This class is a threaded FMM algorithm: it iterates on a tree and call the kernels with good arguments.
* It used the inspector-executor model :
* iterates on the tree and builds an array to work in parallel on this array
* !Warning: this class does not deallocate pointer given in arguments.
*
* Because this is a Target source model you do not need the P2P to be safe.
* You should not write on sources in the P2P method!
*/
template<class OctreeClass, class CellClass, class ContainerClass, class KernelClass, class LeafClass>
class FFmmAlgorithmThreadTsm
{
    /** The octree to work on */
    OctreeClass* const tree;

    /** The kernels */
    KernelClass** kernels;

    /** array of iterators */
    typename OctreeClass::Iterator* iterArray;

    /** maximum number of threads */
    const int MaxThreads;

    /** height of the octree */
    const int OctreeHeight;

protected:
     /** Where to start the work */
    int upperWorkingLevel;

    /** Where to end the work (exclusive) */
    int lowerWorkingLevel;

    /** Height of the tree */
    int nbLevelsInTree;

public:
    /** The constructor need the octree and the kernels used for computation
      * @param inTree the octree to work on
      * @param inKernels the kernels to call
      * An assert is launched if one of the arguments is null
      */
    FFmmAlgorithmThreadTsm(OctreeClass* const inTree, KernelClass* const inKernels)
                      : tree(inTree), kernels(nullptr), iterArray(nullptr),
                      MaxThreads(omp_get_max_threads()), OctreeHeight(tree->getHeight()),
                      upperWorkingLevel(2), lowerWorkingLevel(0), nbLevelsInTree(-1)
        {
        FAssertLF(tree, "tree cannot be null");
        this->kernels = new KernelClass*[MaxThreads];
        #pragma omp parallel num_threads(MaxThreads)  //create the openMP parallel region with MaxThreads threads
        {
            #pragma omp critical (InitFFmmAlgorithmTsm) //code to execute in one thread at a time, InitFFmmAlgorithmTsm is an optional name
            {
                this->kernels[omp_get_thread_num()] = new KernelClass(*inKernels);
            }
        }
        nbLevelsInTree = tree->getHeight();
        lowerWorkingLevel = nbLevelsInTree;
        }

    /** destructor */
    ~FFmmAlgorithmThreadTsm()
        {
        for(int idxThread = 0 ; idxThread < MaxThreads ; ++idxThread) { delete this->kernels[idxThread]; }
        delete [] this->kernels;
        }

    /** \brief Execute the whole fmm. */
    void execute()
        {
        upperWorkingLevel = 2;
        lowerWorkingLevel = nbLevelsInTree;
        FAssertLF(upperWorkingLevel <= lowerWorkingLevel);
        executeCore(FFmmNearAndFarFields);
        }

protected:
    /** To execute the fmm algorithm, call this function to run the complete algorithm */
    void executeCore(const unsigned operationsToProceed)
        {
        // Count leaf
        int numberOfLeafs = 0;
        typename OctreeClass::Iterator octreeIterator(tree);
        octreeIterator.gotoBottomLeft();
        do{
            ++numberOfLeafs;
        } while(octreeIterator.moveRight());
        iterArray = new typename OctreeClass::Iterator[numberOfLeafs];
        FAssertLF(iterArray, "iterArray bad alloc");

        if(operationsToProceed & FFmmP2M) bottomPass();

        if(operationsToProceed & FFmmM2M) upwardPass();

        if(operationsToProceed & FFmmM2L) transferPass();

        if(operationsToProceed & FFmmL2L) downardPass();

        if((operationsToProceed & FFmmP2P) || (operationsToProceed & FFmmL2P)) directPass((operationsToProceed & FFmmP2P),(operationsToProceed & FFmmL2P));

        delete [] iterArray;
        iterArray = nullptr;
        }

    /** P2M */
    void bottomPass(){
        typename OctreeClass::Iterator octreeIterator(tree);
        int numberOfLeafs = 0;
        // Iterate on leafs
        octreeIterator.gotoBottomLeft();
        do{
            iterArray[numberOfLeafs] = octreeIterator;
            ++numberOfLeafs;
        } while(octreeIterator.moveRight());

        const int chunkSize = FMath::Max(1 , numberOfLeafs/(omp_get_max_threads()*omp_get_max_threads())); // if not specified to schedule(), default chunksize is 1
        #pragma omp parallel num_threads(MaxThreads) //create the openMP parallel region with MaxThreads threads
        {
            KernelClass * const myThreadkernels = kernels[omp_get_thread_num()];
            #pragma omp for nowait schedule(dynamic,chunkSize)
            //dynamic means during the runtime the iterations are redistributed among threads, in contrast to static
            //there is a higher overhead than the static scheduling, but the later must be balanced to be efficient
            // auto could be used to let the scheduling decision to the compiler
            for(int idxLeafs = 0 ; idxLeafs < numberOfLeafs ; ++idxLeafs){
                // We need the current cell that represent the leaf and the list of particles
                ContainerClass* const sources = iterArray[idxLeafs].getCurrentListSrc();
                if(sources->getNbParticles()){
                    iterArray[idxLeafs].getCurrentCell()->setSrcChildTrue();
                    myThreadkernels->P2M( iterArray[idxLeafs].getCurrentCell() , sources);
                }
                if(iterArray[idxLeafs].getCurrentListTargets()->getNbParticles()){
                    iterArray[idxLeafs].getCurrentCell()->setTargetsChildTrue();
                }
            }
        }
    }

    /** M2M */
    void upwardPass(){
        // Start from leal level - 1
        typename OctreeClass::Iterator octreeIterator(tree);
        octreeIterator.gotoBottomLeft();
        octreeIterator.moveUp();

        for(int idxLevel = OctreeHeight - 2 ; idxLevel > lowerWorkingLevel-1 ; --idxLevel)
            { octreeIterator.moveUp(); }

        typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);

        // for each levels
        for(int idxLevel = FMath::Min(OctreeHeight - 2, lowerWorkingLevel - 1) ; idxLevel >= upperWorkingLevel ; --idxLevel ) 
            {
            int numberOfCells = 0;
            // for each cells
            do{
                iterArray[numberOfCells] = octreeIterator;
                ++numberOfCells;
            } while(octreeIterator.moveRight());
            avoidGotoLeftIterator.moveUp();
            octreeIterator = avoidGotoLeftIterator;// equal octreeIterator.moveUp(); octreeIterator.gotoLeft();

            const int chunkSize = FMath::Max(1 , numberOfCells/(omp_get_max_threads()*omp_get_max_threads()));
            #pragma omp parallel num_threads(MaxThreads)  //create the openMP parallel region with MaxThreads threads
            {
                KernelClass * const myThreadkernels = kernels[omp_get_thread_num()];
                #pragma omp for nowait schedule(dynamic, chunkSize)
                for(int idxCell = 0 ; idxCell < numberOfCells ; ++idxCell){
                    // We need the current cell and the child
                    // child is an array (of 8 child) that may be null
                    CellClass* potentialChild[8];
                    CellClass** const realChild = iterArray[idxCell].getCurrentChild();
                    CellClass* const currentCell = iterArray[idxCell].getCurrentCell();
                    int nbChildWithSrc = 0;
                    for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
                        potentialChild[idxChild] = nullptr;
                        if(realChild[idxChild]){
                            if(realChild[idxChild]->hasSrcChild()){
                                nbChildWithSrc += 1;
                                potentialChild[idxChild] = realChild[idxChild];
                            }
                            if(realChild[idxChild]->hasTargetsChild()){
                                currentCell->setTargetsChildTrue();
                            }
                        }
                    }
                    if(nbChildWithSrc){
                        currentCell->setSrcChildTrue();
                        myThreadkernels->M2M( currentCell , potentialChild, idxLevel);
                    }
                }
            }
        }
    }

    /** M2L */
    void transferPass(){
            typename OctreeClass::Iterator octreeIterator(tree);
            octreeIterator.moveDown();

            for(int idxLevel = 2 ; idxLevel < upperWorkingLevel ; ++idxLevel)
                { octreeIterator.moveDown(); }

            typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);

            // for each levels
            for(int idxLevel = upperWorkingLevel ; idxLevel < lowerWorkingLevel ; ++idxLevel )
                {
                int numberOfCells(0);
                // for each cells
                do{
                    iterArray[numberOfCells] = octreeIterator;
                    ++numberOfCells;
                } while(octreeIterator.moveRight());
                avoidGotoLeftIterator.moveDown();
                octreeIterator = avoidGotoLeftIterator;

                const int chunkSize = FMath::Max(1 , numberOfCells/(omp_get_max_threads()*omp_get_max_threads()));
                #pragma omp parallel num_threads(MaxThreads)  //create the openMP parallel region with MaxThreads threads
                {
                    KernelClass * const myThreadkernels = kernels[omp_get_thread_num()];
                    const CellClass* neighbors[342];
                    int neighborPositions[342];

                    #pragma omp for nowait schedule(dynamic,chunkSize)
                    for(int idxCell = 0 ; idxCell < numberOfCells ; ++idxCell){
                        CellClass* const currentCell = iterArray[idxCell].getCurrentCell();
                        if(currentCell->hasTargetsChild()){
                            const int counter = tree->getInteractionNeighbors(neighbors, neighborPositions, iterArray[idxCell].getCurrentGlobalCoordinate(), idxLevel);
                            if( counter ){
                                int counterWithSrc = 0;
                                for(int idxRealNeighbors = 0 ; idxRealNeighbors < counter ; ++idxRealNeighbors ){
                                    if(neighbors[idxRealNeighbors]->hasSrcChild()){
                                        neighbors[counterWithSrc] = neighbors[idxRealNeighbors];
                                        neighborPositions[counterWithSrc] = neighborPositions[idxRealNeighbors];
                                        ++counterWithSrc;
                                    }
                                }
                                if(counterWithSrc){
                                    myThreadkernels->M2L( currentCell , neighbors, neighborPositions, counterWithSrc, idxLevel);
                                }
                            }
                        }
                    }
                }
            }
        }

        /** L2L downward pass */
        void downardPass(){
            typename OctreeClass::Iterator octreeIterator(tree);
            octreeIterator.moveDown();

            for(int idxLevel = 2 ; idxLevel < upperWorkingLevel ; ++idxLevel)
                { octreeIterator.moveDown();}

            typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);

            const int heightMinusOne = lowerWorkingLevel - 1;
            // for each levels excepted leaf level
            for(int idxLevel = upperWorkingLevel ; idxLevel < heightMinusOne ; ++idxLevel)
                {
                int numberOfCells = 0;
                // for each cells
                do{
                    iterArray[numberOfCells] = octreeIterator;
                    ++numberOfCells;
                } while(octreeIterator.moveRight());
                avoidGotoLeftIterator.moveDown();
                octreeIterator = avoidGotoLeftIterator;

                const int chunkSize = FMath::Max(1 , numberOfCells/(omp_get_max_threads()*omp_get_max_threads()));
                #pragma omp parallel num_threads(MaxThreads)  //create the openMP parallel region with MaxThreads threads
                {
                    KernelClass * const myThreadkernels = kernels[omp_get_thread_num()];
                    #pragma omp for nowait schedule(dynamic,chunkSize)
                    for(int idxCell = 0 ; idxCell < numberOfCells ; ++idxCell){
                        if( iterArray[idxCell].getCurrentCell()->hasTargetsChild() ){
                            CellClass* potentialChild[8];
                            CellClass** const realChild = iterArray[idxCell].getCurrentChild();
                            CellClass* const currentCell = iterArray[idxCell].getCurrentCell();
                            for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
                                if(realChild[idxChild] && realChild[idxChild]->hasTargetsChild()){
                                    potentialChild[idxChild] = realChild[idxChild];
                                }
                                else{
                                    potentialChild[idxChild] = nullptr;
                                }
                            }
                            myThreadkernels->L2L( currentCell , potentialChild, idxLevel);
                        }
                    }
                }
            }
        }


    /** P2P */
    void directPass(const bool p2pEnabled, const bool l2pEnabled){
        int numberOfLeafs = 0;
        {
            typename OctreeClass::Iterator octreeIterator(tree);
            octreeIterator.gotoBottomLeft();
            // for each leaf
            do
                {
                iterArray[numberOfLeafs] = octreeIterator;
                ++numberOfLeafs;
                } while(octreeIterator.moveRight());
        }

        const int chunkSize = FMath::Max(1 , numberOfLeafs/(omp_get_max_threads()*omp_get_max_threads()));
        const int heightMinusOne = OctreeHeight - 1;
        #pragma omp parallel num_threads(MaxThreads)  //create the openMP parallel region with MaxThreads threads
        {
            KernelClass * const myThreadkernels = kernels[omp_get_thread_num()];
            // There is a maximum of 26 neighbors
            ContainerClass* neighbors[26];
            int neighborPositions[26];

            #pragma omp for nowait schedule(dynamic,chunkSize)
            for(int idxLeafs = 0 ; idxLeafs < numberOfLeafs ; ++idxLeafs){
                if( iterArray[idxLeafs].getCurrentCell()->hasTargetsChild() ){
                    if(l2pEnabled){
                        myThreadkernels->L2P(iterArray[idxLeafs].getCurrentCell(), iterArray[idxLeafs].getCurrentListTargets());
                    }
                    if(p2pEnabled){
                        // need the current particles and neighbors particles
                        if(iterArray[idxLeafs].getCurrentCell()->hasSrcChild()){
                            myThreadkernels->P2P( iterArray[idxLeafs].getCurrentGlobalCoordinate(), iterArray[idxLeafs].getCurrentListTargets(),
                                          iterArray[idxLeafs].getCurrentListSrc() , neighbors, neighborPositions, 0);
                        }
                        const int counter = tree->getLeafsNeighbors(neighbors, neighborPositions, iterArray[idxLeafs].getCurrentGlobalCoordinate(),heightMinusOne);
                        myThreadkernels->P2PRemote( iterArray[idxLeafs].getCurrentGlobalCoordinate(), iterArray[idxLeafs].getCurrentListTargets(),
                                      iterArray[idxLeafs].getCurrentListSrc() , neighbors, neighborPositions, counter);
                    }
                }
            }
        }
    }

};

#endif //FFMMALGORITHMTHREADTSM_HPP
