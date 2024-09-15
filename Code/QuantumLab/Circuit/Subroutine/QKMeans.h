//=============================================================================
// FILENAME : QKMeans.h
// 
// DESCRIPTION:
// This build the state for |n>|v_n>, with an arbitary phase, this is used for QKMeans
// 
// 1. build complex vector (v1, v2, v3, v3) for (y1, y2, y3, y4, y5, y6)
// 2. build state |n>|v_n>
// 3. act inverse v on qubits 0 and 1
// 4. controlled measurement on 0 and 1 to get |00>, now we have <phi_n|psi>|n>
// 5. now we have <phi_n|psi>|n>, do the amplitude amplitfy
// 6. Then apply measurement to get k-value
// 
// This QRam only work with QKMeans, only work with len=4 vectors.
// 
// To build QRam, just use amplitude-encode. list the vectors as [v1, v2, v3, ...].flatten(), 
// the result of amplitude-encode is automatically |n>|v_n>.
//
// REVISION: [dd/mm/yy]
//  [20/12/2022 nbale]
//=============================================================================

#ifndef _QKMEANS_H_
#define _QKMEANS_H_

__BEGIN_NAMESPACE


class QLAPI QLQuantumKmeans
{
public:

    /**
    * pData is not freed in QLQuantumKmeans
    */
    QLQuantumKmeans();
    ~QLQuantumKmeans();

    /**
    * Build <v|u_i>|i>
    *
    *
    */
    static void TestCircuitBuildState(const CCString& sReferenceCSV, const CCString& sAmplitudeSave, const CCString& sMeasureRate, UINT vectorCount, UINT testRepeat, UINT aaRepeat);

    static TArray<QLComplex> TestCircuitBuildState(UINT vectorDim, UINT vectorCount, UINT testRepeat, UINT aaRepeat);

    /**
    * test the circuit to build sum_i <vi|u>|i>, with amplitude amplification
    */
    static void TestCircuitBuildStateOnce(const QLMatrix& hostVi, const QLMatrix& hostU, UINT vectorCount, UINT vectorDim);

    //static void Kmeans2D(const CCString& sPoints, const CCString& sSaveK, 
    //    const CCString& sSaveCenter, const CCString& sRepeat,
    //    BYTE k, UINT iteration, UINT uiMinHit);

    //static void KNN2D(
    //    const CCString& sTrainingPoints, const CCString& sTestPoints,
    //    const CCString& sLoadK, const CCString& sSaveK, const CCString& sRepeat,
    //    UINT kHit, UINT uiMaxCluster);

    //static void KNN3D(
    //    const CCString& sTrainingPoints, const CCString& sTestPoints,
    //    const CCString& sLoadK, const CCString& sSaveK, const CCString& sRepeat,
    //    UINT kHit, UINT uiMaxCluster);

    //static void KNN2DAnsatz(const CCString& sAnsatz, const CCString& sTestPoints, const CCString& sSaveK, const CCString& sRepeat, UINT kHit);

    static void KNNAnsatz(const CCString& sAnsatz, const CCString& sTestPoints, const CCString& sScore, const CCString& sHit,
        BYTE ansatzQubits, BYTE byMeasureQubits, UBOOL bAdaptive, UINT uiAnsatzLevel, ELinkStyle eAnsatzStyle, ESingleLayer eAnsatzSingleLayer, ELinkLayer eAnsatzLayer, UINT uiRepeat);

    static void KNNAnsatzSE(const CCString& sAnsatz, const CCString& sTestPoints, const CCString& sScore, const CCString& sHit,
        BYTE ansatzQubits, BYTE byMeasureQubits, BYTE byEncodeQubit, UBOOL bAdaptive, UINT uiAnsatzLevel, ELinkStyle eAnsatzStyle, ESingleLayer eAnsatzSingleLayer, ELinkLayer eAnsatzLayer, ELinkStyle eSimpleEncodeStyle, ELinkLayer eSimpleencodeLayer, UINT uiRepeat);

    static void KNNAE(const CCString& sTrainingPoints, const CCString& sTestPoints, const CCString& sScore, BYTE byMeasureQubits, UINT uiRepeat);

    static void KNNSE(const CCString& sTrainingPoints, const CCString& sTestPoints, const CCString& sScore, BYTE byMeasureQubits, BYTE byEncodeQubits, ELinkStyle eStyle, ELinkLayer eLayer, UINT uiRepeat);

    //static void QAnomaly2D(const CCString& sReferenceCSV, const CCString& sPointCSV, const CCString& sBuildRate, 
    //    Real minX, Real maxX, Real minY, Real maxY);

    //static void QAnomaly3D(const CCString& sTrainingPoints, const CCString& sTestPoints, const CCString& sSaveScore);

protected:

    /**
    * 
    */
    static void LoadFile(const CCString& sReferenceCSV, Real** targetAbs, Real** targetPhase, UINT& uiDim, UINT& uiCount);

    

    UINT m_uiVectorDim;
    UINT m_uiVectorCount;
    UINT m_uiK;

    QLMatrix m_pOrignalPoints;


    Real* m_pDevicePointAbs;
    Real* m_pDevicePointPhase;

    Real* m_pDeviceCenterAbs;
    Real* m_pDeviceCenterPhase;

    Real* m_pWorkingSpaceAbsBuffer;
    Real* m_pWorkingSpacePhaseBuffer;
    
};



__END_NAMESPACE


#endif //#ifndef _QKMEANS_H_

//=============================================================================
// END OF FILE
//=============================================================================