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
    * We do not use extra CPU to judge, but the CSV is required that, w=Power 2 and each vector is already normalized.
    */
    void LoadFile(const CCString& sReferenceCSV, const CCString& sTestCSV);

    UBOOL ReferenceLoaded() const { return m_uiVectorDim > 0; }

    UBOOL TestLoaded() const { return m_uiTestVectorCount > 0; }

    /**
    * Build <v|u_i>|i>
    *
    *
    */
    static void TestCircuitBuildState(const CCString& sReferenceCSV, const CCString& sAmplitudeSave, const CCString& sMeasureRate, UINT vectorCount, UINT testRepeat);

    /**
    * Build <v|u_i>|i> once
    */
    static void TestCircuitBuildStateOnce(const QLMatrix& hostVi, const QLMatrix& hostU, UINT vectorCount, UINT vectorDim);

protected:

    /**
    * 
    */
    static void LoadFile(const CCString& sReferenceCSV, Real** targetAbs, Real** targetPhase, UINT& uiDim, UINT& uiCount);

    

    UINT m_uiVectorDim;
    UINT m_uiReferenceVectorCount;
    UINT m_uiTestVectorCount;

    Real* m_pDeviceReferenceVectorsAbs;
    Real* m_pDeviceReferenceVectorsPhase;
    Real* m_pDeviceTestVectorsAbs;
    Real* m_pDeviceTestVectorsPhase;

    Real* m_pWorkingSpaceAbsBuffer;
    Real* m_pWorkingSpacePhaseBuffer;
    
};



__END_NAMESPACE


#endif //#ifndef _QKMEANS_H_

//=============================================================================
// END OF FILE
//=============================================================================