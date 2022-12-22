//=============================================================================
// FILENAME : QRamSubroutine.h
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
// REVISION: [dd/mm/yy]
//  [20/12/2022 nbale]
//=============================================================================

#ifndef _QRAMSUBROUTINE_H_
#define _QRAMSUBROUTINE_H_

__BEGIN_NAMESPACE


class QLAPI QLQuantumKmeans
{
public:

    /**
    * pData is not freed in QLQuantumKmeans
    */
    QLQuantumKmeans(BYTE maxK, const Real* pData);
    ~QLQuantumKmeans();

    QLGate CompareCircuit(const Real* hostVectorBuffer);

protected:

    void CalculateDegrees(const Real* hostVectorBuffer);
    QLGate ApplyInitialState() const;
    void ApplyInverseVector(QLGate& gate) const;

    BYTE m_byMaxK;
    BYTE m_byQubit;
    BYTE m_byVectorCount;

    UINT m_uiBlock;
    UINT m_uiThread;

    Real* m_pDeviceVBuffer;
    QLComplex* m_pDeviceCVBuffer;

    Real* m_pDeviceY1Buffer;
    Real* m_pDeviceY2Buffer;
    Real* m_pDeviceZ1Buffer;
    Real* m_pDeviceZ2Buffer;

    Real* m_pHostY1Buffer;
    Real* m_pHostY2Buffer;
    Real* m_pHostZ1Buffer;
    Real* m_pHostZ2Buffer;
};

/**
* hostVectorBuffer is vectors v[i] + w, for example, 64 v[i] + 1 w.
* in this case, vectorCount = 65, qubitCount = 2 + log2(64) = 8
*/
extern QLGate QLAPI QKMeanCompareCircuit(const Real* hostVectorBuffer, UINT vectorCount, BYTE qubitCount);

__END_NAMESPACE


#endif //#ifndef _QRAMSUBROUTINE_H_

//=============================================================================
// END OF FILE
//=============================================================================