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
    QLQuantumKmeans(UINT maxK);
    ~QLQuantumKmeans();

    void Prepare(const CCString& fileName, const CCString& sStartCenterFile = _T(""), UINT n = 0);
    void KMeans(const CCString& sResultFileName, UINT uiStop, UINT uiRepeat, UBOOL bUseCC, UINT uiStep = 0, UBOOL bDebug = FALSE);

    void TestCircuit(const Real* hostVectors);

protected:

    void CalculateDegreesOnlyCenters();
    void CalculateDegrees();
    QLGate ApplyInitialState() const;
    void ApplyInverseVector(QLGate& gate, UBOOL bControlledCollpase) const;
    QLGate CompareCircuit(UINT uiIdx);
    QLGate CompareCircuit();
    UINT Measure(const QLGate& gate, UINT uiRepeat, UINT* count = NULL, UINT* measureCount = NULL, Real* measureProb = NULL) const;
    UINT MeasureWithoutCollapse(const QLGate& gate, UINT uiRepeat, UINT* count = NULL, UINT* measureCount = NULL) const;

    UBOOL CalculateCenters(UBOOL bDebug);
    UINT Reclassify(UINT uiIdx, UINT* measurecount);
    UINT Reclassify(UBOOL bDebug);
    void InitialK(UBOOL bDebug);
    void InitialWithCenterFile();
    void ExportDebugInfo();
    
    UINT m_byMaxK;
    BYTE m_byQubit;
    UINT m_byVectorCount;
    UINT m_uiRepeat;
    UINT m_bControlledCollapse;

    UINT m_uiN;

    UINT m_uiBlock;
    UINT m_uiThread;
    UINT m_uiBlockN;
    UINT m_uiThreadN;
    UINT m_uiBlockC;
    UINT m_uiThreadC;

    Real* m_pDeviceVBuffer;
    QLComplex* m_pDeviceCVBuffer;
    UINT* m_pDeviceKCounts; //How many points are there in one cluster
    UINT* m_pHostKCounts;
    UINT* m_pDeviceTempKBuffer; //used for split

    Real* m_pDeviceY1Buffer;
    Real* m_pDeviceY2Buffer;
    Real* m_pDeviceZ1Buffer;
    Real* m_pDeviceZ2Buffer;

    Real* m_pHostY1Buffer;
    Real* m_pHostY2Buffer;
    Real* m_pHostZ1Buffer;
    Real* m_pHostZ2Buffer;

    Real* m_pDeviceData;
    Real* m_pHostDataY1Buffer;
    Real* m_pHostDataY2Buffer;
    Real* m_pHostDataZ1Buffer;
    Real* m_pHostDataZ2Buffer;
    UINT* m_pHostKValues; //used for count how many point changes
    UINT* m_pDeviceKValues; //save the k-values
    
    Real* m_pDeviceRealWorkingBuffer; //used for reduce
    UINT* m_pDeviceUIntWorkingBuffer;

    BYTE* m_pHostVectorReal;
    BYTE* m_pHostVectorImag;
    LONGLONG m_llVeclen;

    //========= debug use ==========
    Real* m_pHostCenters;
    Real* m_pHostMeasureProbability;
    UINT m_uiStep;
    UINT* m_pMeasureCounts;
    CCString m_sSaveNameHead;

    //========= continue ==========
    CCString m_sStartCenterFile;
};


__END_NAMESPACE


#endif //#ifndef _QKMEANS_H_

//=============================================================================
// END OF FILE
//=============================================================================