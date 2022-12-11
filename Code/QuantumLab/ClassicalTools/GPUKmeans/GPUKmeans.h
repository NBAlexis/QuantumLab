//=============================================================================
// FILENAME : GPUKmeans.h
// 
// DESCRIPTION:
// This is for kmeans to anomally detection
// Before supporing quantum version, we impliment classical version first
//
// REVISION: [dd/mm/yy]
//  [09/12/2022 nbale]
//=============================================================================

#ifndef _GPUKMEANS_H_
#define _GPUKMEANS_H_

__BEGIN_NAMESPACE

class QLAPI QLGPUKmeans
{
public:

    /**
    * the number of data is decided by the file
    */
    QLGPUKmeans(const CCString& fileName, BYTE dim, BYTE k, UBOOL bDebugChange);
    ~QLGPUKmeans();

    void Build(UINT minChange);
    void Save(const CCString& fileName);

protected:

    UBOOL m_bDebugChange;
    BYTE m_byK;
    UINT m_uiN;
    BYTE m_byDim;
    BYTE m_byStride;

    UINT m_uiBlock;
    UINT m_uiThread;

    FLOAT* m_pDeviceData;
    BYTE* m_pDeviceK;
    BYTE* m_pHostK;
    TArray<UINT> m_HostCount;
    TArray<TArray<FLOAT>> m_pHostCenters;
    TArray<FLOAT*> m_pDeviceCenters;
    FLOAT** m_pDeviceDeviceCenters;

    FLOAT* m_pDeviceFLOATWorkingBuffer;
    UINT* m_pDeviceUINTWorkingBuffer;
    BYTE* m_pDeviceTempKBuffer;

    void InitialK();
    UBOOL CalculateCountAndSplit();
    void CalculateCenter();
    UINT ReClassify();

};

__END_NAMESPACE


#endif //#ifndef _GPUKMEANS_H_

//=============================================================================
// END OF FILE
//=============================================================================