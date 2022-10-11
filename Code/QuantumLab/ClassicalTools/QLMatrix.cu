//=============================================================================
// FILENAME : QLMatrix.cu
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [12/09/2022 nbale]
//=============================================================================

#include "QuantumLabPCH.h"

#include "cublas_v2.h"
#include "cusolverDn.h"

__BEGIN_NAMESPACE

#pragma region device functions

__inline__ __device__ QLComplex IExp(const QLComplex& c, Real r)
{
    return __cuCexpf(_make_cuComplex(-c.y * r, c.x * r));
}

__device__ complexfunc devicefunctionExp = __cuCexpf;
__device__ complexfunc devicefunctionSqrt = __cuCsqrtf;

__device__ complexfuncTwoR devicefunctionIExp = IExp;

#pragma endregion



#pragma region kernels

__global__ void _QL_LAUNCH_BOUND
_kernelRandomOne(QLComplex* pDevicePtr, UINT perthread, UINT uiMax)
{
    const UINT uiThread = threadIdx.x;
    for (UINT i = 0; i < perthread; ++i)
    {
        const UINT idx = uiThread * perthread + i;
        if (idx < uiMax)
        {
            pDevicePtr[idx] = _deviceRandomC(uiThread);
        }
    }
}

__global__ void _QL_LAUNCH_BOUND
_kernelTranspose(QLComplex* pDevicePtr, const QLComplex* __restrict__ pSource, UBOOL bConj, UINT uiXLen, UINT uiYLen, UINT uiMax)
{
    const UINT uiX = threadIdx.x + blockDim.x * blockIdx.x;
    const UINT uiY = threadIdx.y + blockDim.y * blockIdx.y;
    const UINT uiID1 = uiYLen * uiX + uiY;
    const UINT uiID2 = uiXLen * uiY + uiX;
    if (uiID1 < uiMax && uiID2 < uiMax)
    {
        pDevicePtr[uiID2] = bConj ? _cuConjf(pSource[uiID1]) : pSource[uiID1];
    }
}

__global__ void _QL_LAUNCH_BOUND
_kernelOpposite(QLComplex* pDevicePtr, UINT uiYLen, UINT uiMax)
{
    const UINT uiX = threadIdx.x + blockDim.x * blockIdx.x;
    const UINT uiY = threadIdx.y + blockDim.y * blockIdx.y;
    const UINT uiID1 = uiYLen * uiX + uiY;
    if (uiID1 < uiMax)
    {
        pDevicePtr[uiID1] = _make_cuComplex(-pDevicePtr[uiID1].x, -pDevicePtr[uiID1].y);
    }
}

__global__ void _QL_LAUNCH_BOUND
_kernelAdd(QLComplex* pDevicePtr, const QLComplex* __restrict__ pSource, UINT uiYLen, UINT uiMax)
{
    const UINT uiX = threadIdx.x + blockDim.x * blockIdx.x;
    const UINT uiY = threadIdx.y + blockDim.y * blockIdx.y;
    const UINT uiID1 = uiYLen * uiX + uiY;
    if (uiID1 < uiMax)
    {
        pDevicePtr[uiID1] = _cuCaddf(pDevicePtr[uiID1], pSource[uiID1]);
    }
}

__global__ void _QL_LAUNCH_BOUND
_kernelAdd(QLComplex* pDevicePtr, QLComplex v, UINT uiYLen, UINT uiMax)
{
    const UINT uiX = threadIdx.x + blockDim.x * blockIdx.x;
    const UINT uiID1 = uiYLen * uiX + uiX;
    if (uiID1 < uiMax)
    {
        pDevicePtr[uiID1] = _cuCaddf(pDevicePtr[uiID1], v);
    }
}

__global__ void _QL_LAUNCH_BOUND
_kernelAdd(QLComplex* pDevicePtr, Real v, UINT uiYLen, UINT uiMax)
{
    const UINT uiX = threadIdx.x + blockDim.x * blockIdx.x;
    const UINT uiID1 = uiYLen * uiX + uiX;
    if (uiID1 < uiMax)
    {
        pDevicePtr[uiID1] = cuCaddf_cr(pDevicePtr[uiID1], v);
    }
}

__global__ void _QL_LAUNCH_BOUND
_kernelMul(QLComplex* pDevicePtr, QLComplex v, UINT uiYLen, UINT uiMax)
{
    const UINT uiX = threadIdx.x + blockDim.x * blockIdx.x;
    const UINT uiY = threadIdx.y + blockDim.y * blockIdx.y;
    const UINT uiID1 = uiYLen * uiX + uiY;
    if (uiID1 < uiMax)
    {
        pDevicePtr[uiID1] = _cuCmulf(pDevicePtr[uiID1], v);
    }
}

__global__ void _QL_LAUNCH_BOUND
_kernelMul(QLComplex* pDevicePtr, Real v, UINT uiYLen, UINT uiMax)
{
    const UINT uiX = threadIdx.x + blockDim.x * blockIdx.x;
    const UINT uiY = threadIdx.y + blockDim.y * blockIdx.y;
    const UINT uiID1 = uiYLen * uiX + uiY;
    if (uiID1 < uiMax)
    {
        pDevicePtr[uiID1] = cuCmulf_cr(pDevicePtr[uiID1], v);
    }
}

__global__ void _QL_LAUNCH_BOUND
_kernelDiv(QLComplex* pDevicePtr, QLComplex v, UINT uiYLen, UINT uiMax)
{
    const UINT uiX = threadIdx.x + blockDim.x * blockIdx.x;
    const UINT uiY = threadIdx.y + blockDim.y * blockIdx.y;
    const UINT uiID1 = uiYLen * uiX + uiY;
    if (uiID1 < uiMax)
    {
        pDevicePtr[uiID1] = _cuCdivf(pDevicePtr[uiID1], v);
    }
}

__global__ void _QL_LAUNCH_BOUND
_kernelDiv(QLComplex* pDevicePtr, Real v, UINT uiYLen, UINT uiMax)
{
    const UINT uiX = threadIdx.x + blockDim.x * blockIdx.x;
    const UINT uiY = threadIdx.y + blockDim.y * blockIdx.y;
    const UINT uiID1 = uiYLen * uiX + uiY;
    if (uiID1 < uiMax)
    {
        pDevicePtr[uiID1] = cuCdivf_cr(pDevicePtr[uiID1], v);
    }
}

__global__ void _QL_LAUNCH_BOUND
_kernelDot(QLComplex* pDevicePtr, const QLComplex* __restrict__ pDevicePtr2, UINT uiYLen, UINT uiMax, UBOOL bConjLeft, UBOOL bConjRight)
{
    const UINT uiX = threadIdx.x + blockDim.x * blockIdx.x;
    const UINT uiY = threadIdx.y + blockDim.y * blockIdx.y;
    const UINT uiID1 = uiYLen * uiX + uiY;
    if (uiID1 < uiMax)
    {
        if (bConjLeft && !bConjRight)
        {
            pDevicePtr[uiID1] = _cuCmulf(_cuConjf(pDevicePtr[uiID1]), pDevicePtr2[uiID1]);
        }
        else if (!bConjLeft && bConjRight)
        {
            pDevicePtr[uiID1] = _cuCmulf(pDevicePtr[uiID1], _cuConjf(pDevicePtr2[uiID1]));
        }
        else if (bConjLeft && bConjRight)
        {
            pDevicePtr[uiID1] = _cuCmulf(_cuConjf(pDevicePtr[uiID1]), _cuConjf(pDevicePtr2[uiID1]));
        }
        else
        {
            pDevicePtr[uiID1] = _cuCmulf(pDevicePtr[uiID1], pDevicePtr2[uiID1]);
        }
    }
}

__global__ void _QL_LAUNCH_BOUND
_kernelReduceComplex(QLComplex* arr, UINT uiJump, UINT uiMax)
{
    //for length 16 array
    //for jump = 1, this is 1->0, 3->2, 5->4, 7->6, 9->10, 11->10, 13->12, 15->14 
    //for jump = 2, this is 2->0, 6->4, 10->8, 14->12 
    //for jump = 4, this is 4->0, 12->8 
    //for jump = 8, this is 8->0, and is finished.

    //id target = idx * (jump << 1)
    //id from = target + jump
    UINT uiIdFrom = (threadIdx.x + blockIdx.x * blockDim.x) * (uiJump << 1) + uiJump;
    if (uiIdFrom < uiMax)
    {
        arr[uiIdFrom - uiJump] = _cuCaddf(arr[uiIdFrom - uiJump], arr[uiIdFrom]);
    }
}

__global__ void _QL_LAUNCH_BOUND
_kernelComplexFunc(QLComplex* pDevicePtr, UINT uiYLen, UINT uiMax, complexfunc func)
{
    const UINT uiX = threadIdx.x + blockDim.x * blockIdx.x;
    const UINT uiY = threadIdx.y + blockDim.y * blockIdx.y;
    const UINT uiID1 = uiYLen * uiX + uiY;
    if (uiID1 < uiMax)
    {
        pDevicePtr[uiID1] = (*func)(pDevicePtr[uiID1]);
    }
}

__global__ void _QL_LAUNCH_BOUND
_kernelComplexFuncTwoR(QLComplex* pDevicePtr, UINT uiYLen, UINT uiMax, complexfuncTwoR func, Real r)
{
    const UINT uiX = threadIdx.x + blockDim.x * blockIdx.x;
    const UINT uiY = threadIdx.y + blockDim.y * blockIdx.y;
    const UINT uiID1 = uiYLen * uiX + uiY;
    if (uiID1 < uiMax)
    {
        pDevicePtr[uiID1] = (*func)(pDevicePtr[uiID1], r);
    }
}

__global__ void _QL_LAUNCH_BOUND
_kernelKronecker(QLComplex* pDevicePtr, 
    const QLComplex* __restrict__ pDevicem1, 
    const QLComplex* __restrict__ pDevicem2,
    UINT uiYm1,
    UINT uiYm2,
    UINT uiXm2,
    UINT uiMax1,
    UINT uiMax2)
{
    const UINT uiX1 = threadIdx.x + blockDim.x * blockIdx.x;
    const UINT uiY1 = threadIdx.y + blockDim.y * blockIdx.y;
    const UINT uiID1 = uiYm1 * uiX1 + uiY1;

    const UINT uiID2 = threadIdx.z + blockDim.z * blockIdx.z;
    const UINT uiX2 = uiID2 / uiYm2;
    const UINT uiY2 = uiID2 - uiX2 * uiYm2;

    if (uiID1 < uiMax1 && uiID2 < uiMax2)
    {
        const UINT uiIDres = ((uiX1 * uiXm2) + uiX2) * uiYm1 * uiYm2 + (uiY1 * uiYm2) + uiY2;
        pDevicePtr[uiIDres] = _cuCmulf(pDevicem1[uiID1], pDevicem2[uiID2]);
    }
}

#pragma endregion

QLMatrix::QLMatrix()
    : m_uiX(0)
    , m_uiY(0)
    , m_pData(NULL)
{

}

QLMatrix::QLMatrix(UINT uiX, UINT uiY)
    : m_uiX(uiX)
    , m_uiY(uiY)
{
    m_pData = new QLMatrixData();
    m_pData->m_nRefs = 1;
    m_pData->m_pData = (QLComplex*)calloc(m_uiX * m_uiY, sizeof(QLComplex));
}

QLMatrix::QLMatrix(UINT uiX, UINT uiY, QLComplex* buffer)
    : m_uiX(uiX)
    , m_uiY(uiY)
{
    m_pData = new QLMatrixData();
    m_pData->m_nRefs = 1; 
    m_pData->m_pData = buffer;
}

QLMatrix QLMatrix::CopyCreate(UINT uiX, UINT uiY, QLComplex* buffer)
{
    QLComplex* pCopyBuffer = reinterpret_cast<QLComplex*>(malloc(sizeof(QLComplex) * uiX * uiY));
    memcpy(pCopyBuffer, buffer, sizeof(QLComplex) * uiX * uiY);
    return QLMatrix(uiX, uiY, pCopyBuffer);
}

QLMatrix::QLMatrix(const QLMatrix& other)
    : m_uiX(other.m_uiX)
    , m_uiY(other.m_uiY)
    , m_pData(other.m_pData)
{
    m_pData->m_nRefs = m_pData->m_nRefs + 1;
}

QLMatrix::~QLMatrix()
{
    if (NULL != m_pData)
    {
        m_pData->m_nRefs = m_pData->m_nRefs - 1;

        if (m_pData->m_nRefs < 1)
        {
            appSafeDelete(m_pData);
        }
    }
}

void QLMatrix::RandomOne()
{
    UINT uiMax = m_uiX * m_uiY;
    UINT uiPerthread = Ceil(uiMax, _QL_LAUNCH_MAX_THREAD);
    QLComplex* deviceBuffer = NULL;

    checkCudaErrors(cudaMalloc((void**)&deviceBuffer, sizeof(QLComplex) * uiMax));

    _kernelRandomOne << <1, _QL_LAUNCH_MAX_THREAD >> > (deviceBuffer, uiPerthread, uiMax);

    OnChangeContent();
    checkCudaErrors(cudaMemcpy(m_pData->m_pData, deviceBuffer, sizeof(QLComplex) * uiMax, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaFree(deviceBuffer));
}

void QLMatrix::Print(const CCString& sName) const
{
    appPushLogDate(FALSE);
    if (!sName.IsEmpty())
    {
        appGeneral("%s=", sName);
    }
    appGeneral("{");
    for (UINT y = 0; y < m_uiY; ++y)
    {
        for (UINT x = 0; x < m_uiX; ++x)
        {
            if (0 == x)
            {
                appGeneral("{");
            }
            else
            {
                appGeneral(", ");
            }
            QLComplex toPrint = Get(x, y);
            appGeneral(appPrintComplex(toPrint.x, toPrint.y));
        }
        if (y == (m_uiY - 1))
        {
            appGeneral("}};\n");
        }
        else
        {
            appGeneral("},\n");
        }
    }
    appPopLogDate();
}

void QLMatrix::QR(QLMatrix& q, QLMatrix& r) const
{
    if (m_uiX > m_uiY)
    {
        QLMatrix sqm = SquareMatrixByAddOne();
        QLMatrix qq, rr;
        sqm._QR(qq, rr);
        q = qq.GetBlock(0, m_uiY, 0, m_uiY);
        r = rr.GetBlock(0, m_uiX, 0, m_uiY);
        return;
    }
    _QR(q, r);
}

//======================================================================================================
//from https://github.com/NVIDIA/CUDALibrarySamples/blob/master/cuSOLVER/orgqr/cusolver_orgqr_example.cu
//Does not support m_uiX > m_uiY
void QLMatrix::_QR(QLMatrix& q, QLMatrix& r) const
{
    cusolverDnHandle_t cusolverH = NULL;
    cublasHandle_t cublasH = NULL;
    cudaStream_t stream{};

    const INT m = m_uiY;
    const INT n = m_uiX;
    const INT k = std::min(m, n);

    QLComplex* pQHost = reinterpret_cast<QLComplex*>(malloc(sizeof(QLComplex) * m * k));
    QLComplex* pRHost = reinterpret_cast<QLComplex*>(calloc(k * n, sizeof(QLComplex)));
    //const std::vector<double> A = { 1.0, 4.0, 2.0, 2.0, 5.0, 1.0 };
    //std::vector<double> Q(lda * n, 0); // orthonormal columns
    //std::vector<double> R(n * n, 0);   // R = I - Q**T*Q

    /* device memory */
    QLComplex* d_A = NULL;
    QLComplex* d_tau = NULL;
    INT* d_info = NULL;
    QLComplex* d_work = NULL;

    INT lwork_geqrf = 0;
    INT lwork_orgqr = 0;
    INT lwork = 0;
    INT info = 0;

    //std::printf("A = (matlab base-1)\n");
    //print_matrix(m, n, A.data(), lda);
    //std::printf("=====\n");

    /* step 1: create cudense/cublas handle */
    checkCudaErrors(cusolverDnCreate(&cusolverH));
    checkCudaErrors(cublasCreate(&cublasH));

    checkCudaErrors(cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking));
    checkCudaErrors(cusolverDnSetStream(cusolverH, stream));
    checkCudaErrors(cublasSetStream(cublasH, stream));

    /* step 2: copy A and B to device */
    checkCudaErrors(cudaMalloc(reinterpret_cast<void**>(&d_A), sizeof(QLComplex) * m * n));
    checkCudaErrors(cudaMalloc(reinterpret_cast<void**>(&d_tau), sizeof(QLComplex) * n));
    checkCudaErrors(cudaMalloc(reinterpret_cast<void**>(&d_info), sizeof(INT)));

    checkCudaErrors(cudaMemcpyAsync(d_A, HostBuffer(), sizeof(QLComplex) * m * n, cudaMemcpyHostToDevice, stream));

    /* step 3: query working space of geqrf and orgqr */
#if _QL_DOUBLEFLOAT
    checkCudaErrors(cusolverDnZgeqrf_bufferSize(cusolverH, m, n, d_A, m, &lwork_geqrf));
    checkCudaErrors(cusolverDnZungqr_bufferSize(cusolverH, m, n, k, d_A, m, d_tau, &lwork_orgqr));
#else
    checkCudaErrors(cusolverDnCgeqrf_bufferSize(cusolverH, m, n, d_A, m, &lwork_geqrf));
    checkCudaErrors(cusolverDnCungqr_bufferSize(cusolverH, m, n, k, d_A, m, d_tau, &lwork_orgqr));
#endif

    lwork = std::max(lwork_geqrf, lwork_orgqr);

    checkCudaErrors(cudaMalloc(reinterpret_cast<void**>(&d_work), sizeof(QLComplex) * lwork));

    /* step 4: compute QR factorization */
#if _QL_DOUBLEFLOAT
    checkCudaErrors(cusolverDnZgeqrf(cusolverH, m, n, d_A, m, d_tau, d_work, lwork, d_info));
#else
    checkCudaErrors(cusolverDnCgeqrf(cusolverH, m, n, d_A, m, d_tau, d_work, lwork, d_info));
#endif

    /* check if QR is successful or not */
    checkCudaErrors(cudaMemcpyAsync(&info, d_info, sizeof(INT), cudaMemcpyDeviceToHost, stream));

    checkCudaErrors(cudaStreamSynchronize(stream));

    if (0 > info) 
    {
        appCrucial("%d-th parameter is wrong \n", -info);
    }

    for (INT i = 0; i < n; ++i)
    {
        checkCudaErrors(cudaMemcpyAsync(pRHost + i * n, d_A + i * m, sizeof(QLComplex) * std::min(i + 1, m), cudaMemcpyDeviceToHost, stream));
    }

    /* step 5: compute Q */
#if _QL_DOUBLEFLOAT
    checkCudaErrors(cusolverDnZungqr(cusolverH, m, n, k, d_A, m, d_tau, d_work, lwork, d_info));
#else
    checkCudaErrors(cusolverDnCungqr(cusolverH, m, n, k, d_A, m, d_tau, d_work, lwork, d_info));
#endif

    /* check if QR is good or not */
    checkCudaErrors(cudaMemcpyAsync(&info, d_info, sizeof(INT), cudaMemcpyDeviceToHost, stream));

    checkCudaErrors(cudaStreamSynchronize(stream));

    if (0 > info) 
    {
        appCrucial("%d-th parameter is wrong \n", -info);
    }

    for (INT i = 0; i < k; ++i)
    {
        checkCudaErrors(cudaMemcpyAsync(pQHost + i * m, d_A + i * m, sizeof(QLComplex) * m, cudaMemcpyDeviceToHost, stream));
    }

    checkCudaErrors(cudaStreamSynchronize(stream));

    //std::printf("Q = (matlab base-1)\n");
    //print_matrix(m, n, Q.data(), lda);


    q = QLMatrix(k, m, pQHost);
    r = QLMatrix(n, k, pRHost);
    
    checkCudaErrors(cudaStreamSynchronize(stream));

    /* free resources */
    checkCudaErrors(cudaFree(d_A));
    checkCudaErrors(cudaFree(d_tau));
    checkCudaErrors(cudaFree(d_info));
    checkCudaErrors(cudaFree(d_work));

    checkCudaErrors(cublasDestroy(cublasH));
    checkCudaErrors(cusolverDnDestroy(cusolverH));
    checkCudaErrors(cudaStreamDestroy(stream));
}

//======================================================================================================
//from https://github.com/NVIDIA/CUDALibrarySamples/blob/master/cuSOLVER/gesvd/cusolver_gesvd_example.cu
void QLMatrix::SVD(QLMatrix& u, QLMatrix& s, QLMatrix& v) const
{
    cusolverDnHandle_t cusolverH = NULL;
    cublasHandle_t cublasH = NULL;
    cudaStream_t stream = NULL;

    const INT m = m_uiY;   
    const INT n = m_uiX;   

    //const std::vector<double> A = { 1.0, 4.0, 2.0, 2.0, 5.0, 1.0 };
    //std::vector<double> U(lda * m, 0);  /* m-by-m unitary matrix, left singular vectors  */
    //std::vector<double> VT(lda * n, 0); /* n-by-n unitary matrix, right singular vectors */
    //std::vector<double> S(n, 0);        /* numerical singular value */
    //std::vector<double> S_exact = { 7.065283497082729,
    //                               1.040081297712078 }; /* exact singular values */

    INT info_gpu = 0;                                  /* host copy of error info */

    QLComplex* d_A = NULL;
    Real* d_S = NULL;  /* singular values */
    QLComplex* d_U = NULL;  /* left singular vectors */
    QLComplex* d_VT = NULL; /* right singular vectors */
    //QLComplex* d_W = NULL;  /* W = S*VT */

    INT* devInfo = NULL;

    INT lwork = 0; /* size of workspace */
    QLComplex* d_work = NULL;
    //QLComplex* d_rwork = NULL;

    //const double h_one = 1;
    //const double h_minus_one = -1;

    //std::printf("A = (matlab base-1)\n");
    //print_matrix(m, n, A.data(), lda);
    //std::printf("=====\n");

    /* step 1: create cusolver handle, bind a stream */
    checkCudaErrors(cusolverDnCreate(&cusolverH));
    checkCudaErrors(cublasCreate(&cublasH));

    checkCudaErrors(cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking));
    checkCudaErrors(cusolverDnSetStream(cusolverH, stream));

    /* step 2: copy A to device */
    checkCudaErrors(cudaMalloc(reinterpret_cast<void**>(&d_A), sizeof(QLComplex) * m * n));
    checkCudaErrors(cudaMalloc(reinterpret_cast<void**>(&d_S), sizeof(Real) * n));
    checkCudaErrors(cudaMalloc(reinterpret_cast<void**>(&d_U), sizeof(QLComplex) * m * m));
    checkCudaErrors(cudaMalloc(reinterpret_cast<void**>(&d_VT), sizeof(QLComplex) * m * n));
    checkCudaErrors(cudaMalloc(reinterpret_cast<void**>(&devInfo), sizeof(INT)));

    checkCudaErrors(
        cudaMemcpyAsync(d_A, HostBuffer(), sizeof(QLComplex) * m * n, cudaMemcpyHostToDevice, stream));

    /* step 3: query working space of SVD */
#if _QL_DOUBLEFLOAT
    checkCudaErrors(cusolverDnZgesvd_bufferSize(cusolverH, m, n, &lwork));
#else
    checkCudaErrors(cusolverDnCgesvd_bufferSize(cusolverH, m, n, &lwork));
#endif

    checkCudaErrors(cudaMalloc(reinterpret_cast<void**>(&d_work), sizeof(QLComplex) * lwork));

    /* step 4: compute SVD*/
    signed char jobu = 'A';  // all m columns of U
    signed char jobvt = 'A'; // all n columns of VT
#if _QL_DOUBLEFLOAT
    checkCudaErrors(cusolverDnZgesvd(cusolverH, jobu, jobvt, m, n, d_A, m, d_S, d_U,
        m, // ldu
        d_VT,
        m, // ldvt,
        d_work, lwork, NULL, devInfo));
#else
    checkCudaErrors(cusolverDnCgesvd(cusolverH, jobu, jobvt, m, n, d_A, m, d_S, d_U,
        m, // ldu
        d_VT,
        m, // ldvt,
        d_work, lwork, NULL, devInfo));
#endif

    QLComplex* pUHost = reinterpret_cast<QLComplex*>(malloc(sizeof(QLComplex) * m * m));
    Real* pSHost0 = reinterpret_cast<Real*>(malloc(sizeof(Real) * n));
    QLComplex* pSHost = reinterpret_cast<QLComplex*>(calloc(m * n, sizeof(QLComplex)));
    QLComplex* pVHost = reinterpret_cast<QLComplex*>(calloc(n * n, sizeof(QLComplex)));

    checkCudaErrors(cudaMemcpyAsync(pUHost, d_U, sizeof(QLComplex) * m * m, cudaMemcpyDeviceToHost, stream));
    //checkCudaErrors(cudaMemcpyAsync(pVHost, d_VT, sizeof(QLComplex) * m * n, cudaMemcpyDeviceToHost,stream));
    checkCudaErrors(cudaMemcpyAsync(pSHost0, d_S, sizeof(Real) * n, cudaMemcpyDeviceToHost, stream));
    checkCudaErrors(cudaMemcpyAsync(&info_gpu, devInfo, sizeof(INT), cudaMemcpyDeviceToHost, stream));

    for (INT i = 0; i < n; ++i)
    {
        pSHost[i * m + i] = _make_cuComplex(pSHost0[i], F(0.0));
        checkCudaErrors(cudaMemcpyAsync(pVHost + i * n, d_VT + i * m, sizeof(QLComplex) * n, cudaMemcpyDeviceToHost, stream));
    }
    

    checkCudaErrors(cudaStreamSynchronize(stream));

    if (0 > info_gpu) 
    {
        appCrucial("%d-th parameter is wrong \n", -info_gpu);
    }
    else if (0 < info_gpu) 
    {
        appWarning("WARNING: info = %d : gesvd does not converge \n", info_gpu);
    }

    u = QLMatrix(m, m, pUHost);
    s = QLMatrix(n, m, pSHost);
    v = QLMatrix(n, n, pVHost);

    appSafeFree(pSHost0);

    /* free resources */
    checkCudaErrors(cudaFree(d_A));
    checkCudaErrors(cudaFree(d_U));
    checkCudaErrors(cudaFree(d_VT));
    checkCudaErrors(cudaFree(d_S));
    checkCudaErrors(cudaFree(devInfo));
    checkCudaErrors(cudaFree(d_work));
    //checkCudaErrors(cudaFree(d_rwork));

    checkCudaErrors(cusolverDnDestroy(cusolverH));
    checkCudaErrors(cublasDestroy(cublasH));

    checkCudaErrors(cudaStreamDestroy(stream));
}

//======================================================================================================
//from https://github.com/NVIDIA/CUDALibrarySamples/blob/master/cuSOLVER/gesvd/cusolver_gesvd_example.cu
void QLMatrix::SVDJ(QLMatrix& u, QLMatrix& s, QLMatrix& v) const
{
    cusolverDnHandle_t cusolverH = NULL;
    cudaStream_t stream = NULL;
    gesvdjInfo_t gesvdj_params = NULL;

    const INT m = m_uiY;                   /* 1 <= m <= 32 */
    const INT n = m_uiX;                   /* 1 <= n <= 32 */
    //const int lda = m;                 /* lda >= m */
    //const int ldu = m;                 /* ldu >= m */
    //const int ldv = n;                 /* ldv >= n */
    const INT minmn = (m < n) ? m : n; /* min(m,n) */



    INT info = 0;                                      /* host copy of error info */

    QLComplex* d_A = NULL;
    Real* d_S = NULL; /* singular values */
    QLComplex* d_U = NULL; /* left singular vectors */
    QLComplex* d_V = NULL; /* right singular vectors */

    INT* d_info = NULL;

    INT lwork = 0;            /* size of workspace */
    QLComplex* d_work = NULL; /* device workspace for getrf */

    /* configuration of gesvdj  */
    const DOUBLE tol = 1.0E-7;
    const INT max_sweeps = 100;
    const cusolverEigMode_t jobz = CUSOLVER_EIG_MODE_VECTOR; // compute eigenvectors.
    const INT econ = 0;                                      /* econ = 1 for economy size */

    /* numerical results of gesvdj  */
    //Real residual = 0;
    //INT executed_sweeps = 0;

    //std::printf("m = %d, n = %d \n", m, n);
    //std::printf("tol = %E, default value is machine zero \n", tol);
    //std::printf("max. sweeps = %d, default value is 100\n", max_sweeps);
    //std::printf("econ = %d \n", econ);

    //std::printf("A = (matlab base-1)\n");
    //print_matrix(m, n, A.data(), lda);
    //std::printf("=====\n");

    /* step 1: create cusolver handle, bind a stream */
    checkCudaErrors(cusolverDnCreate(&cusolverH));

    checkCudaErrors(cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking));
    checkCudaErrors(cusolverDnSetStream(cusolverH, stream));

    /* step 2: configuration of gesvdj */
    checkCudaErrors(cusolverDnCreateGesvdjInfo(&gesvdj_params));

    /* default value of tolerance is machine zero */
    checkCudaErrors(cusolverDnXgesvdjSetTolerance(gesvdj_params, tol));

    /* default value of max. sweeps is 100 */
    checkCudaErrors(cusolverDnXgesvdjSetMaxSweeps(gesvdj_params, max_sweeps));

    /* step 3: copy A to device */
    checkCudaErrors(cudaMalloc(reinterpret_cast<void**>(&d_A), sizeof(QLComplex) * m * n));
    checkCudaErrors(cudaMalloc(reinterpret_cast<void**>(&d_S), sizeof(Real) * minmn));
    checkCudaErrors(cudaMalloc(reinterpret_cast<void**>(&d_U), sizeof(QLComplex) * m * m));
    checkCudaErrors(cudaMalloc(reinterpret_cast<void**>(&d_V), sizeof(QLComplex) * n * n));
    checkCudaErrors(cudaMalloc(reinterpret_cast<void**>(&d_info), sizeof(INT)));

    checkCudaErrors(cudaMemcpyAsync(d_A, HostBuffer(), sizeof(QLComplex) * m * n, cudaMemcpyHostToDevice, stream));

    /* step 4: query working space of SVD */
#if _QL_DOUBLEFLOAT
    checkCudaErrors(cusolverDnZgesvdj_bufferSize(
        cusolverH, jobz, /* CUSOLVER_EIG_MODE_NOVECTOR: compute singular values only */
                         /* CUSOLVER_EIG_MODE_VECTOR: compute singular value and singular vectors */
        econ,            /* econ = 1 for economy size */
        m,               /* nubmer of rows of A, 0 <= m */
        n,               /* number of columns of A, 0 <= n  */
        d_A,             /* m-by-n */
        m,               /* leading dimension of A */
        d_S,             /* min(m,n) */
                         /* the singular values in descending order */
        d_U,             /* m-by-m if econ = 0 */
                         /* m-by-min(m,n) if econ = 1 */
        m,               /* leading dimension of U, ldu >= max(1,m) */
        d_V,             /* n-by-n if econ = 0  */
                         /* n-by-min(m,n) if econ = 1  */
        n,               /* leading dimension of V, ldv >= max(1,n) */
        &lwork, gesvdj_params));
#else
    checkCudaErrors(cusolverDnCgesvdj_bufferSize(
        cusolverH, jobz, /* CUSOLVER_EIG_MODE_NOVECTOR: compute singular values only */
                         /* CUSOLVER_EIG_MODE_VECTOR: compute singular value and singular vectors */
        econ,            /* econ = 1 for economy size */
        m,               /* nubmer of rows of A, 0 <= m */
        n,               /* number of columns of A, 0 <= n  */
        d_A,             /* m-by-n */
        m,               /* leading dimension of A */
        d_S,             /* min(m,n) */
                         /* the singular values in descending order */
        d_U,             /* m-by-m if econ = 0 */
                         /* m-by-min(m,n) if econ = 1 */
        m,               /* leading dimension of U, ldu >= max(1,m) */
        d_V,             /* n-by-n if econ = 0  */
                         /* n-by-min(m,n) if econ = 1  */
        n,               /* leading dimension of V, ldv >= max(1,n) */
        &lwork, gesvdj_params));
#endif

    checkCudaErrors(cudaMalloc(reinterpret_cast<void**>(&d_work), sizeof(QLComplex) * lwork));

    /* step 5: compute SVD*/
#if _QL_DOUBLEFLOAT
    checkCudaErrors(cusolverDnZgesvdj(
        cusolverH, jobz, /* CUSOLVER_EIG_MODE_NOVECTOR: compute singular values only */
                         /* CUSOLVER_EIG_MODE_VECTOR: compute singular value and singular vectors */
        econ,            /* econ = 1 for economy size */
        m,               /* nubmer of rows of A, 0 <= m */
        n,               /* number of columns of A, 0 <= n  */
        d_A,             /* m-by-n */
        m,               /* leading dimension of A */
        d_S,             /* min(m,n)  */
                         /* the singular values in descending order */
        d_U,             /* m-by-m if econ = 0 */
                         /* m-by-min(m,n) if econ = 1 */
        m,               /* leading dimension of U, ldu >= max(1,m) */
        d_V,             /* n-by-n if econ = 0  */
                         /* n-by-min(m,n) if econ = 1  */
        n,               /* leading dimension of V, ldv >= max(1,n) */
        d_work, lwork, d_info, gesvdj_params));
#else
    checkCudaErrors(cusolverDnCgesvdj(
        cusolverH, jobz, /* CUSOLVER_EIG_MODE_NOVECTOR: compute singular values only */
                         /* CUSOLVER_EIG_MODE_VECTOR: compute singular value and singular vectors */
        econ,            /* econ = 1 for economy size */
        m,               /* nubmer of rows of A, 0 <= m */
        n,               /* number of columns of A, 0 <= n  */
        d_A,             /* m-by-n */
        m,               /* leading dimension of A */
        d_S,             /* min(m,n)  */
                         /* the singular values in descending order */
        d_U,             /* m-by-m if econ = 0 */
                         /* m-by-min(m,n) if econ = 1 */
        m,               /* leading dimension of U, ldu >= max(1,m) */
        d_V,             /* n-by-n if econ = 0  */
                         /* n-by-min(m,n) if econ = 1  */
        n,               /* leading dimension of V, ldv >= max(1,n) */
        d_work, lwork, d_info, gesvdj_params));
#endif

    QLComplex* pUHost = reinterpret_cast<QLComplex*>(malloc(sizeof(QLComplex) * m * m));
    Real* pSHost0 = reinterpret_cast<Real*>(malloc(sizeof(Real) * n));
    QLComplex* pSHost = reinterpret_cast<QLComplex*>(calloc(m * n, sizeof(QLComplex)));
    QLComplex* pVHost = reinterpret_cast<QLComplex*>(calloc(n * n, sizeof(QLComplex)));

    checkCudaErrors(cudaMemcpyAsync(pUHost, d_U, sizeof(QLComplex) * m * m, cudaMemcpyDeviceToHost, stream));
    checkCudaErrors(cudaMemcpyAsync(pVHost, d_V, sizeof(QLComplex) * n * n, cudaMemcpyDeviceToHost, stream));
    checkCudaErrors(cudaMemcpyAsync(pSHost0, d_S, sizeof(Real) * minmn, cudaMemcpyDeviceToHost, stream));
    checkCudaErrors(cudaMemcpyAsync(&info, d_info, sizeof(INT), cudaMemcpyDeviceToHost, stream));

    checkCudaErrors(cudaStreamSynchronize(stream));

    if (0 > info) 
    {
        appCrucial("%d-th parameter is wrong \n", -info);
    }
    else if (0 < info) 
    {
        appWarning("WARNING: info = %d : gesvdj does not converge \n", info);
    }

    for (INT i = 0; i < minmn; ++i)
    {
        pSHost[i * m + i] = _make_cuComplex(pSHost0[i], F(0.0));
    }
    

    u = QLMatrix(m, m, pUHost);
    s = QLMatrix(n, m, pSHost);
    v = QLMatrix(n, n, pVHost);

    //std::printf("S = singular values (matlab base-1)\n");
    //print_matrix(minmn, 1, S.data(), minmn);
    //std::printf("=====\n");

    //std::printf("U = left singular vectors (matlab base-1)\n");
    //print_matrix(m, m, U.data(), ldu);
    //std::printf("=====\n");

    //std::printf("V = right singular vectors (matlab base-1)\n");
    //print_matrix(n, n, V.data(), ldv);
    //std::printf("=====\n");

    ///* step 6: measure error of singular value */
    //double ds_sup = 0;
    //for (int j = 0; j < minmn; j++) {
    //    double err = fabs(S[j] - S_exact[j]);
    //    ds_sup = (ds_sup > err) ? ds_sup : err;
    //}
    //std::printf("|S - S_exact|_sup = %E \n", ds_sup);

    //CUSOLVER_CHECK(cusolverDnXgesvdjGetSweeps(cusolverH, gesvdj_params, &executed_sweeps));

    //CUSOLVER_CHECK(cusolverDnXgesvdjGetResidual(cusolverH, gesvdj_params, &residual));

    //std::printf("residual |A - U*S*V**H|_F = %E \n", residual);
    //std::printf("number of executed sweeps = %d \n", executed_sweeps);

    /* free resources */
    appSafeFree(pSHost0);

    checkCudaErrors(cudaFree(d_A));
    checkCudaErrors(cudaFree(d_U));
    checkCudaErrors(cudaFree(d_V));
    checkCudaErrors(cudaFree(d_S));
    checkCudaErrors(cudaFree(d_info));
    checkCudaErrors(cudaFree(d_work));

    checkCudaErrors(cusolverDnDestroyGesvdjInfo(gesvdj_params));

    checkCudaErrors(cusolverDnDestroy(cusolverH));

    checkCudaErrors(cudaStreamDestroy(stream));
}

void QLMatrix::RandomUnitary()
{
    RandomOne();
    QLMatrix q;
    QLMatrix r;

    QR(q, r);

    *this = q;
}

void QLMatrix::Transpose(UBOOL bConjugate)
{
    dim3 blocks;
    dim3 threads;
    GetDim(blocks, threads);
    QLComplex* deviceBuffer1 = NULL;
    QLComplex* deviceBuffer2 = NULL;
    checkCudaErrors(cudaMalloc((void**)&deviceBuffer1, sizeof(QLComplex) * m_uiX * m_uiY));
    checkCudaErrors(cudaMalloc((void**)&deviceBuffer2, sizeof(QLComplex) * m_uiX * m_uiY));
    checkCudaErrors(cudaMemcpy(deviceBuffer1, HostBuffer(), sizeof(QLComplex) * m_uiX * m_uiY, cudaMemcpyHostToDevice));
    _kernelTranspose << <blocks, threads >> > (deviceBuffer2, deviceBuffer1, bConjugate, m_uiX, m_uiY, m_uiX * m_uiY);

    OnChangeContent();
    UINT uiY = m_uiY;
    m_uiY = m_uiX;
    m_uiX = uiY;

    checkCudaErrors(cudaMemcpy(m_pData->m_pData, deviceBuffer2, sizeof(QLComplex) * m_uiX * m_uiY, cudaMemcpyDeviceToHost));

    checkCudaErrors(cudaFree(deviceBuffer1));
    checkCudaErrors(cudaFree(deviceBuffer2));
}


QLMatrix QLMatrix::SquareMatrixByAddOne() const
{
    if (m_uiX == m_uiY)
    {
        return *this;
    }

    if (m_uiX > m_uiY)
    {
        QLComplex* pBuffer = reinterpret_cast<QLComplex*>(calloc(m_uiX * m_uiX, sizeof(QLComplex)));
        for (UINT i = 0; i < m_uiX; ++i)
        {
            memcpy(pBuffer + i * m_uiX, HostBuffer() + i * m_uiY, sizeof(QLComplex) * m_uiY);
        }
        for (UINT i = m_uiY; i < m_uiX; ++i)
        {
            pBuffer[i * m_uiX + i] = _make_cuComplex(F(1.0), F(0.0));
        }
        return QLMatrix(m_uiX, m_uiX, pBuffer);
    }

    QLComplex* pBuffer = reinterpret_cast<QLComplex*>(calloc(m_uiY * m_uiY, sizeof(QLComplex)));
    memcpy(pBuffer, HostBuffer(), sizeof(QLComplex) * m_uiX * m_uiY);
    for (UINT i = m_uiX; i < m_uiY; ++i)
    {
        pBuffer[i * m_uiY + i] = _make_cuComplex(F(1.0), F(0.0));
    }
    return QLMatrix(m_uiY, m_uiY, pBuffer);
}

QLMatrix QLMatrix::ExtendZeros(UINT uiNewX, UINT uiNewY) const
{
    if (uiNewX == m_uiX && uiNewY == m_uiY)
    {
        return *this;
    }

    if (uiNewX < m_uiX || uiNewY < m_uiY)
    {
        appCrucial("become smaller!\n");
        return *this;
    }

    QLComplex* pBuffer = reinterpret_cast<QLComplex*>(calloc(uiNewX * uiNewY, sizeof(QLComplex)));
    for (UINT x = 0; x < m_uiX; ++x)
    {
        memcpy(pBuffer + x * uiNewY, HostBuffer() + x * m_uiY, m_uiY * sizeof(QLComplex));
    }
    return QLMatrix(uiNewX, uiNewY, pBuffer);
}

void QLMatrix::Opposite()
{
    dim3 blocks;
    dim3 threads;
    GetDim(blocks, threads);

    QLComplex* deviceBuffer = NULL;
    checkCudaErrors(cudaMalloc((void**)&deviceBuffer, sizeof(QLComplex) * m_uiX * m_uiY));
    checkCudaErrors(cudaMemcpy(deviceBuffer, HostBuffer(), sizeof(QLComplex) * m_uiX * m_uiY, cudaMemcpyHostToDevice));
    _kernelOpposite << <blocks, threads >> > (deviceBuffer, m_uiY, m_uiX * m_uiY);

    OnChangeContent();
    checkCudaErrors(cudaMemcpy(m_pData->m_pData, deviceBuffer, sizeof(QLComplex) * m_uiX * m_uiY, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaFree(deviceBuffer));
}

void QLMatrix::Add(const QLMatrix& other)
{
    if (m_uiX != other.m_uiX || m_uiY != other.m_uiY)
    {
        appCrucial("Add of two matrix not much!\n");
        return;
    }

    dim3 blocks;
    dim3 threads;
    GetDim(blocks, threads);
    QLComplex* deviceBuffer1 = NULL;
    QLComplex* deviceBuffer2 = NULL;
    checkCudaErrors(cudaMalloc((void**)&deviceBuffer1, sizeof(QLComplex) * m_uiX * m_uiY));
    checkCudaErrors(cudaMalloc((void**)&deviceBuffer2, sizeof(QLComplex) * m_uiX * m_uiY));
    checkCudaErrors(cudaMemcpy(deviceBuffer1, HostBuffer(), sizeof(QLComplex) * m_uiX * m_uiY, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(deviceBuffer2, other.HostBuffer(), sizeof(QLComplex) * m_uiX * m_uiY, cudaMemcpyHostToDevice));
    _kernelAdd << <blocks, threads >> > (deviceBuffer1, deviceBuffer2, m_uiY, m_uiX * m_uiY);

    OnChangeContent();

    checkCudaErrors(cudaMemcpy(m_pData->m_pData, deviceBuffer1, sizeof(QLComplex) * m_uiX * m_uiY, cudaMemcpyDeviceToHost));

    checkCudaErrors(cudaFree(deviceBuffer1));
    checkCudaErrors(cudaFree(deviceBuffer2));
}

void QLMatrix::Add(const QLComplex& other)
{
    dim3 blocks;
    dim3 threads;
    GetDimDiag(blocks, threads);

    QLComplex* deviceBuffer = NULL;
    checkCudaErrors(cudaMalloc((void**)&deviceBuffer, sizeof(QLComplex) * m_uiX * m_uiY));
    checkCudaErrors(cudaMemcpy(deviceBuffer, HostBuffer(), sizeof(QLComplex) * m_uiX * m_uiY, cudaMemcpyHostToDevice));
    _kernelAdd << <blocks, threads >> > (deviceBuffer, other, m_uiY, m_uiX * m_uiY);

    OnChangeContent();
    checkCudaErrors(cudaMemcpy(m_pData->m_pData, deviceBuffer, sizeof(QLComplex) * m_uiX * m_uiY, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaFree(deviceBuffer));
}

void QLMatrix::Add(const Real& other)
{
    dim3 blocks;
    dim3 threads;
    GetDimDiag(blocks, threads);

    QLComplex* deviceBuffer = NULL;
    checkCudaErrors(cudaMalloc((void**)&deviceBuffer, sizeof(QLComplex) * m_uiX * m_uiY));
    checkCudaErrors(cudaMemcpy(deviceBuffer, HostBuffer(), sizeof(QLComplex) * m_uiX * m_uiY, cudaMemcpyHostToDevice));
    _kernelAdd << <blocks, threads >> > (deviceBuffer, other, m_uiY, m_uiX * m_uiY);

    OnChangeContent();
    checkCudaErrors(cudaMemcpy(m_pData->m_pData, deviceBuffer, sizeof(QLComplex) * m_uiX * m_uiY, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaFree(deviceBuffer));
}

//======================================================================================================
// C = alpha A.B + beta C
// If beta==0, C does not have to be a valid input
//from https://github.com/NVIDIA/CUDALibrarySamples/blob/master/cuBLAS/Level-3/gemm/cublas_gemm_example.cu
QLMatrix QLMatrix::Mul(const QLMatrix& other, cublasOperation_t left, cublasOperation_t right) const
{
    if (left == CUBLAS_OP_N && right == CUBLAS_OP_N && m_uiX != other.m_uiY)
    {
        appCrucial("me (%d x %d) cannot mult other (%d x %d)!\n", m_uiX, m_uiY, other.m_uiX, other.m_uiY);
        return QLMatrix();
    }
    if (left != CUBLAS_OP_N && right == CUBLAS_OP_N && m_uiY != other.m_uiY)
    {
        appCrucial("me^+ (%d x %d) cannot mult other (%d x %d)!\n", m_uiX, m_uiY, other.m_uiX, other.m_uiY);
        return QLMatrix();
    }
    if (left == CUBLAS_OP_N && right != CUBLAS_OP_N && m_uiX != other.m_uiX)
    {
        appCrucial("me (%d x %d) cannot mult other^+ (%d x %d)!\n", m_uiX, m_uiY, other.m_uiX, other.m_uiY);
        return QLMatrix();
    }
    if (left != CUBLAS_OP_N && right != CUBLAS_OP_N && m_uiY != other.m_uiX)
    {
        appCrucial("me^+ (%d x %d) cannot mult other^+ (%d x %d)!\n", m_uiX, m_uiY, other.m_uiX, other.m_uiY);
        return QLMatrix();
    }

    cublasHandle_t cublasH = NULL;
    cudaStream_t stream = NULL;

    //This is after dagger
    //me = m*k
    //other = k*n
    //res = m*n
    //y first
    INT m = m_uiY; //number of rows of matrix op(A) and C.
    INT k = m_uiX; //number of columns of op(A) and rows of op(B).
    INT n = other.m_uiX; //number of columns of matrix op(B) and C.
    if (left != CUBLAS_OP_N && right == CUBLAS_OP_N)
    {
        m = m_uiX;
        k = m_uiY;
    }
    if (left != CUBLAS_OP_N && right != CUBLAS_OP_N)
    {
        m = m_uiX;
        k = m_uiY;
        n = other.m_uiY;
    }
    if (left == CUBLAS_OP_N && right != CUBLAS_OP_N)
    {
        n = other.m_uiY;
    }

    //This is before dagger
    const INT lda = m_uiY;
    const INT ldb = other.m_uiY;
    //This is uiY of output matrix
    INT ldc = m_uiY;
    if (left != CUBLAS_OP_N)
    {
        ldc = m_uiX;
    }
    
    /*
     *   A = | 1.0 | 2.0 |
     *       | 3.0 | 4.0 |
     *
     *   B = | 5.0 | 6.0 |
     *       | 7.0 | 8.0 |
     */

    //const std::vector<data_type> A = { 1.0, 2.0, 3.0, 4.0 };
    //const std::vector<data_type> B = { 5.0, 6.0, 7.0, 8.0 };
    //std::vector<data_type> C(m * n);
    const QLComplex alpha = _make_cuComplex(F(1.0), F(0.0));
    const QLComplex beta = _make_cuComplex(F(0.0), F(0.0));

    QLComplex* d_A = NULL;
    QLComplex* d_B = NULL;
    QLComplex* d_C = NULL;

    cublasOperation_t transa = left;
    cublasOperation_t transb = right;


    /* step 1: create cublas handle, bind a stream */
    checkCudaErrors(cublasCreate(&cublasH));

    checkCudaErrors(cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking));
    checkCudaErrors(cublasSetStream(cublasH, stream));

    /* step 2: copy data to device */
    checkCudaErrors(cudaMalloc(reinterpret_cast<void**>(&d_A), sizeof(QLComplex) * m_uiX * m_uiY));
    checkCudaErrors(cudaMalloc(reinterpret_cast<void**>(&d_B), sizeof(QLComplex) * other.m_uiX * other.m_uiY));
    checkCudaErrors(cudaMalloc(reinterpret_cast<void**>(&d_C), sizeof(QLComplex) * m * n));

    checkCudaErrors(cudaMemcpyAsync(d_A, HostBuffer(), sizeof(QLComplex) * m_uiX * m_uiY, cudaMemcpyHostToDevice, stream));
    checkCudaErrors(cudaMemcpyAsync(d_B, other.HostBuffer(), sizeof(QLComplex) * other.m_uiX * other.m_uiY, cudaMemcpyHostToDevice, stream));

    /* step 3: compute */
#if _QL_DOUBLEFLOAT
    checkCudaErrors(cublasZgemm(cublasH, transa, transb, m, n, k, &alpha, d_A, lda, d_B, ldb, &beta, d_C, ldc));
#else
    checkCudaErrors(cublasCgemm(cublasH, transa, transb, m, n, k, &alpha, d_A, lda, d_B, ldb, &beta, d_C, ldc));
#endif

    /* step 4: copy data to host */
    QLComplex* cbuffer = reinterpret_cast<QLComplex*>(malloc(m * n * sizeof(QLComplex)));
    checkCudaErrors(cudaMemcpyAsync(cbuffer, d_C, m* n * sizeof(QLComplex), cudaMemcpyDeviceToHost, stream));
    QLMatrix ret(n, m, cbuffer);

    checkCudaErrors(cudaStreamSynchronize(stream));

    /*
     *   C = | 23.0 | 31.0 |
     *       | 34.0 | 46.0 |
     */

    //printf("C\n");
    //print_matrix(m, n, C.data(), ldc);
    //printf("=====\n");

    /* free resources */
    checkCudaErrors(cudaFree(d_A));
    checkCudaErrors(cudaFree(d_B));
    checkCudaErrors(cudaFree(d_C));

    checkCudaErrors(cublasDestroy(cublasH));

    return ret;
}

void QLMatrix::Mul(const QLComplex& other)
{
    dim3 blocks;
    dim3 threads;
    GetDim(blocks, threads);

    QLComplex* deviceBuffer = NULL;
    checkCudaErrors(cudaMalloc((void**)&deviceBuffer, sizeof(QLComplex) * m_uiX * m_uiY));
    checkCudaErrors(cudaMemcpy(deviceBuffer, HostBuffer(), sizeof(QLComplex) * m_uiX * m_uiY, cudaMemcpyHostToDevice));
    _kernelMul << <blocks, threads >> > (deviceBuffer, other, m_uiY, m_uiX * m_uiY);

    OnChangeContent();
    checkCudaErrors(cudaMemcpy(m_pData->m_pData, deviceBuffer, sizeof(QLComplex) * m_uiX * m_uiY, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaFree(deviceBuffer));
}

void QLMatrix::Mul(const Real& other)
{
    dim3 blocks;
    dim3 threads;
    GetDim(blocks, threads);

    QLComplex* deviceBuffer = NULL;
    checkCudaErrors(cudaMalloc((void**)&deviceBuffer, sizeof(QLComplex) * m_uiX * m_uiY));
    checkCudaErrors(cudaMemcpy(deviceBuffer, HostBuffer(), sizeof(QLComplex) * m_uiX * m_uiY, cudaMemcpyHostToDevice));
    _kernelMul << <blocks, threads >> > (deviceBuffer, other, m_uiY, m_uiX * m_uiY);

    OnChangeContent();
    checkCudaErrors(cudaMemcpy(m_pData->m_pData, deviceBuffer, sizeof(QLComplex) * m_uiX * m_uiY, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaFree(deviceBuffer));
}

void QLMatrix::Div(const QLComplex& other)
{
    dim3 blocks;
    dim3 threads;
    GetDim(blocks, threads);

    QLComplex* deviceBuffer = NULL;
    checkCudaErrors(cudaMalloc((void**)&deviceBuffer, sizeof(QLComplex) * m_uiX * m_uiY));
    checkCudaErrors(cudaMemcpy(deviceBuffer, HostBuffer(), sizeof(QLComplex) * m_uiX * m_uiY, cudaMemcpyHostToDevice));
    _kernelDiv << <blocks, threads >> > (deviceBuffer, other, m_uiY, m_uiX * m_uiY);

    OnChangeContent();
    checkCudaErrors(cudaMemcpy(m_pData->m_pData, deviceBuffer, sizeof(QLComplex) * m_uiX * m_uiY, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaFree(deviceBuffer));
}

void QLMatrix::Div(const Real& other)
{
    dim3 blocks;
    dim3 threads;
    GetDim(blocks, threads);

    QLComplex* deviceBuffer = NULL;
    checkCudaErrors(cudaMalloc((void**)&deviceBuffer, sizeof(QLComplex) * m_uiX * m_uiY));
    checkCudaErrors(cudaMemcpy(deviceBuffer, HostBuffer(), sizeof(QLComplex) * m_uiX * m_uiY, cudaMemcpyHostToDevice));
    _kernelDiv << <blocks, threads >> > (deviceBuffer, other, m_uiY, m_uiX * m_uiY);

    OnChangeContent();
    checkCudaErrors(cudaMemcpy(m_pData->m_pData, deviceBuffer, sizeof(QLComplex) * m_uiX * m_uiY, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaFree(deviceBuffer));
}

QLMatrix QLMatrix::GetBlock(UINT uiXStart, UINT uiXLen, UINT uiYStart, UINT uiYLen) const
{
    if (uiXStart + uiXLen > m_uiX)
    {
        appCrucial("Block larger than the matrix\n");
    }

    if (uiYStart + uiYLen > m_uiY)
    {
        appCrucial("Block larger than the matrix\n");
    }

    QLComplex* pBuffer = reinterpret_cast<QLComplex*>(malloc(uiXLen * uiYLen * sizeof(QLComplex)));
    for (UINT uiX = 0; uiX < uiXLen; ++uiX)
    {
        memcpy(pBuffer + uiYLen * uiX, HostBuffer() + m_uiY * (uiXStart + uiX) + uiYStart, uiYLen * sizeof(QLComplex));
    }

    return QLMatrix(uiXLen, uiYLen, pBuffer);
}

void QLMatrix::SetBlock(UINT uiXStart, UINT uiXLen, UINT uiYStart, UINT uiYLen, const QLComplex* content)
{
    if (uiXStart + uiXLen > m_uiX)
    {
        appCrucial("Block larger than the matrix\n");
        return;
    }

    if (uiYStart + uiYLen > m_uiY)
    {
        appCrucial("Block larger than the matrix\n");
        return;
    }

    OnChangeContent();
    for (UINT uiX = 0; uiX < uiXLen; ++uiX)
    {
        memcpy(m_pData->m_pData + m_uiY * (uiXStart + uiX) + uiYStart, content + uiYLen * uiX, uiYLen * sizeof(QLComplex));
    }
}

QLMatrix QLMatrix::GetDiag() const
{
    QLComplex* pBuffer = reinterpret_cast<QLComplex*>(calloc(m_uiX * m_uiY, sizeof(QLComplex)));
    UINT uiMinXY = std::min(m_uiX, m_uiY);
    for (UINT i = 0; i < uiMinXY; ++i)
    {
        pBuffer[i * m_uiY + i] = HostBuffer()[i * m_uiY + i];
    }
    return QLMatrix(m_uiX, m_uiY, pBuffer);
}

QLMatrix QLMatrix::BlockAdd(const QLMatrix& m, BYTE addtype) const
{
    UINT uiLengthAllX = (m_uiX + m.m_uiX);
    UINT uiLengthAllY = (m_uiY + m.m_uiY);
    QLComplex* pBuffer = NULL;
    UINT maxY = std::max(m_uiY, m.m_uiY);

    switch (addtype)
    {
    case 0:
        pBuffer = reinterpret_cast<QLComplex*>(calloc(uiLengthAllX * uiLengthAllY, sizeof(QLComplex)));
        for (UINT i = 0; i < m_uiX; ++i)
        {
            memcpy(pBuffer + i * uiLengthAllY, HostBuffer() + i * m_uiY, m_uiY * sizeof(QLComplex));
        }
        for (UINT i = 0; i < m.m_uiX; ++i)
        {
            memcpy(pBuffer + m_uiY + (i + m_uiX) * uiLengthAllY, m.HostBuffer() + i * m.m_uiY, m.m_uiY * sizeof(QLComplex));
        }
        return QLMatrix(uiLengthAllX, uiLengthAllY, pBuffer);
    case 1:
        pBuffer = reinterpret_cast<QLComplex*>(calloc(uiLengthAllX * maxY, sizeof(QLComplex)));
        for (UINT i = 0; i < m_uiX; ++i)
        {
            memcpy(pBuffer + i * maxY, HostBuffer() + i * m_uiY, m_uiY * sizeof(QLComplex));
        }
        for (UINT i = 0; i < m.m_uiX; ++i)
        {
            memcpy(pBuffer + (i + m_uiX) * maxY, m.HostBuffer() + i * m.m_uiY, m.m_uiY * sizeof(QLComplex));
        }
        return QLMatrix(uiLengthAllX, maxY, pBuffer);
    default:
        break;
    }

    UINT maxX = std::max(m_uiX, m.m_uiX);
    pBuffer = reinterpret_cast<QLComplex*>(calloc(maxX * uiLengthAllY, sizeof(QLComplex)));
    for (UINT i = 0; i < m_uiX; ++i)
    {
        memcpy(pBuffer + i * uiLengthAllY, HostBuffer() + i * m_uiY, m_uiY * sizeof(QLComplex));
    }
    for (UINT i = 0; i < m.m_uiX; ++i)
    {
        memcpy(pBuffer + i * uiLengthAllY + m_uiY, m.HostBuffer() + i * m.m_uiY, m.m_uiY * sizeof(QLComplex));
    }
    return QLMatrix(maxX, uiLengthAllY, pBuffer);
}

void QLMatrix::CSD2BY1(QLMatrix& u1, QLMatrix& u2, QLMatrix& c, QLMatrix& s, QLMatrix& v, UINT uiSep) const
{
    if (uiSep < 1 || uiSep > (m_uiY - 1))
    {
        appWarning("seperation of the dimension failed\n");
        return;
    }

    QLMatrix q1 = GetBlock(0, m_uiX, 0, uiSep);
    QLMatrix q2 = GetBlock(0, m_uiX, uiSep, (m_uiY - uiSep));

    q2.SVDJ(u2, s, v);
    QLMatrix x = q1 * v;
    QLMatrix r;
    x.QR(u1, r);
    c = r.GetDiag();
}

void QLMatrix::CSD(QLMatrix& u1, QLMatrix& u2, QLMatrix& c, QLMatrix& s, QLMatrix& v1, QLMatrix& v2, UINT uiXSep, UINT uiYSep) const
{
    QLMatrix leftBlock = GetBlock(0, uiXSep, 0, m_uiY);
    leftBlock.CSD2BY1(u1, u2, c, s, v1, uiYSep);

    assert(u1.m_uiX == uiYSep);
    assert(u1.m_uiY == uiYSep);
    assert(u2.m_uiX == m_uiY - uiYSep);
    assert(u2.m_uiY == m_uiY - uiYSep);

    assert(c.m_uiX == uiXSep);
    assert(s.m_uiX == uiXSep);

    assert(v1.m_uiX == uiXSep);
    assert(v1.m_uiY == uiXSep);

    //c is uiXSep * uiYSep, with uiXSep <= uiYSep
    //s is uiXSep * (m_uiY - uiYSep), with uiXSep <= (m_uiY - uiYSep)
    //so both c and s are uiXSep * uiXSep
    if (c.m_uiY > c.m_uiX)
    {
        c = c.GetBlock(0, c.m_uiX, 0, c.m_uiX);
    }
    if (s.m_uiY > s.m_uiX)
    {
        s = s.GetBlock(0, s.m_uiX, 0, s.m_uiX);
    }

    QLMatrix q12 = GetBlock(uiXSep, m_uiX - uiXSep, 0, uiYSep);
    QLMatrix q22 = GetBlock(uiXSep, m_uiX - uiXSep, uiYSep, m_uiY - uiYSep);

    //since q12^+ can dot u11, we need u11.x = q12.x, which is uiXSep
    assert(uiXSep <= uiYSep);
    QLMatrix u11 = u1.GetBlock(0, uiXSep, 0, uiYSep);
    //since q22^+ can dot u21, we need u21.x = q22.x, which is uiXSep
    assert(uiXSep <= m_uiY - uiYSep);
    QLMatrix u21 = u2.GetBlock(0, uiXSep, 0, m_uiY - uiYSep);

    //now, q12^+.u11 = q12.x * u11.x which is xy = (uiXSep, m_uiX - uiXSep) 
    //now, q12^+.u22 = q22.x * u21.x which is xy = (uiXSep, m_uiX - uiXSep) 

    //we need c,s to be (m_uiX - uiXSep)*uiXSep, it is now uiXSep*uiXSep
    QLMatrix mulc = c;
    QLMatrix muls = s;
    if (m_uiX - uiXSep < uiXSep)
    {
        mulc = c.GetBlock(0, m_uiX - uiXSep, 0, uiXSep);
        muls = s.GetBlock(0, m_uiX - uiXSep, 0, uiXSep);
    }
    if (m_uiX - uiXSep > uiXSep)
    {
        mulc = c.ExtendZeros(m_uiX - uiXSep, uiXSep);
        muls = s.ExtendZeros(m_uiX - uiXSep, uiXSep);
    }

    v2 = q22.Mul(u21, CUBLAS_OP_C) * mulc - q12.Mul(u11, CUBLAS_OP_C) * muls;
}

QLComplex QLMatrix::VectorDot(const QLMatrix& other, UBOOL bConjL, UBOOL bConjR) const
{
    UINT uiL = m_uiX * m_uiY;
    if (uiL != other.m_uiX * other.m_uiY)
    {
        appWarning(_T("matrix size not same!\n"));
        const UINT uiL2 = other.m_uiX * other.m_uiY;
        if (uiL2 < uiL)
        {
            uiL = uiL2;
        }
    }
    dim3 blocks;
    dim3 threads;
    GetDim(blocks, threads);
    QLComplex* deviceBuffer1 = NULL;
    QLComplex* deviceBuffer2 = NULL;
    checkCudaErrors(cudaMalloc((void**)&deviceBuffer1, sizeof(QLComplex) * uiL));
    checkCudaErrors(cudaMalloc((void**)&deviceBuffer2, sizeof(QLComplex) * uiL));
    checkCudaErrors(cudaMemcpy(deviceBuffer1, HostBuffer(), sizeof(QLComplex) * uiL, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(deviceBuffer2, other.HostBuffer(), sizeof(QLComplex) * uiL, cudaMemcpyHostToDevice));

    _kernelDot<<<blocks, threads>>>(deviceBuffer1, deviceBuffer2, m_uiY, uiL, bConjL, bConjR);

    QLComplex res = ReduceSum(deviceBuffer1, uiL);

    checkCudaErrors(cudaFree(deviceBuffer1));
    checkCudaErrors(cudaFree(deviceBuffer2));

    return res;
}

void QLMatrix::CombinedCSD(QLMatrix& u, QLMatrix& cs, QLMatrix& v, const QLMatrix& u1, const QLMatrix& u2, const QLMatrix& c, const QLMatrix& s, const QLMatrix& v1, const QLMatrix& v2, UINT uiXSep, UINT uiYSep)
{
    u = u1.BlockAdd(u2, 0);
    v = v1.BlockAdd(v2, 0);
    
    //const UINT outX = v.m_uiY;
    //const UINT outY = u.m_uiY;
    const UINT csX = v.m_uiX;
    const UINT csY = u.m_uiX;

    assert(c.m_uiX == c.m_uiY);
    assert(s.m_uiX == s.m_uiY);
    assert(c.m_uiX == s.m_uiX);

    assert(c.m_uiX == uiXSep);
    assert(uiXSep * 2 <= csX);
    assert(c.m_uiY <= uiYSep);
    assert(c.m_uiY + uiYSep <= csY);

    // cs = | C | -S 0 0 |
    //      | 0 |  0 I 0 |
    //      --------------
    //      | S |  C 0 0 |
    //      | 0 |  0 0 I |

    QLMatrix cup = c.ExtendZeros(c.m_uiX, uiYSep);
    QLMatrix cdown = c.ExtendZeros(c.m_uiX, csY - uiYSep);
    QLMatrix sup = s.ExtendZeros(s.m_uiX, uiYSep);
    QLMatrix sdown = s.ExtendZeros(s.m_uiX, csY - uiYSep);
    sup.Opposite();

    cs = cup.BlockAdd(sdown, 2);
    QLMatrix csright = sup.BlockAdd(cdown, 2);
    cs = cs.BlockAdd(csright, 1);

    if (c.m_uiY < uiYSep)
    {
        csright = QLMatrix(uiYSep - c.m_uiY, c.m_uiY);
        csright = csright.BlockAdd(CreateEye(uiYSep - c.m_uiY, uiYSep - c.m_uiY), 2);
        csright = csright.BlockAdd(QLMatrix(uiYSep - c.m_uiY, csY - uiYSep), 2);
        cs = cs.BlockAdd(csright, 1);
    }

    if (c.m_uiY < csY - uiYSep)
    {
        const UINT toAdd = csY - uiYSep - c.m_uiY;
        csright = QLMatrix(csY - uiYSep - c.m_uiY, uiYSep + c.m_uiY);
        csright = csright.BlockAdd(CreateEye(toAdd, toAdd), 2);
        cs = cs.BlockAdd(csright, 1);
    }
}

/**
* see https://github.com/NVIDIA/CUDALibrarySamples/blob/master/cuSOLVER/Xsyevd/cusolver_Xsyevd_example.cu
*/
void QLMatrix::EVD(QLMatrix& vm, QLMatrix& lm) const
{
    if (m_uiX != m_uiY)
    {
        appCrucial(_T("EVD only supports squared matrix!\n"));
        return;
    }

    cusolverDnHandle_t cusolverH = NULL;
    cudaStream_t stream = NULL;

    const INT m = static_cast<INT>(m_uiX);
    const INT lda = static_cast<INT>(m_uiY);
    /*
     *       | 3.5 0.5 0.0 |
     *   A = | 0.5 3.5 0.0 |
     *       | 0.0 0.0 2.0 |
     *
     */
    //const std::vector<data_type> A = { 3.5, 0.5, 0.0, 0.5, 3.5, 0.0, 0.0, 0.0, 2.0 };
    //const std::vector<data_type> lambda = { 2.0, 3.0, 4.0 };

    //std::vector<data_type> V(lda * m, 0); // eigenvectors
    //std::vector<data_type> W(m, 0);       // eigenvalues

    QLComplex* d_A = NULL;
    Real* d_W = NULL;
    INT* d_info = NULL;

    INT info = 0;

    INT d_lwork = 0;     /* size of workspace */
    QLComplex* d_work = NULL; /* device workspace */

    //std::printf("A = (matlab base-1)\n");
    //print_matrix(m, m, A.data(), lda);
    //std::printf("=====\n");

    /* step 1: create cusolver handle, bind a stream */
    checkCudaErrors(cusolverDnCreate(&cusolverH));

    checkCudaErrors(cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking));
    checkCudaErrors(cusolverDnSetStream(cusolverH, stream));

    checkCudaErrors(cudaMalloc(reinterpret_cast<void**>(&d_A), sizeof(QLComplex) * m * lda));
    checkCudaErrors(cudaMalloc(reinterpret_cast<void**>(&d_W), sizeof(Real) * m));
    checkCudaErrors(cudaMalloc(reinterpret_cast<void**>(&d_info), sizeof(INT)));

    checkCudaErrors(cudaMemcpyAsync(d_A, HostBuffer(), sizeof(QLComplex) * m * lda, cudaMemcpyHostToDevice, stream));

    // step 3: query working space of syevd
    cusolverEigMode_t jobz = CUSOLVER_EIG_MODE_VECTOR; // compute eigenvalues and eigenvectors.
    cublasFillMode_t uplo = CUBLAS_FILL_MODE_LOWER;

#if _QL_DOUBLEFLOAT
    checkCudaErrors(cusolverDnZheevd_bufferSize(cusolverH, jobz, uplo, m, d_A, lda, d_W, &d_lwork));
#else
    checkCudaErrors(cusolverDnCheevd_bufferSize(cusolverH, jobz, uplo, m, d_A, lda, d_W, &d_lwork));
#endif

    checkCudaErrors(cudaMalloc(reinterpret_cast<void**>(&d_work), sizeof(QLComplex) * d_lwork));

    // step 4: compute spectrum
#if _QL_DOUBLEFLOAT
    checkCudaErrors(cusolverDnZheevd(cusolverH, jobz, uplo, m, d_A, lda, d_W, d_work, d_lwork, d_info));
#else
    checkCudaErrors(cusolverDnCheevd(cusolverH, jobz, uplo, m, d_A, lda, d_W, d_work, d_lwork, d_info));
#endif


    QLComplex* v = reinterpret_cast<QLComplex*>(malloc(sizeof(QLComplex) * lda * m));
    Real* w = reinterpret_cast<Real*>(malloc(sizeof(Real) * m));

    checkCudaErrors(cudaMemcpyAsync(v, d_A, sizeof(QLComplex) * lda * m, cudaMemcpyDeviceToHost, stream));
    checkCudaErrors(cudaMemcpyAsync(w, d_W, sizeof(Real) * m, cudaMemcpyDeviceToHost, stream));
    checkCudaErrors(cudaMemcpyAsync(&info, d_info, sizeof(INT), cudaMemcpyDeviceToHost, stream));

    checkCudaErrors(cudaStreamSynchronize(stream));

    if (0 > info) 
    {
        appCrucial("%d-th parameter is wrong \n", -info);
        checkCudaErrors(cudaFree(d_A));
        checkCudaErrors(cudaFree(d_W));
        checkCudaErrors(cudaFree(d_info));
        checkCudaErrors(cudaFree(d_work));

        appSafeFree(v);
        appSafeFree(w);

        return;
    }

    /* free resources */
    QLComplex* cw = reinterpret_cast<QLComplex*>(malloc(sizeof(QLComplex) * m));
    for (INT i = 0; i < m; ++i)
    {
        cw[i] = _make_cuComplex(w[i], F(0.0));
    }

    checkCudaErrors(cudaFree(d_A));
    checkCudaErrors(cudaFree(d_W));
    checkCudaErrors(cudaFree(d_info));
    checkCudaErrors(cudaFree(d_work));
    appSafeFree(w);

    checkCudaErrors(cusolverDnDestroy(cusolverH));
    checkCudaErrors(cudaStreamDestroy(stream));

    vm = QLMatrix(m_uiY, m_uiX, v);
    lm = QLMatrix(m_uiX, 1, cw);
}

void QLMatrix::MatrixFunction(complexfunc func)
{
    if (m_uiX != m_uiY)
    {
        appCrucial(_T("EVD only supports squared matrix!\n"));
        return;
    }

    QLMatrix v, w;
    EVD(v, w);

    w.ElementWiseFunction(func);

    OnChangeContent();
    memset(m_pData->m_pData, 0, sizeof(QLComplex) * m_uiX * m_uiY);
    for (UINT i = 0; i < m_uiX; ++i)
    {
        m_pData->m_pData[i * m_uiY + i] = w.Get(i, 0);
    }
    Print("w");

    *this = v * (*this);
    v.Dagger();
    *this = (*this) * v;
}

void QLMatrix::MatrixFunctionTwoR(complexfuncTwoR func, Real r)
{
    if (m_uiX != m_uiY)
    {
        appCrucial(_T("EVD only supports squared matrix!\n"));
        return;
    }

    QLMatrix v, w;
    EVD(v, w);
    w.ElementWiseFunctionTwoR(func, r);

    OnChangeContent();
    memset(m_pData->m_pData, 0, sizeof(QLComplex) * m_uiX * m_uiY);
    for (UINT i = 0; i < m_uiX; ++i)
    {
        m_pData->m_pData[i * m_uiY + i] = w.Get(i, 0);
    }

    *this = v * (*this);
    v.Dagger();
    *this = (*this) * v;
}

void QLMatrix::ElementWiseFunction(complexfunc func)
{
    QLComplex* deviceBuffer;
    checkCudaErrors(cudaMalloc(reinterpret_cast <void**>(&deviceBuffer), sizeof(QLComplex) * m_uiX * m_uiY));
    checkCudaErrors(cudaMemcpy(deviceBuffer, HostBuffer(), sizeof(QLComplex) * m_uiX * m_uiY, cudaMemcpyHostToDevice));

    dim3 blocks;
    dim3 threads;
    GetDim(blocks, threads);
    _kernelComplexFunc << <blocks, threads >> > (deviceBuffer, m_uiY, m_uiX * m_uiY, func);
    checkCudaErrors(cudaGetLastError());
    OnChangeContent();
    checkCudaErrors(cudaMemcpy(m_pData->m_pData, deviceBuffer, sizeof(QLComplex) * m_uiX * m_uiY, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaFree(deviceBuffer));
}

void QLMatrix::ElementWiseFunctionTwoR(complexfuncTwoR func, Real r)
{
    QLComplex* deviceBuffer;
    checkCudaErrors(cudaMalloc(reinterpret_cast <void**>(&deviceBuffer), sizeof(QLComplex) * m_uiX * m_uiY));
    checkCudaErrors(cudaMemcpy(deviceBuffer, HostBuffer(), sizeof(QLComplex) * m_uiX * m_uiY, cudaMemcpyHostToDevice));

    dim3 blocks;
    dim3 threads;
    GetDim(blocks, threads);
    _kernelComplexFuncTwoR << <blocks, threads >> > (deviceBuffer, m_uiY, m_uiX * m_uiY, func, r);
    checkCudaErrors(cudaGetLastError());
    OnChangeContent();
    checkCudaErrors(cudaMemcpy(m_pData->m_pData, deviceBuffer, sizeof(QLComplex) * m_uiX * m_uiY, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaFree(deviceBuffer));
}

void QLMatrix::MatrixIExp(Real t)
{
    complexfuncTwoR host_func;
    checkCudaErrors(cudaMemcpyFromSymbol(&host_func, devicefunctionIExp, sizeof(complexfuncTwoR)));

    MatrixFunctionTwoR(host_func, t);
}

void QLMatrix::MatrixExp()
{
    complexfunc host_func;
    checkCudaErrors(cudaMemcpyFromSymbol(&host_func, devicefunctionExp, sizeof(complexfunc)));

    MatrixFunction(host_func);
}

void QLMatrix::MatrixSqrt()
{
    complexfunc host_func;
    checkCudaErrors(cudaMemcpyFromSymbol(&host_func, devicefunctionSqrt, sizeof(complexfunc)));

    MatrixFunction(host_func);
}

QLMatrix QLMatrix::CreateEye(UINT uiX, UINT uiY)
{
    QLComplex* pBuffer = (QLComplex*)calloc(uiX * uiY, sizeof(QLComplex));
    UINT minXY = std::min(uiX, uiY);
    for (UINT i = 0; i < minXY; ++i)
    {
        pBuffer[i * uiY + i] = _make_cuComplex(F(1.0), F(0.0));
    }
    return QLMatrix(uiX, uiY, pBuffer);
}

QLComplex QLMatrix::ReduceSum(QLComplex* deviceBuffer, UINT uiLength)
{
    const UINT iRequiredDim = (uiLength + 1) >> 1;
    const UINT iPower = MostSignificantPowerTwo(iRequiredDim);
    for (UINT i = 0; i <= iPower; ++i)
    {
        UINT iJump = 1 << i;
        UINT iThreadNeeded = 1 << (iPower - i);
        UINT iBlock = iThreadNeeded > _QL_LAUNCH_MAX_THREAD ? iThreadNeeded / _QL_LAUNCH_MAX_THREAD : 1;
        UINT iThread = iThreadNeeded > _QL_LAUNCH_MAX_THREAD ? _QL_LAUNCH_MAX_THREAD : iThreadNeeded;
        _kernelReduceComplex << <iBlock, iThread >> > (deviceBuffer, iJump, uiLength);
    }
    QLComplex result[1];
    checkCudaErrors(cudaMemcpy(result, deviceBuffer, sizeof(QLComplex), cudaMemcpyDeviceToHost));
    return result[0];
}

QLMatrix QLMatrix::VectorFFT(UBOOL bForward) const
{
    cufftHandle plan;
    INT n[1] = { static_cast<INT>(m_uiX * m_uiY) };

    const cufftResult planRes = cufftPlanMany(&plan, 1, n,
        n, 1, 1,
        n, 1, 1,
        CUFFT_Z2Z, 1);

    if (CUFFT_SUCCESS != planRes)
    {
        appCrucial(_T("cufftPlanMany failed! %d\n"), planRes);
        return QLMatrix();
    }

    QLComplex* buffer = NULL;
    checkCudaErrors(cudaMalloc((void**)&buffer, sizeof(QLComplex) * m_uiX * m_uiY));
    checkCudaErrors(cudaMemcpy(buffer, HostBuffer(), sizeof(QLComplex) * m_uiX * m_uiY, cudaMemcpyHostToDevice));

#if _QL_DOUBLEFLOAT
    const cufftResult res = cufftExecZ2Z(plan, buffer, buffer, bForward ? CUFFT_FORWARD : CUFFT_INVERSE);
#else
    const cufftResult res = cufftExecC2C(plan, buffer, buffer, bForward ? CUFFT_FORWARD : CUFFT_INVERSE);
#endif

    if (CUFFT_SUCCESS != res)
    {
        appCrucial(_T("cufftResult failed! %d\n"), res);
        checkCudaErrors(cudaFree(buffer));
        return QLMatrix();
    }

    QLComplex* hostBuffer = reinterpret_cast<QLComplex*>(malloc(sizeof(QLComplex) * m_uiX * m_uiY));
    checkCudaErrors(cudaMemcpy(hostBuffer, buffer, sizeof(QLComplex) * m_uiX * m_uiY, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaFree(buffer));
    return QLMatrix(m_uiX, m_uiY, hostBuffer);
}

TArray<QLComplex> QLMatrix::ToVector() const
{
    const QLComplex* buffer = HostBuffer();
    TArray<QLComplex> ret;
    for (UINT i = 0; i < m_uiX * m_uiY; ++i)
    {
        ret.AddItem(buffer[i]);
    }
    return ret;
}

QLMatrix QLMatrix::KroneckerProduct(const QLMatrix& m2) const
{
    const UINT resX = m_uiX * m2.m_uiX;
    const UINT resY = m_uiY * m2.m_uiY;

    QLComplex* resbuffer;
    QLComplex* m1buffer;
    QLComplex* m2buffer;
    checkCudaErrors(cudaMalloc((void**)&resbuffer, sizeof(QLComplex) * resX * resY));
    checkCudaErrors(cudaMalloc((void**)&m1buffer, sizeof(QLComplex) * m_uiX * m_uiY));
    checkCudaErrors(cudaMalloc((void**)&m2buffer, sizeof(QLComplex) * m2.m_uiX * m2.m_uiY));

    checkCudaErrors(cudaMemcpy(m1buffer, HostBuffer(), sizeof(QLComplex) * m_uiX * m_uiY, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(m2buffer, m2.HostBuffer(), sizeof(QLComplex) * m2.m_uiX * m2.m_uiY, cudaMemcpyHostToDevice));

    dim3 blocks;
    dim3 threads;
    GetDim(blocks, threads);

    const UINT uiMaxThreadAllowed = Ceil(_QL_LAUNCH_MAX_THREAD, threads.x * threads.y);
    const UINT uiThreadNeeded = m2.m_uiX * m2.m_uiY;
    UINT uiThreadZ = (uiThreadNeeded >= uiMaxThreadAllowed) ? uiMaxThreadAllowed : uiThreadNeeded;
    UINT uiBlockZ = (uiThreadZ > 1) ? Ceil(uiThreadNeeded, uiThreadZ) : uiThreadNeeded;
    if (uiBlockZ > 1)
    {
        uiThreadZ = Ceil(uiThreadNeeded, uiBlockZ);
    }

    blocks.z = uiBlockZ;
    threads.z = uiThreadZ;

    _kernelKronecker << <blocks, threads >> > (resbuffer, m1buffer, m2buffer, m_uiY, m2.m_uiY, m2.m_uiX, m_uiX * m_uiY, uiThreadNeeded);

    QLComplex* hostBuffer = reinterpret_cast<QLComplex*>(malloc(sizeof(QLComplex) * resX * resY));
    checkCudaErrors(cudaMemcpy(hostBuffer, resbuffer, sizeof(QLComplex) * resX * resY, cudaMemcpyDeviceToHost));

    checkCudaErrors(cudaFree(resbuffer));
    checkCudaErrors(cudaFree(m1buffer));
    checkCudaErrors(cudaFree(m2buffer));

    return QLMatrix(resX, resY, hostBuffer);
}

/**
* 
*/
QLMatrix QLMatrix::GELS(const QLMatrix& y) const
{
    cusolverDnHandle_t cusolverH = NULL;
    cudaStream_t stream = NULL;

    const INT m = m_uiX;

    INT info = 0;

    QLComplex* d_A = NULL;
    QLComplex* d_B = NULL;
    QLComplex* d_X = NULL; 
    INT* d_info = NULL; 

    SIZE_T d_lwork = 0; 
    QLComplex* d_work = NULL;

    checkCudaErrors(cusolverDnCreate(&cusolverH));

    checkCudaErrors(cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking));
    checkCudaErrors(cusolverDnSetStream(cusolverH, stream));

    checkCudaErrors(cudaMalloc(reinterpret_cast<void**>(&d_A), sizeof(QLComplex) * m_uiX * m_uiY));
    checkCudaErrors(cudaMalloc(reinterpret_cast<void**>(&d_B), sizeof(QLComplex) * m_uiX));
    checkCudaErrors(cudaMalloc(reinterpret_cast<void**>(&d_X), sizeof(QLComplex) * m_uiX));
    checkCudaErrors(cudaMalloc(reinterpret_cast<void**>(&d_info), sizeof(INT)));

    checkCudaErrors(cudaMemcpyAsync(d_A, HostBuffer(), sizeof(QLComplex) * m_uiX * m_uiY, cudaMemcpyHostToDevice, stream));
    checkCudaErrors(cudaMemcpyAsync(d_B, y.HostBuffer(), sizeof(QLComplex) * m_uiX, cudaMemcpyHostToDevice, stream));

#if _QL_DOUBLEFLOAT
    checkCudaErrors(cusolverDnZZgels_bufferSize(cusolverH, m, m, 1, d_A, m, d_B, m, d_X, m, d_work, &d_lwork));
#else
    checkCudaErrors(cusolverDnCCgels_bufferSize(cusolverH, m, m, 1, d_A, m, d_B, m, d_X, m, d_work, &d_lwork));
#endif

    checkCudaErrors(cudaMalloc(reinterpret_cast<void**>(&d_work), sizeof(QLComplex) * d_lwork));

    INT niters;
#if _QL_DOUBLEFLOAT
    checkCudaErrors(cusolverDnZZgels(cusolverH, m, m, 1, d_A, m, d_B, m, d_X, m, d_work, d_lwork, &niters, d_info));
#else
    checkCudaErrors(cusolverDnCCgels(cusolverH, m, m, 1, d_A, m, d_B, m, d_X, m, d_work, d_lwork, &niters, d_info));
#endif

    checkCudaErrors(cudaMemcpyAsync(&info, d_info, sizeof(INT), cudaMemcpyDeviceToHost, stream));
    checkCudaErrors(cudaStreamSynchronize(stream));

    if (0 > info) 
    {
        appCrucial(_T("%d-th parameter is wrong \n"), -info);
        return QLMatrix();
    }

    if (niters < 0)
    {
        appWarning(_T("iteration failed: %d\n"), niters);
    }

    QLComplex* hostX = reinterpret_cast<QLComplex*>(malloc(sizeof(QLComplex) * m_uiX));
    checkCudaErrors(cudaMemcpyAsync(hostX, d_X, sizeof(QLComplex) * m_uiX, cudaMemcpyDeviceToHost, stream));
    checkCudaErrors(cudaStreamSynchronize(stream));

    checkCudaErrors(cudaFree(d_A));
    checkCudaErrors(cudaFree(d_B));
    checkCudaErrors(cudaFree(d_X));
    checkCudaErrors(cudaFree(d_info));
    checkCudaErrors(cudaFree(d_work));

    checkCudaErrors(cusolverDnDestroy(cusolverH));

    checkCudaErrors(cudaStreamDestroy(stream));


    return QLMatrix(m_uiX, 1, hostX);
}

static QLComplex _hadamardmtr[4] = { 
    _make_cuComplex(InvSqrt2, F(0.0)), 
    _make_cuComplex(InvSqrt2, F(0.0)), 
    _make_cuComplex(InvSqrt2, F(0.0)), 
    _make_cuComplex(-InvSqrt2, F(0.0)) 
};

static QLComplex _PauliXmtr[4] = {
    _make_cuComplex(F(0.0), F(0.0)),
    _make_cuComplex(F(1.0), F(0.0)),
    _make_cuComplex(F(1.0), F(0.0)),
    _make_cuComplex(F(0.0), F(0.0))
};

static QLComplex _PauliYmtr[4] = {
    _make_cuComplex(F(0.0), F(0.0)),
    _make_cuComplex(F(0.0), F(1.0)),
    _make_cuComplex(F(0.0), F(-1.0)),
    _make_cuComplex(F(0.0), F(0.0))
};

static QLComplex _PauliZmtr[4] = {
    _make_cuComplex(F(1.0), F(0.0)),
    _make_cuComplex(F(0.0), F(0.0)),
    _make_cuComplex(F(0.0), F(0.0)),
    _make_cuComplex(F(-1.0), F(0.0))
};

static QLComplex _I2mtr[4] = {
    _make_cuComplex(F(1.0), F(0.0)),
    _make_cuComplex(F(0.0), F(0.0)),
    _make_cuComplex(F(0.0), F(0.0)),
    _make_cuComplex(F(1.0), F(0.0))
};

const QLAPI QLMatrix _hadamard = QLMatrix::CopyCreate(2U, 2U, _hadamardmtr);
const QLAPI QLMatrix _PauliX = QLMatrix::CopyCreate(2U, 2U, _PauliXmtr);
const QLAPI QLMatrix _PauliY = QLMatrix::CopyCreate(2U, 2U, _PauliYmtr);
const QLAPI QLMatrix _PauliZ = QLMatrix::CopyCreate(2U, 2U, _PauliZmtr);
const QLAPI QLMatrix _I2 = QLMatrix::CopyCreate(2U, 2U, _I2mtr);

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================