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

#pragma region kernels

__global__ void _QL_LAUNCH_BOUND
_kernalRandomOne(QLComplex* pDevicePtr, UINT perthread, UINT uiMax)
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

#pragma endregion

QLMatrix::QLMatrix()
    : m_uiX(0)
    , m_uiY(0)
    , m_pHostBuffer(NULL)
{

}

QLMatrix::QLMatrix(UINT uiX, UINT uiY)
    : m_uiX(uiX)
    , m_uiY(uiY)
{
    m_pHostBuffer = (QLComplex*)calloc(m_uiX * m_uiY, sizeof(QLComplex));
}

QLMatrix::QLMatrix(UINT uiX, UINT uiY, QLComplex* buffer)
    : m_uiX(uiX)
    , m_uiY(uiY)
{
    m_pHostBuffer = (QLComplex*)malloc(m_uiX * m_uiY * sizeof(QLComplex));
    if (NULL == m_pHostBuffer)
    {
        return;
    }
    memcpy(m_pHostBuffer, buffer, m_uiX * m_uiY * sizeof(QLComplex));
}

QLMatrix::QLMatrix(UINT uiX, UINT uiY, QLComplex** theArray)
    : m_uiX(uiX)
    , m_uiY(uiY)
{
    m_pHostBuffer = (QLComplex*)malloc(m_uiX * m_uiY * sizeof(QLComplex));
    if (NULL == m_pHostBuffer)
    {
        return;
    }
    for (UINT i = 0; i < m_uiX; ++i)
    {
        memcpy(m_pHostBuffer + m_uiY * i, theArray[i], m_uiX * sizeof(QLComplex));
    }
}

QLMatrix::~QLMatrix()
{
    appSafeFree(m_pHostBuffer);
}

void QLMatrix::RandomOne()
{
    UINT uiMax = m_uiX * m_uiY;
    UINT uiPerthread = 1 + uiMax / _QL_LAUNCH_MAX_THREAD;
    QLComplex* deviceBuffer = NULL;

    checkCudaErrors(cudaMalloc((void**)&deviceBuffer, sizeof(QLComplex) * uiMax));

    _kernalRandomOne << <1, _QL_LAUNCH_MAX_THREAD >> > (deviceBuffer, uiPerthread, uiMax);

    checkCudaErrors(cudaMemcpy(m_pHostBuffer, deviceBuffer, sizeof(QLComplex) * uiMax, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaFree(deviceBuffer));
}

CCString QLMatrix::Print() const
{
    CCString ret;
    ret += "{\n";
    for (UINT y = 0; y < m_uiY; ++y)
    {
        for (UINT x = 0; x < m_uiX; ++x)
        {
            if (0 == x)
            {
                ret += "{";
            }
            else
            {
                ret += ", ";
            }
            INT xy[2] = { static_cast<INT>(x), static_cast<INT>(y) };
            QLComplex toPrint = (*this)[xy];
            ret += appPrintComplex(toPrint.x, toPrint.y);
        }
        if (y == (m_uiY - 1))
        {
            ret += "}\n}\n";
        }
        else
        {
            ret += "},\n";
        }
    }
    return ret;
}

void QLMatrix::QR(QLMatrix** q, QLMatrix** r)
{
    cusolverDnHandle_t cusolverH = NULL;
    cublasHandle_t cublasH = NULL;
    cudaStream_t stream{};

    //======================================================================================================
    //from https://github.com/NVIDIA/CUDALibrarySamples/blob/master/cuSOLVER/orgqr/cusolver_orgqr_example.cu

    const INT m = m_uiY;
    const INT n = m_uiX;

    QLComplex* pQHost = reinterpret_cast<QLComplex*>(malloc(sizeof(QLComplex) * m * n));
    QLComplex* pRHost = reinterpret_cast<QLComplex*>(calloc(n * n, sizeof(QLComplex)));
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

    checkCudaErrors(cudaMemcpyAsync(d_A, m_pHostBuffer, sizeof(QLComplex) * m * n, cudaMemcpyHostToDevice, stream));

    /* step 3: query working space of geqrf and orgqr */
#if _QL_DOUBLEFLOAT
    checkCudaErrors(cusolverDnZgeqrf_bufferSize(cusolverH, m, n, d_A, m, &lwork_geqrf));
    checkCudaErrors(cusolverDnZungqr_bufferSize(cusolverH, m, n, n, d_A, m, d_tau, &lwork_orgqr));
#else
    checkCudaErrors(cusolverDnCgeqrf_bufferSize(cusolverH, m, n, d_A, m, &lwork_geqrf));
    checkCudaErrors(cusolverDnCungqr_bufferSize(cusolverH, m, n, n, d_A, m, d_tau, &lwork_orgqr));
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
        checkCudaErrors(cudaMemcpyAsync(pRHost + i * n, d_A + i * m, sizeof(QLComplex) * (i + 1), cudaMemcpyDeviceToHost, stream));
    }

    /* step 5: compute Q */
#if _QL_DOUBLEFLOAT
    checkCudaErrors(cusolverDnZungqr(cusolverH, m, n, n, d_A, m, d_tau, d_work, lwork, d_info));
#else
    checkCudaErrors(cusolverDnCungqr(cusolverH, m, n, n, d_A, m, d_tau, d_work, lwork, d_info));
#endif

    /* check if QR is good or not */
    checkCudaErrors(cudaMemcpyAsync(&info, d_info, sizeof(INT), cudaMemcpyDeviceToHost, stream));

    checkCudaErrors(cudaStreamSynchronize(stream));

    if (0 > info) 
    {
        appCrucial("%d-th parameter is wrong \n", -info);
    }

    checkCudaErrors(cudaMemcpyAsync(pQHost, d_A, sizeof(QLComplex) * m * n, cudaMemcpyDeviceToHost, stream));

    checkCudaErrors(cudaStreamSynchronize(stream));

    //std::printf("Q = (matlab base-1)\n");
    //print_matrix(m, n, Q.data(), lda);


    *q = new QLMatrix(n, m, pQHost);
    *r = new QLMatrix(n, n, pRHost);
    
    checkCudaErrors(cudaStreamSynchronize(stream));

    /* free resources */
    checkCudaErrors(cudaFree(d_A));
    checkCudaErrors(cudaFree(d_tau));
    checkCudaErrors(cudaFree(d_info));
    checkCudaErrors(cudaFree(d_work));

    checkCudaErrors(cublasDestroy(cublasH));
    checkCudaErrors(cusolverDnDestroy(cusolverH));
    checkCudaErrors(cudaStreamDestroy(stream));

    appSafeFree(pQHost);
    appSafeFree(pRHost);
}

void QLMatrix::RandomUnitary()
{
    RandomOne();
    QLMatrix* q;
    QLMatrix* r;

    QR(&q, &r);

    memcpy(m_pHostBuffer, q->m_pHostBuffer, sizeof(QLComplex) * m_uiX * m_uiY);

    appSafeDelete(q);
    appSafeDelete(r);
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================