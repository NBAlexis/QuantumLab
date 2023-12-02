//=============================================================================
// FILENAME : QLMatrix.h
// 
// DESCRIPTION:
// note that in cuda
// m = {1,2,3,4,5,6,7,8}
// is 1  5
//    2  6
//    3  7
//    4  8
//
// REVISION: [dd/mm/yy]
//  [12/09/2022 nbale]
//=============================================================================

#ifndef _QLMATRIX_H_
#define _QLMATRIX_H_

__BEGIN_NAMESPACE

typedef QLComplex(*complexfunc)(const QLComplex& c);
typedef QLComplex(*complexfuncTwo)(const QLComplex& c1, const QLComplex& c2);
typedef QLComplex(*complexfuncTwoR)(const QLComplex& c, Real r);

class QLMatrixData
{
    QLMatrixData()
        : m_nRefs(0)
        , m_pData(NULL)

    {

    }

    ~QLMatrixData()
    {
        appSafeFree(m_pData);
    }

private:

    std::atomic<INT> m_nRefs;

public:

    friend class QLMatrix;
    QLComplex* m_pData;
};

class QLAPI QLMatrix
{
private:

    QLMatrixData* m_pData;

    void OnChangeContent()
    {
        if (NULL != m_pData && m_pData->m_nRefs <= 1)
        {
            return;
        }

        if (NULL != m_pData && m_pData->m_nRefs > 1)
        {
            m_pData->m_nRefs = m_pData->m_nRefs - 1;
        }
        QLMatrixData* pData = new QLMatrixData();
        pData->m_nRefs = 1;
        if (NULL == m_pData)
        {
            pData->m_pData = (QLComplex*)calloc(m_uiX * m_uiY, sizeof(QLComplex));
        }
        else
        {
            pData->m_pData = (QLComplex*)malloc(m_uiX * m_uiY * sizeof(QLComplex));
            assert(NULL != pData->m_pData);
            if (NULL != pData->m_pData)
            {
                memcpy(pData->m_pData, m_pData->m_pData, m_uiX * m_uiY * sizeof(QLComplex));
            }
        }
        m_pData = pData;
    }

public:

    QLMatrix();

    QLMatrix(UINT uiX, UINT uiY);

    //do NOT release the buffer
    QLMatrix(UINT uiX, UINT uiY, QLComplex* buffer);

    QLMatrix(const QLMatrix& other);

    static QLMatrix CopyCreate(UINT uiX, UINT uiY, QLComplex* buffer);

    const QLMatrix& operator=(const QLMatrix& other)
    {
        if (m_uiX * m_uiY != other.m_uiX * other.m_uiY)
        {
            assert(m_pData != other.m_pData);
            if (NULL != m_pData)
            {
                if (m_pData->m_nRefs > 1)
                {
                    m_pData->m_nRefs = m_pData->m_nRefs - 1;
                }
                else
                {
                    appSafeDelete(m_pData);
                }
            }
            m_uiX = other.m_uiX;
            m_uiY = other.m_uiY;
            m_pData = other.m_pData;
            m_pData->m_nRefs = m_pData->m_nRefs + 1;
            return *this;
        }

        if (m_pData != other.m_pData)
        {
            if (NULL != m_pData)
            {
                if (m_pData->m_nRefs > 1)
                {
                    m_pData->m_nRefs = m_pData->m_nRefs - 1;
                }
                else
                {
                    appSafeDelete(m_pData);
                }
            }
            m_uiX = other.m_uiX;
            m_uiY = other.m_uiY;
            m_pData = other.m_pData;
            m_pData->m_nRefs = m_pData->m_nRefs + 1;
        }

        return *this;
    }

    virtual ~QLMatrix();

    QLComplex Get(INT x, INT y) const
    {
        return HostBuffer()[m_uiY * x + y];
    }

    void Set(INT x, INT y, const QLComplex& v)
    {
        OnChangeContent();
        m_pData->m_pData[m_uiY * x + y] = v;
    }

    //void SetHostBuffer(QLComplex* pBuffer)
    //{

    //}

    inline const QLComplex* HostBuffer() const
    {
        return (const QLComplex*)(m_pData->m_pData);
    }

    //random(-1, 1) + random(-1, 1) I
    void RandomOne();

    void RandomOneReal();

    //same as scipy, m^dagger m = 1
    void RandomUnitary();

    //M = Q.R, where Q^+ Q =1, R is upper traingular
    void QR(QLMatrix& q, QLMatrix& r) const;

    void ReShape(UINT uiX, UINT uiY)
    {
        if (uiX * uiY == m_uiX * m_uiY)
        {
            m_uiX = uiX;
            m_uiY = uiY;
        }
    }

protected:
    void _QR(QLMatrix& q, QLMatrix& r) const;

public:

    //M = u.s.v, where s is (upper) diagnal and in fact real, u^+ u = 1, v^+ v = 1. 
    //Note that, M = u.s.v instead of u.s.v^+
    //SVD only supports Y > X matrix
    void SVD(QLMatrix& u, QLMatrix& s, QLMatrix& v) const;

    //Same as SVD but using Jacob 
    //Note that, M = u.s.v^+ instead of u.s.v
    void SVDJ(QLMatrix& u, QLMatrix& s, QLMatrix& v) const;

    void Transpose(UBOOL bConjugate = FALSE);
    void Dagger() { Transpose(TRUE); }
    void Opposite();

    void Add(const QLMatrix& other);
    void Add(const QLComplex& other);
    void Add(const Real& other);

    void Sub(const QLMatrix& other);
    void Sub(const QLComplex& other);
    void Sub(const Real& other);

    //cublasOperation_t = CUBLAS_OP_N (none), CUBLAS_OP_T (transpose) CUBLAS_OP_C (dagger)
    QLMatrix Mul(const QLMatrix& other, cublasOperation_t left = CUBLAS_OP_N, cublasOperation_t right = CUBLAS_OP_N) const;
    void Mul(const QLComplex& other);
    void Mul(const Real& other);
    void Div(const QLComplex& other);
    void Div(const Real& other);

    void ElementAbs();
    void ElementExp();
    void ElementIExp();
    void ElementAbsSq();
    void ElementSqrt();
    void ElementMul(const QLMatrix& other);
    void ElementDiv(const QLMatrix& other);

    QLAPI friend QLMatrix operator+(const QLMatrix& m1, const QLMatrix& m2) { QLMatrix ret = m1; ret.Add(m2); return ret; }
    QLAPI friend QLMatrix operator+(const QLMatrix& m, const QLComplex& v) { QLMatrix ret = m; ret.Add(v); return ret; }
    QLAPI friend QLMatrix operator+(const QLMatrix& m, const Real& v) { QLMatrix ret = m; ret.Add(v); return ret; }
    QLAPI friend QLMatrix operator+(const QLComplex& v, const QLMatrix& m) { return m + v; }
    QLAPI friend QLMatrix operator+(const Real& v, const QLMatrix& m) { return m + v; }

    QLAPI friend QLMatrix operator-(const QLMatrix& m1, const QLMatrix& m2) { QLMatrix ret = m2; ret.Opposite(); return m1 + ret; }
    QLAPI friend QLMatrix operator-(const QLMatrix& m, const QLComplex& v) { return m + _make_cuComplex(-v.x, -v.y); }
    QLAPI friend QLMatrix operator-(const QLMatrix& m, const Real& v) { return m + (-v); }
    QLAPI friend QLMatrix operator-(const QLComplex& v, const QLMatrix& m) { QLMatrix ret = m; ret.Opposite(); return ret + v; }
    QLAPI friend QLMatrix operator-(const Real& v, const QLMatrix& m) { QLMatrix ret = m; ret.Opposite(); return ret + v; }

    QLAPI friend QLMatrix operator*(const QLMatrix& m1, const QLMatrix& m2) { return m1.Mul(m2); }
    QLAPI friend QLMatrix operator*(const QLMatrix& m, const QLComplex& v) { QLMatrix ret = m; ret.Mul(v); return ret; }
    QLAPI friend QLMatrix operator*(const QLMatrix& m, const Real& v) { QLMatrix ret = m; ret.Mul(v); return ret; }
    QLAPI friend QLMatrix operator*(const Real& v, const QLMatrix& m) { return m * v; }
    QLAPI friend QLMatrix operator*(const QLComplex& v, const QLMatrix& m) { return m * v; }

    QLAPI friend QLMatrix operator/(const QLMatrix& m, const QLComplex& v) { QLMatrix ret = m; ret.Div(v); return ret; }
    QLAPI friend QLMatrix operator/(const QLMatrix& m, const Real& v) { QLMatrix ret = m; ret.Div(v); return ret; }

    void Print(const CCString& sName = __EmptyString) const;

    QLMatrix GetBlock(UINT uiXStart, UINT uiXLen, UINT uiYStart, UINT uiYLen) const;
    void SetBlock(UINT uiXStart, UINT uiXLen, UINT uiYStart, UINT uiYLen, const QLComplex* content);
    QLMatrix GetDiag() const;
    QLMatrix SquareMatrixByAddOne() const;
    QLMatrix ExtendZeros(UINT uiNewX, UINT uiNewY) const;

    /**
    * type 0:
    * C = A 0
    *     0 B
    * 
    * type 1:
    * C = A B
    * 
    * type 2:
    * C = A
    *     B
    */
    QLMatrix BlockAdd(const QLMatrix& m, BYTE addtype) const;

    /**
    * 
    * M = | U1 0  | | C | V^+
    *     | 0  U2 | | S | 
    * C, S are in fact diagonal of real numbers.
    * if M is unitary, C^2+S^2 = 1
    * It assumes that: 
    * M = | Q1 |
    *     | Q2 |
    * where Q1, Q2 have m_uiY >= m_uiX, C and S has same size as Q1, Q2
    */
    void CSD2BY1(QLMatrix& u1, QLMatrix& u2, QLMatrix& c, QLMatrix& s, QLMatrix& v, UINT uiSep) const;

    /**
    * M = | U1 0  | | C -S| |V1  0|^+
    *     | 0  U2 | | S  C| |0  V2|
    * NOTE~!!! it only tested for square matrix, with half decompose!
    */
    void CSD(QLMatrix& u1, QLMatrix& u2, QLMatrix& c, QLMatrix& s, QLMatrix& v1, QLMatrix& v2, UINT uiXSep, UINT uiYSep) const;

    /**
    * explain the result of csd.
    * here u, v are simply 
    * U = | U1 0  | V = |V1  0|
    *     | 0  U2 |     |0  V2|
    * 
    * cs = | C | -S 0 0 |
    *      | 0 |  0 I 0 |
    *      --------------
    *      | S |  C 0 0 |
    *      | 0 |  0 0 I |
    * 
    */
    static void CombinedCSD(QLMatrix& u, QLMatrix& cs, QLMatrix& v, const QLMatrix& u1, const QLMatrix& u2, const QLMatrix& c, const QLMatrix& s, const QLMatrix& v1, const QLMatrix& v2, UINT uiXSep, UINT uiYSep);

    /**
    * M V = V L
    * M must be a nxn Hermitian matrix
    * where L is a nxn diagonal matrix, filled with eignvalues in ascending order
    * V is eigenvector
    */
    void EVD(QLMatrix& v, QLMatrix& l) const;

    QLMatrix KroneckerProduct(const QLMatrix& m2) const;

    /**
    * solve Ax=y
    */
    QLMatrix GELS(const QLMatrix& y) const;

protected:

    /**
    * func(M)
    * M must be a nxn Hermitian matrix
    * func must be a static device function
    * the function pointer must be on device, therefore this method is protected.
    */
    void MatrixFunction(complexfunc func);
    void MatrixFunctionTwoR(complexfuncTwoR func, Real r);
    void ElementWiseFunction(complexfunc func);
    void ElementWiseFunctionTwo(complexfuncTwo func, const QLMatrix& other);
    void ElementWiseFunctionTwoR(complexfuncTwoR func, Real r);

public:
    
    void MatrixIExp(Real t);
    void MatrixExp();
    void MatrixSqrt();

    static QLMatrix CreateEye(UINT uiX, UINT uiY);

protected:

    static inline INT Ceil(ULONGLONG a, ULONGLONG b)
    {
        const ULONGLONG c = a / b;
        return static_cast<INT>((a == b * c) ? c : (c + 1));
    }

    void GetDimDiag(dim3& blocks, dim3& threads) const
    {
        if (m_uiX < m_uiY)
        {
            UINT uiThreadX = _QL_LAUNCH_MAX_THREAD;
            UINT uiBlockX = Ceil(m_uiX, _QL_LAUNCH_MAX_THREAD);

            if (uiBlockX > 1)
            {
                uiThreadX = Ceil(m_uiX, uiBlockX);
                assert(uiThreadX <= _QL_LAUNCH_MAX_THREAD);
            }

            blocks = dim3(uiBlockX, 1, 1);
            threads = dim3(uiThreadX, 1, 1);
        }
        else
        {
            UINT uiThreadY = _QL_LAUNCH_MAX_THREAD;
            UINT uiBlockY = Ceil(m_uiY, _QL_LAUNCH_MAX_THREAD);

            if (uiBlockY > 1)
            {
                uiThreadY = Ceil(m_uiY, uiBlockY);
                assert(uiThreadY <= _QL_LAUNCH_MAX_THREAD);
            }

            blocks = dim3(uiBlockY, 1, 1);
            threads = dim3(uiThreadY, 1, 1);
        }
    }

    void GetDim(dim3& blocks, dim3& threads) const
    {
        UINT uiThreadX = (m_uiX >= _QL_LAUNCH_MAX_THREAD_SQRT) ? _QL_LAUNCH_MAX_THREAD_SQRT : m_uiX;
        UINT uiThreadY = (m_uiY >= _QL_LAUNCH_MAX_THREAD_SQRT) ? _QL_LAUNCH_MAX_THREAD_SQRT : m_uiY;
        UINT uiBlockX = Ceil(m_uiX, uiThreadX);
        UINT uiBlockY = Ceil(m_uiY, uiThreadY);

        if (uiBlockX > 1)
        {
            uiThreadX = Ceil(m_uiX, uiBlockX);
            assert(uiThreadX <= _QL_LAUNCH_MAX_THREAD_SQRT);
        }
        if (uiBlockY > 1)
        {
            uiThreadY = Ceil(m_uiY, uiBlockY);
            assert(uiThreadY <= _QL_LAUNCH_MAX_THREAD_SQRT);
        }

        blocks = dim3(uiBlockX, uiBlockY, 1);
        threads = dim3(uiThreadX, uiThreadY, 1);
    }
        
    UINT m_uiX;
    UINT m_uiY;

public:

    UBOOL operator==(const QLMatrix& other) const
    {
        appCrucial(_T("== not support for QLMatrix\n"));
        return FALSE;
    }

    UINT X() const { return m_uiX; }
    UINT Y() const { return m_uiY; }

    TArray<QLComplex> ToVector() const;
    TArray<Real> ToVectorRe() const;
    QLComplex VectorDot(const QLMatrix& other, UBOOL bConjL = TRUE, UBOOL bConjR = FALSE) const;
    Real Norm2() const
    {
        return sqrt(_cuCabsf(this->VectorDot(*this)));
    }
    QLMatrix VectorFFT(UBOOL bForward = TRUE) const;

    static QLComplex ReduceSum(QLComplex* deviceBuffer, UINT uiLen);
};

extern const QLAPI QLMatrix _hadamard;
extern const QLAPI QLMatrix _PauliX;
extern const QLAPI QLMatrix _PauliY;
extern const QLAPI QLMatrix _PauliZ;
extern const QLAPI QLMatrix _I2;


__END_NAMESPACE


#endif //#ifndef _QLMATRIX_H_

//=============================================================================
// END OF FILE
//=============================================================================