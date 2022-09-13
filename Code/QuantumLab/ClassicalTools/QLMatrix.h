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

//TODO MAKE DATA REF
class QLAPI QLMatrix
{

public:

    QLMatrix();

    QLMatrix(UINT uiX, UINT uiY);

    QLMatrix(UINT uiX, UINT uiY, QLComplex* buffer);

    QLMatrix(UINT uiX, UINT uiY, QLComplex** theArray);

    virtual ~QLMatrix();

    QLComplex Get(INT x, INT y) const
    {
        return m_pHostBuffer[m_uiY * x + y];
    }

    void Set(INT x, INT y, const QLComplex& v)
    {
        m_pHostBuffer[m_uiY * x + y] = v;
    }

    inline const QLComplex& operator[](INT xy[]) const
    {
        return m_pHostBuffer[m_uiY * xy[0] + xy[1]];
    }

    inline QLComplex& operator[](INT xy[])
    {
        return m_pHostBuffer[m_uiY * xy[0] + xy[1]];
    }

    inline const QLComplex& operator[](BYTE xy[]) const
    {
        INT ixy[2] = {xy[0], xy[1]};
        return (*this)[ixy];
    }

    inline QLComplex& operator[](BYTE xy[])
    {
        INT ixy[2] = { xy[0], xy[1] };
        return (*this)[ixy];
    }

    inline const QLComplex& operator[](UINT xy[]) const
    {
        INT ixy[2] = { static_cast<INT>(xy[0]), static_cast<INT>(xy[1]) };
        return (*this)[ixy];
    }

    inline QLComplex& operator[](UINT xy[])
    {
        INT ixy[2] = { static_cast<INT>(xy[0]), static_cast<INT>(xy[1]) };
        return (*this)[ixy];
    }

    inline const QLComplex& operator[](LONGLONG xy[]) const
    {
        INT ixy[2] = { static_cast<INT>(xy[0]), static_cast<INT>(xy[1]) };
        return (*this)[ixy];
    }

    inline QLComplex& operator[](LONGLONG xy[])
    {
        INT ixy[2] = { static_cast<INT>(xy[0]), static_cast<INT>(xy[1]) };
        return (*this)[ixy];
    }

    inline const QLComplex* GetData() const
    {
        return (const QLComplex*)(m_pHostBuffer);
    }

    inline QLComplex* GetData()
    {
        return (QLComplex*)(m_pHostBuffer);
    }

    //random(-1, 1) + random(-1, 1) I
    void RandomOne();

    //same as scipy, m^dagger m = 1
    void RandomUnitary();

    void QR(QLMatrix** q, QLMatrix** r);

    CCString Print() const;

protected:

    UINT m_uiX;
    UINT m_uiY;
    QLComplex* m_pHostBuffer;
};




__END_NAMESPACE


#endif //#ifndef _QLMATRIX_H_

//=============================================================================
// END OF FILE
//=============================================================================