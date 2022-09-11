//=============================================================================
// FILENAME : OtherComplexFunction.h
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [12/09/2022 nbale]
//=============================================================================

#ifndef _OTHER_COMPLEX_FUNCTION_H_
#define _OTHER_COMPLEX_FUNCTION_H_

__BEGIN_NAMESPACE

#if defined(__cplusplus)
extern "C" {
#endif /* __cplusplus */

    /**
    * arg(c)
    */
    __device__ __host__ static __inline__ Real __cuCargf(const QLComplex& c)
    {
        return _atan2(_cuCimagf(c), _cuCrealf(c));
    }

    /**
    * |c|^2
    */
    __device__ __host__ static __inline__ Real __cuCabsSqf(const QLComplex& c)
    {
        return c.x * c.x + c.y * c.y;
    }

#if !_QL_DOUBLEFLOAT
    __device__ __host__ static __inline__ DOUBLE __cuCabsSqfd(const cuDoubleComplex& c)
    {
        return c.x * c.x + c.y * c.y;
    }

    __device__ __host__ static __inline__ DOUBLE __cuCabsSqd(const cuDoubleComplex& c)
    {
        return c.x * c.x + c.y * c.y;
    }
#endif

    /**
    * c^p
    */
    __device__ static __inline__ QLComplex __cuCpowerf(const QLComplex& c, Real p)
    {
        const Real fArg = __cuCargf(c) * p;
        const Real fAbs = _pow(_cuCabsf(c), p);
        return _make_cuComplex(_cos(fArg) * fAbs, _sin(fArg) * fAbs);
    }

    /**
    * exp(c)
    */
    __device__ static __inline__ QLComplex __cuCexpf(const QLComplex& c)
    {
        const Real factor = _exp(c.x);
        return _make_cuComplex(factor * _cos(c.y), factor * _sin(c.y));
    }

    /**
     * log(c)
     */
    __device__ static __inline__ QLComplex __cuClogf(const QLComplex& c)
    {
        const Real fArg = __cuCargf(c);
        return _make_cuComplex(_log(_cuCabsf(c)), fArg > PI ? fArg - PI2 : fArg);
    }

    /**
    * _sqrt(c)
    */
    __device__ static __inline__ QLComplex __cuCsqrtf(const QLComplex& c)
    {
        const Real fRadius = _cuCabsf(c);
        const Real fCosA = __div(c.x, fRadius);
        QLComplex out;
        out.x = _sqrt(F(0.5) * fRadius * (fCosA + F(1.0)));
        out.y = _sqrt(F(0.5) * fRadius * (F(1.0) - fCosA));
        // signbit should be false if x.y is negative
        //if (signbit(c.y))
        //    out.y *= -F(1.0);
        if (c.y < F(0.0)) //same as Mathematica
            out.y *= -F(1.0);

        return out;
    }

    __device__ static __inline__ QLComplex cuCaddf_cr(const QLComplex& x, Real y)
    {
        return _make_cuComplex(x.x + y, x.y);
    }
    __device__ static __inline__ QLComplex cuCaddf_rc(Real y, const QLComplex& x)
    {
        return _make_cuComplex(x.x + y, x.y);
    }

    __device__ static __inline__ QLComplex cuCdivf_cr(const QLComplex& x, Real y)
    {
        return _make_cuComplex(__div(x.x, y), __div(x.y, y));
    }
#if !_QL_DOUBLEFLOAT
    __host__ static __inline__ cuDoubleComplex cuCdivf_cd_host(const cuDoubleComplex& x, DOUBLE y)
    {
        return make_cuDoubleComplex(x.x / y, x.y / y);
    }
    __device__ static __inline__ cuDoubleComplex cuCdivf_cd(const cuDoubleComplex& x, DOUBLE y)
    {
        return make_cuDoubleComplex(__div(x.x, y), __div(x.y, y));
    }
#endif

    __host__ static __inline__ QLComplex cuCdivf_cr_host(const QLComplex& x, Real y)
    {
        return _make_cuComplex(x.x / y, x.y / y);
    }

    __device__ __host__ static __inline__ QLComplex cuCmulf_cr(const QLComplex& x, Real y)
    {
        return _make_cuComplex(x.x * y, x.y * y);
    }

#if !_QL_DOUBLEFLOAT

    __device__ __host__ static __inline__ cuDoubleComplex cuCmulf_cd(const cuDoubleComplex& x, DOUBLE y)
    {
        return make_cuDoubleComplex(x.x * y, x.y * y);
    }
#endif

    __device__ static __inline__ QLComplex cuCmulf_rc(Real y, const QLComplex& x)
    {
        return _make_cuComplex(x.x * y, x.y * y);
    }

    __device__ static __inline__ QLComplex cuCsubf_cr(const QLComplex& x, Real y)
    {
        return _make_cuComplex(x.x - y, x.y);
    }
    __device__ static __inline__ QLComplex cuCsubf_rc(Real y, const QLComplex& x)
    {
        return _make_cuComplex(y - x.x, x.y);
    }

#if defined(__cplusplus)
}

__host__ __device__ __inline__
UBOOL operator==(const cuComplex& c1, const cuComplex& c2)
{
    return c1.x == c2.x && c1.y == c2.y;
}

__host__ __device__ __inline__
UBOOL operator==(const cuDoubleComplex& c1, const cuDoubleComplex& c2)
{
    return c1.x == c2.x && c1.y == c2.y;
}

__host__ __device__ __inline__
UBOOL operator!=(const cuComplex& c1, const cuComplex& c2)
{
    return !(c1 == c2);
}

__host__ __device__ __inline__
UBOOL operator!=(const cuDoubleComplex& c1, const cuDoubleComplex& c2)
{
    return !(c1 == c2);
}

#endif /* __cplusplus */

__END_NAMESPACE

#endif //#ifndef _OTHER_COMPLEX_FUNCTION_H_

//=============================================================================
// END OF FILE
//=============================================================================