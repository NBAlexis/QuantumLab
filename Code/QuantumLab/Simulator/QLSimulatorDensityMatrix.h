//=============================================================================
// FILENAME : QLSimulatorDensityMatrix.h
// 
// DESCRIPTION:
// This is the file for simulation with density matrix (with noise)
//
// REVISION: [dd/mm/yy]
//  [03/01/2024 nbale]
//=============================================================================

#ifndef _QLSIMULATORDENSITYMATRIX_H_
#define _QLSIMULATORDENSITYMATRIX_H_

__BEGIN_NAMESPACE

class QLAPI QLSimulatorParametersDensityMatrix : public QLSimulatorParameters
{
public:

    QLSimulatorParametersDensityMatrix()
        : m_MasterGate()
        , m_byQubitCount(0)
        , m_bPrint(TRUE)

        , m_fDampingAfterGate(F(0.0))
        , m_fDepolarisingAfterGate(F(0.0))
        , m_fDephaseAfterGate(F(0.0))
        , m_fMixPauliAfterGate(F(0.0))
        , m_fTwoDepolarisingAfterGate(F(0.0))
        , m_fTwoDephaseAfterGate(F(0.0))

        , m_fDampingBeforeMeasure(F(0.0))
        , m_fDepolarisingBeforeMeasure(F(0.0))
        , m_fMixPauliBeforeMeasure(F(0.0))
        , m_iMeasureTimes(-1)
        , m_bMeasurePurity(FALSE)
        , m_bMeasureFidelity(FALSE)
    {
        
    }

    UINT BuildZeroStart(BYTE byQubit, Real* pReal, Real* pImag) const override
    {
        UINT uiVectorLength = 1U << (byQubit << 1U);

        memset(pReal, 0, sizeof(Real) * uiVectorLength);
        memset(pImag, 0, sizeof(Real) * uiVectorLength);
        pReal[0] = F(1.0);
        return uiVectorLength;
    }

    //to be changed to 'circuit' which including measurement (if measurement can be viewed as matrix)
    QLGate m_MasterGate;
    BYTE m_byQubitCount;
    UBOOL m_bPrint;

    Real m_fDampingAfterGate;
    Real m_fDepolarisingAfterGate;
    Real m_fDephaseAfterGate;
    Real m_fMixPauliAfterGate;
    Real m_fTwoDepolarisingAfterGate;
    Real m_fTwoDephaseAfterGate;

    Real m_fDampingBeforeMeasure;
    Real m_fDepolarisingBeforeMeasure;
    Real m_fDephaseBeforeMeasure;
    Real m_fMixPauliBeforeMeasure;

    TArray<BYTE> m_lstMeasureQubits;

    /**
    * if set to -1, calculate exact result
    */
    INT m_iMeasureTimes;

    UBOOL m_bMeasurePurity;
    UBOOL m_bMeasureFidelity;
};

class QLAPI QLSimulatorOutputDensityMatrix : public QLSimulatorOutput
{
public:

    /**
    * if the number measure qubits is 3,
    * it contains the probabilities of:
    * 000, 001, 010, 011, ..., 111
    * 
    */
    TArray<Real> m_lstMeasureOutcomes;

    /**
    * If measure purity was true, calculate purity
    */
    TArray<Real> m_lstPurity;

    TArray<Real> m_lstFidelity;

    Real AveragePurity() const
    {
        if (0 == m_lstPurity.Num())
        {
            return F(0.0);
        }
        if (1 == m_lstPurity.Num())
        {
            return m_lstPurity[0];
        }
        Real fsum = F(0.0);
        for (INT i = 0; i < m_lstPurity.Num(); ++i)
        {
            fsum += m_lstPurity[i];
        }
        return fsum / static_cast<Real>(m_lstPurity.Num());
    }

    Real AverageFidelity() const
    {
        if (0 == m_lstFidelity.Num())
        {
            return F(0.0);
        }
        if (1 == m_lstFidelity.Num())
        {
            return m_lstFidelity[0];
        }
        Real fsum = F(0.0);
        for (INT i = 0; i < m_lstFidelity.Num(); ++i)
        {
            fsum += m_lstFidelity[i];
        }
        return fsum / static_cast<Real>(m_lstFidelity.Num());
    }
};

class QLAPI QLSimulatorDensityMatrix : public QLSimulator
{

protected:

public:

    QLSimulatorDensityMatrix()
        : QLSimulator()
    {

    }

    void Simulate(QLSimulatorParameters* params, QLSimulatorOutput* output = NULL) const override;

protected:

    static void MeasurePurity(const QLSimulatorParametersDensityMatrix* params, QLSimulatorOutputDensityMatrix* output, const Qureg* v1);
    static void MeasureFidelity(const QLSimulatorParametersDensityMatrix* params, QLSimulatorOutputDensityMatrix* output, const Qureg* v1, const Qureg* v2);
};


__END_NAMESPACE


#endif //#ifndef _QLSIMULATORDENSITYMATRIX_H_

//=============================================================================
// END OF FILE
//=============================================================================