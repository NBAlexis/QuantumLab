//=============================================================================
// FILENAME : CSDDecompose.h
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [30/09/2022 nbale]
//=============================================================================

#ifndef _CSDDECOMPOSE_H_
#define _CSDDECOMPOSE_H_

__BEGIN_NAMESPACE

enum class ECSDMatrix : UINT
{
    NotInitialed,
    AMatrix,
    BlockedUMatrix,
    DiagonalUMatrix,
    FRZMatrix,
    PartialFRZMatrix,
    GlobalPhase,
};

class QLAPI CSDMatrix
{
public:

    ECSDMatrix m_eType;
    INT m_iLevel;
    INT m_iTotLevel;
    TArray<Real> m_lstDegreeContent;
    TArray<QLMatrix> m_lstSubMatrix;

    CSDMatrix()
        : m_eType(ECSDMatrix::NotInitialed)
        , m_iLevel(0)
        , m_iTotLevel(0)
    {

    }

    void SetAsUMatrix(INT level, INT totLevel, const TArray<QLMatrix>& u);
    void SetAsAMatrix(INT level, INT totLevel, const TArray<QLMatrix>& a);
    void SetAsRzMatrix(INT level, INT totLevel, const TArray<Real>& degrees);
    void SetAsGlobalPhase(Real fPhase, INT totLevel);
    void SetAsPartialRZ(INT level, INT totLevel, const TArray<Real>& degrees);
    void SetAsDiagonalUMatrix(INT level, INT totLevel, const TArray<Real>& degrees);

protected:

    QLMatrix PrintAsAMatrix() const;
    QLMatrix PrintAsBlockedUMatrix() const;
    QLMatrix PrintAsDiagonalUMatrix() const;
    QLMatrix PrintAsFRZMatrix() const;
    QLMatrix PrintAsPartialFRZMatrix() const;
    QLMatrix PrintAsGlobalPhase() const;

    void GateAsAMatrix(QLGate& gate) const;
    void GateAsFRZMatrix(QLGate& gate) const;
    void GateAsPartialFRZMatrix(QLGate& gate) const;
    void GateAsGlobalPhase(QLGate& gate) const;

public:

    QLMatrix PrintMe() const
    {
        switch (m_eType)
        {
        case ECSDMatrix::AMatrix:
            return PrintAsAMatrix();
        case ECSDMatrix::BlockedUMatrix:
            return PrintAsBlockedUMatrix();
        case ECSDMatrix::DiagonalUMatrix:
            return PrintAsDiagonalUMatrix();
        case ECSDMatrix::FRZMatrix:
            return PrintAsFRZMatrix();
        case ECSDMatrix::PartialFRZMatrix:
            return PrintAsPartialFRZMatrix();
        case ECSDMatrix::GlobalPhase:
            return PrintAsGlobalPhase();
        }
        appCrucial(_T("something wrong! PrintMe\n"));
        return QLMatrix();
    }

    void CombinedYZGate(QLGate& gate, const TArray<Real>& zDegrees) const;

    void Gate(QLGate& gate) const
    {
        switch (m_eType)
        {
        case ECSDMatrix::AMatrix:
            return GateAsAMatrix(gate);
        case ECSDMatrix::FRZMatrix:
            return GateAsFRZMatrix(gate);
        case ECSDMatrix::PartialFRZMatrix:
            return GateAsPartialFRZMatrix(gate);
        case ECSDMatrix::GlobalPhase:
            return GateAsGlobalPhase(gate);
        }
        appCrucial(_T("something wrong! Gate\n"));
    }

    TArray<CSDMatrix> DecomposeU() const;

    /**
    * iLevel is the level of A matrix to commute
    */
    TArray<Real> CalculateUDegrees(INT iLevel);

    /**
    * iLevel is the level of A matrix to commute
    */
    void PhaseShiftU(const TArray<Real>& toShift, INT iLevel);

    TArray<CSDMatrix> DecomposeUIntoRZ() const;

    static TArray<CSDMatrix> RecursiveUDecompose(const TArray<CSDMatrix>& mtrList);

    static void PhaseShiftUAll(TArray<CSDMatrix>& mtrList);

    /**
    * iLevel is number of qubits
    */
    static TArray<CSDMatrix> CSDecomposeMatrix(const QLMatrix& u, INT iLevel);
};

extern QLGate QLAPI CSDDecompose(const QLMatrix& u, INT iLevel);

__END_NAMESPACE


#endif //#ifndef _CSDDECOMPOSE_H_

//=============================================================================
// END OF FILE
//=============================================================================