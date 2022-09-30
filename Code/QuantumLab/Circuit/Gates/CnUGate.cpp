//=============================================================================
// FILENAME : CnUGate.cpp
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [23/09/2022 nbale]
//=============================================================================

#include "QuantumLabPCH.h"

__BEGIN_NAMESPACE

void AppendCNOTSqrt(QLGate& gate, BYTE controller, BYTE target, UINT step, UBOOL bDagger)
{
    if (bDagger)
    {
        TArray<BYTE> bits;
        bits.AddItem(controller);
        bits.AddItem(target);

        TArray<BYTE> bits2;
        bits2.AddItem(controller);

        QLGate crx(EBasicOperation::EBO_CRX, -PI / (1ULL << step));
        QLGate p(EBasicOperation::EBO_P, -PI / (1ULL << (step + 1)));
        gate.AppendGate(p, bits2);
        gate.AppendGate(crx, bits);
    }
    else
    {
        TArray<BYTE> bits;
        bits.AddItem(controller);
        bits.AddItem(target);

        TArray<BYTE> bits2;
        bits2.AddItem(controller);

        QLGate crx(EBasicOperation::EBO_CRX, PI / (1ULL << step));
        QLGate p(EBasicOperation::EBO_P, PI / (1ULL << (step + 1)));
        gate.AppendGate(crx, bits);
        gate.AppendGate(p, bits2);
    }
}

void CNOTStep(QLGate& gate, const TArray<BYTE>& controller, BYTE target, UINT step)
{
    if (1 == controller.Num())
    {
        AppendCNOTSqrt(gate, controller[0], target, step, FALSE);
        return;
    }

    TArray<BYTE> newController = controller;
    BYTE lastController = newController.Pop();
    AppendCNOTSqrt(gate, lastController, target, step + 1, FALSE);

    QLGate subCnNot = CreateCnNot(controller.Num());

    gate.AppendGate(subCnNot, controller);

    AppendCNOTSqrt(gate, lastController, target, step + 1, TRUE);

    gate.AppendGate(subCnNot, controller);

    CNOTStep(gate, newController, target, step + 1);
}

QLGate QLAPI CreateCnNot(BYTE numOfController)
{
    if (1 == numOfController)
    {
        //this is a cnot
        return QLGate(EBasicOperation::EBO_CX);
    }

    QLGate retGate;
    retGate.AddQubits(numOfController + 1);
    retGate.m_sName = _T("CnNOT");

    TArray<BYTE> controllerBits;
    for (BYTE bit = 0; bit < numOfController; ++bit)
    {
        controllerBits.AddItem(bit);
    }
    TArray<BYTE> newControllerBits = controllerBits;
    newControllerBits.Pop();

    AppendCNOTSqrt(retGate, numOfController - 1, numOfController, 1, FALSE);

    QLGate subCnNot = CreateCnNot(numOfController - 1);

    retGate.AppendGate(subCnNot, controllerBits);

    AppendCNOTSqrt(retGate, numOfController - 1, numOfController, 1, TRUE);

    retGate.AppendGate(subCnNot, controllerBits);

    CNOTStep(retGate, newControllerBits, numOfController, 1);

    return retGate;
}

QLGate QLAPI CreateCnU(BYTE numOfController, const QLMatrix& mtr)
{
    if (1 == numOfController)
    {
        //this is a c-u
        return CreateControlledZYZGate(mtr);
    }

    QLMatrix v = Sqrt2by2(mtr);

    QLGate retGate;
    retGate.AddQubits(numOfController + 1);
    retGate.m_sName = _T("CnU");

    TArray<BYTE> controllerBits;
    for (BYTE bit = 0; bit < numOfController; ++bit)
    {
        controllerBits.AddItem(bit);
    }
    TArray<BYTE> newControllerBits = controllerBits;
    newControllerBits.Pop();
    TArray<BYTE> subCnUBits = newControllerBits;
    subCnUBits.AddItem(numOfController);

    TArray<BYTE> czyzBits;
    czyzBits.AddItem(numOfController - 1);
    czyzBits.AddItem(numOfController);

    QLGate cv = CreateControlledZYZGate(v);
    cv.m_sName = _T("CV");
    retGate.AppendGate(cv, czyzBits);

    QLGate subCnNot = CreateCnNot(numOfController - 1);
    retGate.AppendGate(subCnNot, controllerBits);

    cv.Dagger();
    cv.m_sName = _T("CVd");
    retGate.AppendGate(cv, czyzBits);

    retGate.AppendGate(subCnNot, controllerBits);

    QLGate subCnU = CreateCnU(numOfController - 1, v);
    subCnU.m_sName = _T("CnV");
    retGate.AppendGate(subCnU, subCnUBits);

    return retGate;
}

QLGate QLAPI CreateCnRX(BYTE numOfController, Real fDegree)
{
    if (1 == numOfController)
    {
        //this is a c-u
        return QLGate(EBasicOperation::EBO_CRX, fDegree);
    }

    QLGate retGate;
    retGate.AddQubits(numOfController + 1);
    retGate.m_sName = _T("CnRx");

    TArray<BYTE> controllerBits;
    for (BYTE bit = 0; bit < numOfController; ++bit)
    {
        controllerBits.AddItem(bit);
    }
    TArray<BYTE> newControllerBits = controllerBits;
    newControllerBits.Pop();
    TArray<BYTE> subCnUBits = newControllerBits;
    subCnUBits.AddItem(numOfController);

    TArray<BYTE> czyzBits;
    czyzBits.AddItem(numOfController - 1);
    czyzBits.AddItem(numOfController);

    QLGate cv = QLGate(EBasicOperation::EBO_CRX, fDegree / 2);
    retGate.AppendGate(cv, czyzBits);

    QLGate subCnNot = CreateCnNot(numOfController - 1);
    retGate.AppendGate(subCnNot, controllerBits);

    cv.Dagger();
    retGate.AppendGate(cv, czyzBits);

    retGate.AppendGate(subCnNot, controllerBits);

    QLGate subCnU = CreateCnRX(numOfController - 1, fDegree / 2);
    retGate.AppendGate(subCnU, subCnUBits);

    return retGate;
}

QLGate QLAPI CreateCnRY(BYTE numOfController, Real fDegree)
{
    if (1 == numOfController)
    {
        //this is a c-u
        return QLGate(EBasicOperation::EBO_CRY, fDegree);
    }

    QLGate retGate;
    retGate.AddQubits(numOfController + 1);
    retGate.m_sName = _T("CnRy");

    TArray<BYTE> controllerBits;
    for (BYTE bit = 0; bit < numOfController; ++bit)
    {
        controllerBits.AddItem(bit);
    }
    TArray<BYTE> newControllerBits = controllerBits;
    newControllerBits.Pop();
    TArray<BYTE> subCnUBits = newControllerBits;
    subCnUBits.AddItem(numOfController);

    TArray<BYTE> czyzBits;
    czyzBits.AddItem(numOfController - 1);
    czyzBits.AddItem(numOfController);

    QLGate cv = QLGate(EBasicOperation::EBO_CRY, fDegree / 2);
    retGate.AppendGate(cv, czyzBits);

    QLGate subCnNot = CreateCnNot(numOfController - 1);
    retGate.AppendGate(subCnNot, controllerBits);

    cv.Dagger();
    retGate.AppendGate(cv, czyzBits);

    retGate.AppendGate(subCnNot, controllerBits);

    QLGate subCnU = CreateCnRY(numOfController - 1, fDegree / 2);
    retGate.AppendGate(subCnU, subCnUBits);

    return retGate;
}

QLGate QLAPI CreateCnRZ(BYTE numOfController, Real fDegree)
{
    if (1 == numOfController)
    {
        //this is a c-u
        return QLGate(EBasicOperation::EBO_CRZ, fDegree);
    }

    QLGate retGate;
    retGate.AddQubits(numOfController + 1);
    retGate.m_sName = _T("CnRz");

    TArray<BYTE> controllerBits;
    for (BYTE bit = 0; bit < numOfController; ++bit)
    {
        controllerBits.AddItem(bit);
    }
    TArray<BYTE> newControllerBits = controllerBits;
    newControllerBits.Pop();
    TArray<BYTE> subCnUBits = newControllerBits;
    subCnUBits.AddItem(numOfController);

    TArray<BYTE> czyzBits;
    czyzBits.AddItem(numOfController - 1);
    czyzBits.AddItem(numOfController);

    QLGate cv = QLGate(EBasicOperation::EBO_CRZ, fDegree / 2);
    retGate.AppendGate(cv, czyzBits);

    QLGate subCnNot = CreateCnNot(numOfController - 1);
    retGate.AppendGate(subCnNot, controllerBits);

    cv.Dagger();
    retGate.AppendGate(cv, czyzBits);

    retGate.AppendGate(subCnNot, controllerBits);

    QLGate subCnU = CreateCnRZ(numOfController - 1, fDegree / 2);
    retGate.AppendGate(subCnU, subCnUBits);

    return retGate;
}

QLGate QLAPI CreateCnP(BYTE numOfController, Real fDegree)
{
    if (1 == numOfController)
    {
        //this is a c-u
        return QLGate(EBasicOperation::EBO_CP, fDegree);
    }

    QLGate retGate;
    retGate.AddQubits(numOfController + 1);
    retGate.m_sName = _T("CnP");

    TArray<BYTE> controllerBits;
    for (BYTE bit = 0; bit < numOfController; ++bit)
    {
        controllerBits.AddItem(bit);
    }
    TArray<BYTE> newControllerBits = controllerBits;
    newControllerBits.Pop();
    TArray<BYTE> subCnUBits = newControllerBits;
    subCnUBits.AddItem(numOfController);

    TArray<BYTE> czyzBits;
    czyzBits.AddItem(numOfController - 1);
    czyzBits.AddItem(numOfController);

    QLGate cv = QLGate(EBasicOperation::EBO_CP, fDegree / 2);
    retGate.AppendGate(cv, czyzBits);

    QLGate subCnNot = CreateCnNot(numOfController - 1);
    retGate.AppendGate(subCnNot, controllerBits);

    cv.Dagger();
    retGate.AppendGate(cv, czyzBits);

    retGate.AppendGate(subCnNot, controllerBits);

    QLGate subCnU = CreateCnP(numOfController - 1, fDegree / 2);
    retGate.AppendGate(subCnU, subCnUBits);

    return retGate;
}

QLGate QLAPI CreateCnPh(BYTE numOfController, Real fDegree)
{
    if (1 == numOfController)
    {
        //this is a c-u
        return QLGate(EBasicOperation::EBO_Phase, fDegree).CreateControlled();
    }

    QLGate retGate;
    retGate.AddQubits(numOfController + 1);
    retGate.m_sName = _T("CnPh");

    TArray<BYTE> controllerBits;
    for (BYTE bit = 0; bit < numOfController; ++bit)
    {
        controllerBits.AddItem(bit);
    }
    TArray<BYTE> newControllerBits = controllerBits;
    newControllerBits.Pop();
    TArray<BYTE> subCnUBits = newControllerBits;
    subCnUBits.AddItem(numOfController);

    TArray<BYTE> czyzBits;
    czyzBits.AddItem(numOfController - 1);
    czyzBits.AddItem(numOfController);

    QLGate cv = QLGate(EBasicOperation::EBO_Phase, fDegree / 2).CreateControlled();
    retGate.AppendGate(cv, czyzBits);

    QLGate subCnNot = CreateCnNot(numOfController - 1);
    retGate.AppendGate(subCnNot, controllerBits);

    cv.Dagger();
    retGate.AppendGate(cv, czyzBits);

    retGate.AppendGate(subCnNot, controllerBits);

    QLGate subCnU = CreateCnPh(numOfController - 1, fDegree / 2);
    retGate.AppendGate(subCnU, subCnUBits);

    return retGate;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================