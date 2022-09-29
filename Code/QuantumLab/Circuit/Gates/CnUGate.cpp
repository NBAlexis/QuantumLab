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
    retGate.AppendGate(cv, czyzBits);

    QLGate subCnNot = CreateCnNot(numOfController - 1);

    retGate.AppendGate(subCnNot, controllerBits);

    cv.Dagger();
    retGate.AppendGate(cv, czyzBits);

    retGate.AppendGate(subCnNot, controllerBits);

    QLGate subCnU = CreateCnU(numOfController - 1, v);
    retGate.AppendGate(subCnU, subCnUBits);

    return retGate;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================