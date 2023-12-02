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

    QLGate subCnNot = CreateCnNot(static_cast<BYTE>(newController.Num()));

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


/**
* it is arranged as:
* c * n + b * (n - 2) - T
* where b is borrowed
* c is controlled
* T is target
* 
* nController = 5
* 0 ---------------o-----------------------o
*                  |                       |
* 1 ---------------o-----------------------o
*                  |                       |
* 5 -----------o---+---o---------------o---+---o
*              |       |               |       |
* 2 -----------o-------o---------------o-------o
*              |       |               |       |
* 6 -------o---+-------+---o-------o---+-------+---o
*          |               |       |               |
* 3 -------o---------------o-------o---------------o
*          |               |       |               |
* 7 ---o---+---------------+---o---+---------------+
*      |                       |
* 4 ---o-----------------------o
*      |                       |
* 8 ---+-----------------------+
* 
* 
*/
QLGate CreateCnNotWithNM2BorrowedAncilla(BYTE nControlledQubits)
{
    assert(nControlledQubits > 2);
    if (0 == nControlledQubits)
    {
        return QLGate(EBasicOperation::EBO_X);
    }
    if (1 == nControlledQubits)
    {
        return QLGate(EBasicOperation::EBO_CX);
    }
    if (2 == nControlledQubits)
    {
        return CreateToffoliGate();
    }

    //TArray<BYTE> allQubits;
    //allQubits.Append(ByteSequnce, 2 * nControlledQubits - 1);

    QLGate toffoli = CreateToffoliGate();
    TArray<BYTE> toAdd;
    toAdd.Append(ByteSequnce, 3);
    QLGate retGate;
    retGate.AddQubits(2 * nControlledQubits - 1);

    //Assume nControlledQubits = 5
    //c = 0,1,2,3,4  
    //b = 5,6,7
    //t = 8
    const BYTE target = 2 * nControlledQubits - 2; // t=8

    BYTE loop = nControlledQubits - 2; //loop = 3

    //it has n-2 b-c-T toffoli up
    toAdd[1] = target;
    for (BYTE start = 0; start < loop; ++start)
    {
        //borrowed start from target - 1,
        //controller start from n - 1,
        toAdd[0] = nControlledQubits - 1 - start;
        //target is always the last used ancilla
        toAdd[2] = toAdd[1];
        toAdd[1] = target - 1 - start;
        retGate.AppendGate(toffoli, toAdd);

        //first loop:  4,7 - 8
        //second loop: 3,6 - 7
        //third loop:  2,5 - 6
    }

    //and then 1 cc-T toffoli  0,1 - 5
    toAdd[0] = 0;
    toAdd[2] = toAdd[1];
    toAdd[1] = 1;
    retGate.AppendGate(toffoli, toAdd);

    //again n-1 b-c-T toffoli down
    for (BYTE start = 0; start < loop; ++start)
    {
        //borrowed start from n,
        //controller start from 2,
        toAdd[0] = 2 + start;
        toAdd[1] = nControlledQubits + start;
        //target is always the next ancilla except for the last one
        toAdd[2] = nControlledQubits + start + 1;

        retGate.AppendGate(toffoli, toAdd);

        //first loop:  2,5 - 6
        //second loop: 3,6 - 7
        //third loop:  4,7 - 8
    }

    //n-2 b-c-T toffoli up
    for (BYTE start = 1; start < loop; ++start)
    {
        //borrowed start from target - 1,
        //controller start from n - 1,
        toAdd[0] = nControlledQubits - 1 - start;
        toAdd[1] = target - 1 - start;
        toAdd[2] = target - start;
        retGate.AppendGate(toffoli, toAdd);

        //first loop:  3,6 - 7
        //second loop: 2,5 - 6
    }
    
    //and then 1 cc-T toffoli  0,1 - 5
    toAdd[0] = 0;
    toAdd[1] = 1;
    toAdd[2] = nControlledQubits;
    retGate.AppendGate(toffoli, toAdd);

    //again n-2 b-c-T toffoli down
    loop = loop - 1;
    for (BYTE start = 0; start < loop; ++start)
    {
        //borrowed start from n,
        //controller start from 2,
        toAdd[0] = 2 + start;
        toAdd[1] = nControlledQubits + start;
        toAdd[2] = nControlledQubits + start + 1;
        retGate.AppendGate(toffoli, toAdd);

        //first loop:  2,5 - 6
        //second loop: 3,6 - 7
    }

    return retGate;
}

/**
* divide into:
* 2/3/4 Cn/2-NOT based on type of eAT
* 
* If no-equal decompose can be archived, for
* For Zeroed, the first one is larger
* Otherwise, the second one is larger
* 
* 
*/
QLGate QLAPI CreateCnNotWithAncilla(BYTE numOfController, EAncillaType eAT)
{
    if (0 == numOfController)
    {
        return QLGate(EBasicOperation::EBO_X);
    }
    if (1 == numOfController)
    {
        return QLGate(EBasicOperation::EBO_CX);
    }
    if (2 == numOfController)
    {
        return CreateToffoliGate();
    }

    TArray<BYTE> firstOneQubits;
    TArray<BYTE> secondOneQubits;
    BYTE firstOneQubitCount = 0;
    BYTE secondOneQubitCount = 0;
    QLGate ret;
    ret.AddQubits(numOfController + 2);
    if (3 == numOfController)
    {
        //0 2 - A
        //1 A - T
        firstOneQubits.AddItem(0);
        firstOneQubits.AddItem(2);
        firstOneQubits.AddItem(4);

        secondOneQubits.AddItem(1);
        secondOneQubits.AddItem(4);
        secondOneQubits.AddItem(3);

        QLGate firstOne = CreateToffoliGate();
        QLGate secondOne = CreateToffoliGate();

        switch (eAT)
        {
        case EAncillaType::EAT_Zeroed:
            {
                ret.AppendGate(firstOne, firstOneQubits);
                ret.AppendGate(secondOne, secondOneQubits);
                ret.AppendGate(firstOne, firstOneQubits);
            }
            break;
        case EAncillaType::EAT_Garbage:
            {
                ret.AppendGate(secondOne, secondOneQubits);
                ret.AppendGate(firstOne, firstOneQubits);
                ret.AppendGate(secondOne, secondOneQubits);
            }
            break;
        case EAncillaType::EAT_Borrowed:
            {
                ret.AppendGate(firstOne, firstOneQubits);
                ret.AppendGate(secondOne, secondOneQubits);
                ret.AppendGate(firstOne, firstOneQubits);
                ret.AppendGate(secondOne, secondOneQubits);
            }
            break;
        case EAncillaType::EAT_Burnable:
        default:
            {
                ret.AppendGate(firstOne, firstOneQubits);
                ret.AppendGate(secondOne, secondOneQubits);
            }
            break;
        }
    }
    else if (4 == numOfController)
    {
        //one toffoli and the other is cn-NOT with borrowed ancilla
        if (EAncillaType::EAT_Zeroed == eAT)
        {
            //first toffoli 0 2 - A
            //second cn-NOT 1 3 A - T, ancilla:2
            firstOneQubits.AddItem(0);
            firstOneQubits.AddItem(2);
            firstOneQubits.AddItem(5);

            secondOneQubits.AddItem(1);
            secondOneQubits.AddItem(3);
            secondOneQubits.AddItem(5);
            secondOneQubits.AddItem(2);
            secondOneQubits.AddItem(4);

            QLGate firstOne = CreateToffoliGate();
            QLGate secondOne = CreateCnNotWithNM2BorrowedAncilla(3);

            ret.AppendGate(firstOne, firstOneQubits);
            ret.AppendGate(secondOne, secondOneQubits);
            ret.AppendGate(firstOne, firstOneQubits);
        }
        else
        {
            //first cn-NOT 0 2 3 - A, ancilla : 1
            //second toffoli 1 A - T
            firstOneQubits.AddItem(0);
            firstOneQubits.AddItem(2);
            firstOneQubits.AddItem(3);
            firstOneQubits.AddItem(1);
            firstOneQubits.AddItem(5);

            secondOneQubits.AddItem(1);
            secondOneQubits.AddItem(5);
            secondOneQubits.AddItem(4);

            QLGate firstOne = CreateCnNotWithNM2BorrowedAncilla(3);
            QLGate secondOne = CreateToffoliGate();

            switch (eAT)
            {
            case EAncillaType::EAT_Garbage:
                {
                    ret.AppendGate(secondOne, secondOneQubits);
                    ret.AppendGate(firstOne, firstOneQubits);
                    ret.AppendGate(secondOne, secondOneQubits);
                }
                break;
            case EAncillaType::EAT_Borrowed:
                {
                    ret.AppendGate(firstOne, firstOneQubits);
                    ret.AppendGate(secondOne, secondOneQubits);
                    ret.AppendGate(firstOne, firstOneQubits);
                    ret.AppendGate(secondOne, secondOneQubits);
                }
                break;
            case EAncillaType::EAT_Burnable:
                default:
                {
                    ret.AppendGate(firstOne, firstOneQubits);
                    ret.AppendGate(secondOne, secondOneQubits);
                }
                break;
            }
        }
    }
    else
    {
        if (EAncillaType::EAT_Zeroed == eAT)
        {
            //===============
            // the first one must target to A
            // the second one is as large as possible, but not more than first one - 2
            // the second one can use all controller of the first one
            // the second one cannot use T as ancilla
            // 
            // If n=5, 2,4
            //      6, 3,4
            //      7, 3,5
            //      8, 4,5
            // 
            // second one is (n+3)/2
            // first one is n-b2+1
            // 
            // controller:
            //  first:  0,1,..,b1-1,
            //  second: b1,b1+1,...,(numOfController - 1),A
            // 
            // ancilla:
            //  first:  first b1-2 of second,
            //  second: first b2-2 of first, since b2>=b1+2, always have enough ancilla

            if (5 == numOfController)
            {
                //T=5, A=6, C=0,1,2,3,4
                firstOneQubits.AddItem(0);
                firstOneQubits.AddItem(2);
                firstOneQubits.AddItem(6);

                secondOneQubits.AddItem(1);
                secondOneQubits.AddItem(3);
                secondOneQubits.AddItem(4);
                secondOneQubits.AddItem(6);

                secondOneQubits.AddItem(0);
                secondOneQubits.AddItem(2);

                secondOneQubits.AddItem(5);

                QLGate firstOne = CreateToffoliGate();
                QLGate secondOne = CreateCnNotWithNM2BorrowedAncilla(4);

                ret.AppendGate(firstOne, firstOneQubits);
                ret.AppendGate(secondOne, secondOneQubits);
                ret.AppendGate(firstOne, firstOneQubits);
            }
            else
            {
                BYTE b2 = (numOfController + 3) / 2;
                BYTE b1 = numOfController - b2 + 1;

                for (BYTE controllers = 0; controllers < numOfController; ++controllers)
                {
                    if (controllers < b1)
                    {
                        firstOneQubits.AddItem(controllers);
                    }
                    else
                    {
                        secondOneQubits.AddItem(controllers);
                    }
                }

                //T=numOfController, A=numOfController+1
                secondOneQubits.AddItem(numOfController + 1);
                for (BYTE anci = 0; anci < b1 - 2; ++anci)
                {
                    firstOneQubits.AddItem(secondOneQubits[anci]);
                }
                for (BYTE anci = 0; anci < b2 - 2; ++anci)
                {
                    secondOneQubits.AddItem(firstOneQubits[anci]);
                }

                firstOneQubits.AddItem(numOfController + 1);
                secondOneQubits.AddItem(numOfController);

                QLGate firstOne = CreateCnNotWithNM2BorrowedAncilla(b1);
                QLGate secondOne = CreateCnNotWithNM2BorrowedAncilla(b2);

                ret.AppendGate(firstOne, firstOneQubits);
                ret.AppendGate(secondOne, secondOneQubits);
                ret.AppendGate(firstOne, firstOneQubits);
            }
        }
        else if (EAncillaType::EAT_Garbage == eAT)
        {
            //===============
            // the first one must target to A and can use T as ancilla
            // the first one is as large as possible, but not more than second one - 2 (because T can be ancilla for the first one)
            // 
            // If n=5, 4,2
            //      6, 4,3
            //      7, 5,3
            //      8, 5,4

            if (5 == numOfController)
            {
                //0,1,2,3,4
                //0,1,2,3-A (4,T)
                //4,A-T
                firstOneQubits.AddItem(0);
                firstOneQubits.AddItem(1);
                firstOneQubits.AddItem(2);
                firstOneQubits.AddItem(3);

                firstOneQubits.AddItem(4);
                firstOneQubits.AddItem(5);

                firstOneQubits.AddItem(6);

                secondOneQubits.AddItem(4);
                secondOneQubits.AddItem(6);
                secondOneQubits.AddItem(5);

                QLGate firstOne = CreateCnNotWithNM2BorrowedAncilla(4); 
                QLGate secondOne = CreateToffoliGate();

                ret.AppendGate(secondOne, secondOneQubits);
                ret.AppendGate(firstOne, firstOneQubits);
                ret.AppendGate(secondOne, secondOneQubits);
            }
            else
            {
                BYTE b1 = (numOfController + 3) / 2;
                BYTE b2 = numOfController - b1 + 1;

                for (BYTE controllers = 0; controllers < numOfController; ++controllers)
                {
                    if (controllers < b1)
                    {
                        firstOneQubits.AddItem(controllers);
                    }
                    else
                    {
                        secondOneQubits.AddItem(controllers);
                    }
                }

                //T=numOfController, A=numOfController+1
                secondOneQubits.AddItem(numOfController + 1);
                for (BYTE anci = 0; anci < b1 - 3; ++anci)
                {
                    firstOneQubits.AddItem(secondOneQubits[anci]);
                }
                firstOneQubits.AddItem(numOfController);

                for (BYTE anci = 0; anci < b2 - 2; ++anci)
                {
                    secondOneQubits.AddItem(firstOneQubits[anci]);
                }

                firstOneQubits.AddItem(numOfController + 1);
                secondOneQubits.AddItem(numOfController);

                QLGate firstOne = CreateCnNotWithNM2BorrowedAncilla(b1);
                QLGate secondOne = CreateCnNotWithNM2BorrowedAncilla(b2);

                ret.AppendGate(secondOne, secondOneQubits);
                ret.AppendGate(firstOne, firstOneQubits);
                ret.AppendGate(secondOne, secondOneQubits);
            }
        }
        else
        {
            if (numOfController & 1)
            {
                //When T = 7
                //decompse into: 0 2 4 6 - A,  ancilla: 3,5
                //               1 3 5 A - T,  ancilla: 4,6
                //always add the last ((n - 1)/2) - 2 as the borrowed ancilla

                firstOneQubitCount = (numOfController + 1) / 2; //4
                secondOneQubitCount = firstOneQubitCount;
                for (BYTE start = 0; start < firstOneQubitCount - 1; ++start)
                {
                    firstOneQubits.AddItem(start * 2);       //0, 2, 4
                    secondOneQubits.AddItem(start * 2 + 1);  //1, 3, 5
                }
                firstOneQubits.AddItem((firstOneQubitCount - 1) * 2); // 6
                secondOneQubits.AddItem(numOfController + 1);         // A

                for (BYTE start = 0; start < firstOneQubitCount - 2; ++start)
                {
                    firstOneQubits.AddItem(secondOneQubits[1 + start]); //3, 5
                    secondOneQubits.AddItem(firstOneQubits[2 + start]); //4, 6
                }

                firstOneQubits.AddItem(numOfController + 1); // A
                secondOneQubits.AddItem(numOfController);    // T
            }
            else
            {
                //When T = 8
                //decompse into: 0 2 4 6 7 - A, ancilla: 1,3,5
                //               1 3 5 A   - T, ancilla: 4,6

                firstOneQubitCount = 1 + (numOfController / 2); //5
                secondOneQubitCount = firstOneQubitCount - 1;   //4
                for (BYTE start = 0; start < secondOneQubitCount - 1; ++start)
                {
                    firstOneQubits.AddItem(start * 2); //0,2,4
                    secondOneQubits.AddItem(start * 2 + 1); //1,3,5
                }
                firstOneQubits.AddItem(2 * secondOneQubitCount - 2); //6
                firstOneQubits.AddItem(numOfController - 1);         //7
                secondOneQubits.AddItem(numOfController + 1);        //A

                for (BYTE start = 0; start < secondOneQubitCount - 2; ++start)
                {
                    firstOneQubits.AddItem(secondOneQubits[start]);     //1,3
                    secondOneQubits.AddItem(firstOneQubits[2 + start]); //4,6
                }
                firstOneQubits.AddItem(secondOneQubits[secondOneQubitCount - 2]); //5

                firstOneQubits.AddItem(numOfController + 1); //A
                secondOneQubits.AddItem(numOfController);    //T
            }

            QLGate firstOne = CreateCnNotWithNM2BorrowedAncilla(firstOneQubitCount);
            QLGate secondOne = CreateCnNotWithNM2BorrowedAncilla(secondOneQubitCount);

            if (EAncillaType::EAT_Borrowed == eAT)
            {
                ret.AppendGate(firstOne, firstOneQubits);
                ret.AppendGate(secondOne, secondOneQubits);
                ret.AppendGate(firstOne, firstOneQubits);
                ret.AppendGate(secondOne, secondOneQubits);
            }
            else
            {
                ret.AppendGate(firstOne, firstOneQubits);
                ret.AppendGate(secondOne, secondOneQubits);
            }
        }
    }
    return ret;
}

/**
* -----o-------o-----
*      |       |
* -V-------U-------U-
*  |   |   |   |   |
* -o---+---o---+---o-
* 
* V=U^-2
* 
* Only Zeroed is calculated
*/
QLGate QLAPI CreateCnUWithAncilla(BYTE numOfController, const QLGate& cu, const QLGate& cui2, EAncillaType eAT)
{
    assert(numOfController > 2);
    if (numOfController <= 2)
    {
        appCrucial(_T("CreateCnUWithAncilla support for only controller >= 3!\n"));
        return QLGate();
    }

    if (EAncillaType::EAT_Burnable == eAT || EAncillaType::EAT_Garbage == eAT)
    {
        appCrucial(_T("CreateCnUWithAncilla only tested with Zeroed and Borrowed!\n"));
        return QLGate();
    }

    TArray<BYTE> cuqubits;
    cuqubits.AddItem(numOfController + 1);
    cuqubits.AddItem(numOfController);

    TArray<BYTE> cnotqubits;
    cnotqubits.Append(ByteSequnce, numOfController);
    cnotqubits.AddItem(numOfController + 1);
    cnotqubits.AddItem(numOfController);

    QLGate cnNOT = CreateCnNotWithAncilla(numOfController, EAncillaType::EAT_Borrowed);
    QLGate ret;
    ret.AddQubits(numOfController + 2);

    if (EAncillaType::EAT_Borrowed == eAT)
    {
        ret.AppendGate(cui2, cuqubits);
        ret.AppendGate(cnNOT, cnotqubits);
        ret.AppendGate(cu, cuqubits);
        ret.AppendGate(cnNOT, cnotqubits);
        ret.AppendGate(cu, cuqubits);
    }
    else
    {
        ret.AppendGate(cnNOT, cnotqubits);
        ret.AppendGate(cu, cuqubits);
        ret.AppendGate(cnNOT, cnotqubits);
    }

    return ret;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================