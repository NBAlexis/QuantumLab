//=============================================================================
// FILENAME : QLGate.cpp
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION: [dd/mm/yy]
//  [09/09/2022 nbale]
//=============================================================================

#include "QuantumLabPCH.h"

__BEGIN_NAMESPACE

QLGate::QLGate(EBasicOperation eOp, Real fParam)
    : m_eOp(eOp)
	, m_bDagger(FALSE)
	, m_fClassicalParameter(fParam)
{
	SBasicOperationInGate op_h;
	op_h.m_eOperation = eOp;

	switch (eOp)
	{
	case EBasicOperation::EBO_H:
		m_lstQubits.AddItem(0);
		m_lstOperations.AddItem(op_h);
		m_sName = "H";
		break;
	case EBasicOperation::EBO_X:
		m_lstQubits.AddItem(0);
		m_lstOperations.AddItem(op_h);
		m_sName = "X";
		break;
	case EBasicOperation::EBO_Y:
		m_lstQubits.AddItem(0);
		m_lstOperations.AddItem(op_h);
		m_sName = "Y";
		break;
	case EBasicOperation::EBO_Z:
		m_lstQubits.AddItem(0);
		m_lstOperations.AddItem(op_h);
		m_sName = "Z";
		break;
	case EBasicOperation::EBO_P:
		m_lstQubits.AddItem(0);
		m_lstOperations.AddItem(op_h);
		m_sName = "p";
		break;
	case EBasicOperation::EBO_Phase:
	{
		m_lstQubits.AddItem(0);
		SBasicOperationInGate op_1;
		op_1.m_eOperation = EBasicOperation::EBO_P;
		SBasicOperationInGate op_2;
		op_2.m_eOperation = EBasicOperation::EBO_X;
		SBasicOperationInGate op_3;
		op_3.m_eOperation = EBasicOperation::EBO_P;
		SBasicOperationInGate op_4;
		op_4.m_eOperation = EBasicOperation::EBO_X;
		m_lstOperations.AddItem(op_1);
		m_lstOperations.AddItem(op_2);
		m_lstOperations.AddItem(op_3);
		m_lstOperations.AddItem(op_4);
		m_sName = "phase";
	}
		break;

	case EBasicOperation::EBO_RX:
		m_lstQubits.AddItem(0);
		m_lstOperations.AddItem(op_h);
		m_sName = "rx";
		break;
	case EBasicOperation::EBO_RY:
		m_lstQubits.AddItem(0);
		m_lstOperations.AddItem(op_h);
		m_sName = "ry";
		break;
	case EBasicOperation::EBO_RZ:
		m_lstQubits.AddItem(0);
		m_lstOperations.AddItem(op_h);
		m_sName = "rz";
		break;

	case EBasicOperation::EBO_CX:
		m_lstQubits.AddItem(0);
		m_lstQubits.AddItem(1);
		m_lstOperations.AddItem(op_h);
		m_sName = "CNOT";
		break;
	case EBasicOperation::EBO_CY:
		m_lstQubits.AddItem(0);
		m_lstQubits.AddItem(1);
		m_lstOperations.AddItem(op_h);
		m_sName = "CY";
		break;
	case EBasicOperation::EBO_CZ:
		m_lstQubits.AddItem(0);
		m_lstQubits.AddItem(1);
		m_lstOperations.AddItem(op_h);
		m_sName = "CZ";
		break;

	case EBasicOperation::EBO_CP:
		m_lstQubits.AddItem(0);
		m_lstQubits.AddItem(1);
		m_lstOperations.AddItem(op_h);
		m_sName = "cp";
		break;
	case EBasicOperation::EBO_CRX:
		m_lstQubits.AddItem(0);
		m_lstQubits.AddItem(1);
		m_lstOperations.AddItem(op_h);
		m_sName = "crx";
		break;
	case EBasicOperation::EBO_CRY:
		m_lstQubits.AddItem(0);
		m_lstQubits.AddItem(1);
		m_lstOperations.AddItem(op_h);
		m_sName = "crx";
		break;
	case EBasicOperation::EBO_CRZ:
		m_lstQubits.AddItem(0);
		m_lstQubits.AddItem(1);
		m_lstOperations.AddItem(op_h);
		m_sName = "crx";
		break;

	case EBasicOperation::EBO_CC:
		m_lstQubits.AddItem(0);
		m_fClassicalParameter = m_fClassicalParameter + 0.00001; // to avoid 1.00000 been treated as 0.999999
		m_lstOperations.AddItem(op_h);
		m_sName = "collapse";
		break;
	default:
		break;
	}
}

QLGate::~QLGate()
{
	
}

QLGate::QLGate(const QLGate& other)
	: m_eOp(other.m_eOp)
	, m_bDagger(other.m_bDagger)
	, m_sName(other.m_sName)
	, m_lstQubits(other.m_lstQubits)
	, m_lstSubGates(other.m_lstSubGates)
	, m_fClassicalParameter(other.m_fClassicalParameter)
	, m_lstOperations(other.m_lstOperations)
{

}

const QLGate& QLGate::operator=(const QLGate& other)
{
	m_eOp = other.m_eOp;
	m_bDagger = other.m_bDagger;
	m_sName = other.m_sName;
	m_lstQubits = other.m_lstQubits;
	m_lstSubGates = other.m_lstSubGates;
	m_fClassicalParameter = other.m_fClassicalParameter;
	m_lstOperations = other.m_lstOperations;
	return *this;
}

void QLGate::Dagger()
{
	m_bDagger = !m_bDagger;
	if (EBasicOperation::EBO_Composite == m_eOp)
	{
		INT gateSize = m_lstSubGates.Num();
		for (INT i = 0; i < gateSize; ++i)
		{
			m_lstSubGates[i].Dagger();
		}
	}
}

QLGate QLGate::CreateControlled() const
{
	if (EBasicOperation::EBO_Composite != m_eOp)
	{
		TArray<BYTE> combinedQubits;
		combinedQubits.AddItem(0);
		for (INT i = 0; i < m_lstQubits.Num(); ++i)
		{
			combinedQubits.AddItem(m_lstQubits[i] + 1);
		}
		
		switch (m_eOp)
		{
		case EBasicOperation::EBO_H:
		{
			QLGate ret = CreateControlledZYZGate(_hadamard, FALSE);
			ret.ApplyOnQubits(combinedQubits);
			ret.m_sName = _T("c") + m_sName;
			ret.m_bDagger = m_bDagger;
			return ret;
		}
		case EBasicOperation::EBO_X:
		{
			QLGate ret(EBasicOperation::EBO_CX);
			ret.ApplyOnQubits(combinedQubits);
			ret.m_sName = _T("c") + m_sName;
			ret.m_bDagger = m_bDagger;
			return ret;
		}
		case EBasicOperation::EBO_Y:
		{
			QLGate ret(EBasicOperation::EBO_CY);
			ret.ApplyOnQubits(combinedQubits);
			ret.m_sName = _T("c") + m_sName;
			ret.m_bDagger = m_bDagger;
			return ret;
		}
		case EBasicOperation::EBO_Z:
		{
			QLGate ret(EBasicOperation::EBO_CZ);
			ret.ApplyOnQubits(combinedQubits);
			ret.m_sName = _T("c") + m_sName;
			ret.m_bDagger = m_bDagger;
			return ret;
		}
		case EBasicOperation::EBO_P:
		{
			QLGate ret(EBasicOperation::EBO_CP, m_fClassicalParameter);
			ret.ApplyOnQubits(combinedQubits);
			ret.m_sName = _T("c") + m_sName;
			ret.m_bDagger = m_bDagger;
			return ret;
		}
		case EBasicOperation::EBO_Phase:
		{
			TArray<BYTE> combinedQubits2;
			combinedQubits2.AddItem(0);
			QLGate ret(EBasicOperation::EBO_P, m_fClassicalParameter);
			ret.ApplyOnQubits(combinedQubits2);
			ret.m_sName = _T("c") + m_sName;
			ret.m_bDagger = m_bDagger;
			return ret;
		}

		case EBasicOperation::EBO_RX:
		{
			QLGate ret(EBasicOperation::EBO_CRX, m_fClassicalParameter);
			ret.ApplyOnQubits(combinedQubits);
			ret.m_sName = _T("c") + m_sName;
			ret.m_bDagger = m_bDagger;
			return ret;
		}
		case EBasicOperation::EBO_RY:
		{
			QLGate ret(EBasicOperation::EBO_CRY, m_fClassicalParameter);
			ret.ApplyOnQubits(combinedQubits);
			ret.m_sName = _T("c") + m_sName;
			ret.m_bDagger = m_bDagger;
			return ret;
		}
		case EBasicOperation::EBO_RZ:
		{
			QLGate ret(EBasicOperation::EBO_CRZ, m_fClassicalParameter);
			ret.ApplyOnQubits(combinedQubits);
			ret.m_sName = _T("c") + m_sName;
			ret.m_bDagger = m_bDagger;
			return ret;
		}

		case EBasicOperation::EBO_CX:
		{
			QLGate ret = CreateToffoliGate();
			ret.ApplyOnQubits(combinedQubits);
			ret.m_bDagger = m_bDagger;
			return ret;
		}
		case EBasicOperation::EBO_CY:
		{
			QLGate ret = CreateCnU(2, _PauliY);
			ret.m_sName = _T("CCY");
			ret.ApplyOnQubits(combinedQubits);
			ret.m_bDagger = m_bDagger;
			return ret;
		}
		case EBasicOperation::EBO_CZ:
		{
			QLGate ret = CreateCnU(2, _PauliZ);
			ret.m_sName = _T("CCZ");
			ret.ApplyOnQubits(combinedQubits);
			ret.m_bDagger = m_bDagger;
			return ret;
		}
		case EBasicOperation::EBO_CP:
		{
			QLGate ret = CreateCnP(2, m_fClassicalParameter);
			ret.ApplyOnQubits(combinedQubits);
			ret.m_bDagger = m_bDagger;
			return ret;
		}
		case EBasicOperation::EBO_CRX:
		{
			QLGate ret = CreateCnRX(2, m_fClassicalParameter);
			ret.ApplyOnQubits(combinedQubits);
			ret.m_bDagger = m_bDagger;
			return ret;
		}
		case EBasicOperation::EBO_CRY:
		{
			QLGate ret = CreateCnRY(2, m_fClassicalParameter);
			ret.ApplyOnQubits(combinedQubits);
			ret.m_bDagger = m_bDagger;
			return ret;
		}
		case EBasicOperation::EBO_CRZ:
		{
			QLGate ret = CreateCnRZ(2, m_fClassicalParameter);
			ret.ApplyOnQubits(combinedQubits);
			ret.m_bDagger = m_bDagger;
			return ret;
		}
		case EBasicOperation::EBO_CC:
			appCrucial("Not supported!");
			break;
		default:
			break;
		}

		appCrucial("something wrong!\n");
		return QLGate();
	}

	QLGate ret;
	ret.m_lstQubits.AddItem(0);
	for (INT i = 0; i < m_lstQubits.Num(); ++i)
	{
		ret.m_lstQubits.AddItem(m_lstQubits[i] + 1);
	}

	for (INT i = 0; i < m_lstSubGates.Num(); ++i)
	{
		ret.m_lstSubGates.AddItem(m_lstSubGates[i].CreateControlled());
	}
	ret.m_bDagger = m_bDagger;
	ret.m_sName = _T("c") + m_sName;
	return ret;
}

QLGate QLGate::Controlled(BYTE controlledQubit, const TArray<BYTE>& lstMappingQubits) const
{
	return QLGate();
}


void QLGate::AppendGate(QLGate toAppend, const TArray<BYTE>& lstMappingQubits)
{
	if (EBasicOperation::EBO_Composite != m_eOp)
	{
		printf("Basic operation gate not allowed to append gates\n");
		return;
	}

	toAppend.ApplyOnQubits(lstMappingQubits);
	m_lstSubGates.AddItem(toAppend);
}

TArray<SBasicOperation> QLGate::GetOperation(const TArray<BYTE>& lstMappingQubits) const
{
	TArray<SBasicOperation> ret;
	if (EBasicOperation::EBO_Composite != m_eOp)
	{
		INT opsize = m_lstOperations.Num();
		for (INT i = 0; i < opsize; ++i)
		{
			INT opIndex = m_bDagger ? (opsize - i - 1)  : i;
			SBasicOperation newone;
			newone.m_eOperation = m_lstOperations[opIndex].m_eOperation;
			INT qubitsize = m_lstQubits.Num();
			for (INT j = 0; j < qubitsize; ++j)
			{
				newone.m_lstQubits.AddItem(lstMappingQubits[m_lstQubits[j]]);
			}
			newone.m_fClassicalParameter = m_bDagger ? -m_fClassicalParameter : m_fClassicalParameter;
			if (EBasicOperation::EBO_CC == newone.m_eOperation && m_bDagger)
			{
				printf("dagger of controlled collapse not support!\n");
				newone.m_fClassicalParameter = m_fClassicalParameter;
			}

			ret.AddItem(newone);
		}
	}
	else
	{
		INT gatesize = m_lstSubGates.Num();
		for (INT i = 0; i < gatesize; ++i)
		{
			INT gateIndex = m_bDagger ? (gatesize - i - 1) : i;
			//TArray<BYTE> subgateQubits = ExchangeQubits(lstMappingQubits);
			TArray<SBasicOperation> oplist = m_lstSubGates[gateIndex].GetOperation(lstMappingQubits);
			ret.Append(oplist);
		}
	}
	return ret;
}

void QLGate::ApplyOnQubits(const TArray<BYTE>& lstMappingQubits)
{
	m_lstQubits = ExchangeQubits(lstMappingQubits); //thess changes are enough for basic gates
	if (EBasicOperation::EBO_Composite == m_eOp)
	{
		INT gatesize = m_lstSubGates.Num();
		for (INT i = 0; i < gatesize; ++i)
		{
			m_lstSubGates[i].ApplyOnQubits(lstMappingQubits);
		}
	}
}

TArray<BYTE> QLGate::ExchangeQubits(const TArray<BYTE>& lstMappingQubits) const
{
	TArray<BYTE> ret;
	INT qubitSize = m_lstQubits.Num();
	for (INT i = 0; i < qubitSize; ++i)
	{
		ret.AddItem(lstMappingQubits[m_lstQubits[i]]);
	}
	return ret;
}

void QLGate::PerformBasicOperation(const Qureg& reg, const SBasicOperation& op)
{
	switch (op.m_eOperation)
	{
	case EBasicOperation::EBO_H:
		hadamard(reg, op.m_lstQubits[0]);
		break;
	case EBasicOperation::EBO_X:
		pauliX(reg, op.m_lstQubits[0]);
		break;
	case EBasicOperation::EBO_Y:
		pauliY(reg, op.m_lstQubits[0]);
		break;
	case EBasicOperation::EBO_Z:
		pauliZ(reg, op.m_lstQubits[0]);
		break;
	case EBasicOperation::EBO_P:
		phaseShift(reg, op.m_lstQubits[0], op.m_fClassicalParameter);
		break;
	case EBasicOperation::EBO_RX:
		rotateX(reg, op.m_lstQubits[0], op.m_fClassicalParameter);
		break;
	case EBasicOperation::EBO_RY:
		rotateY(reg, op.m_lstQubits[0], op.m_fClassicalParameter);
		break;
	case EBasicOperation::EBO_RZ:
		rotateZ(reg, op.m_lstQubits[0], op.m_fClassicalParameter);
		break;

	case EBasicOperation::EBO_CX:
		controlledNot(reg, op.m_lstQubits[0], op.m_lstQubits[1]);
		break;
	case EBasicOperation::EBO_CP:
		controlledPhaseShift(reg, op.m_lstQubits[0], op.m_lstQubits[1], op.m_fClassicalParameter);
		break;
	case EBasicOperation::EBO_CRX:
		controlledRotateX(reg, op.m_lstQubits[0], op.m_lstQubits[1], op.m_fClassicalParameter);
		break;
	case EBasicOperation::EBO_CRY:
		controlledRotateY(reg, op.m_lstQubits[0], op.m_lstQubits[1], op.m_fClassicalParameter);
		break;
	case EBasicOperation::EBO_CRZ:
		controlledRotateZ(reg, op.m_lstQubits[0], op.m_lstQubits[1], op.m_fClassicalParameter);
		break;
	case EBasicOperation::EBO_CY:
		controlledPauliY(reg, op.m_lstQubits[0], op.m_lstQubits[1]);
		break;
	case EBasicOperation::EBO_CZ:
		controlledPhaseFlip(reg, op.m_lstQubits[0], op.m_lstQubits[1]);
		break;

	case EBasicOperation::EBO_CC:
		collapseToOutcome(reg, op.m_lstQubits[0], static_cast<INT>(op.m_fClassicalParameter));
		break;
	default:
		printf("not supported type!\n");
		break;
	}
}

void QLGate::DebugPrint(INT iDepth) const
{
	BYTE maxQ = 0;
	for (INT i = 0; i < m_lstQubits.Num(); ++i)
	{
		if (m_lstQubits[i] > maxQ)
		{
			maxQ = m_lstQubits[i];
		}
	}

	appPushLogDate(FALSE);

	PrintQubits(maxQ + 1);
	PrintOneGate(iDepth, maxQ + 1);
	PrintEmptyLine(maxQ + 1);

	appPopLogDate();
}

void QLGate::PrintQubits(BYTE qubitCount)
{
	PrintGateName(_T(""));
	for (INT i = 0; i < qubitCount; ++i)
	{
		if (i < 10)
		{
			appGeneral(_T("  %d  "), i);
		}
		else if (i < 100)
		{
			appGeneral(_T("  %d "), i);
		}
		else
		{
			appGeneral(_T(" %d "), i);
		}
	}
	appGeneral(_T("\n"));
}

void QLGate::PrintOneGate(INT iDepth, BYTE qubitCount) const
{
	if (EBasicOperation::EBO_Composite == m_eOp)
	{
		if (0 == iDepth)
		{
			PrintOneCompositeGate(qubitCount);
		}
		else
		{
			INT iNextDepth = iDepth > 0 ? (iDepth - 1) : iDepth;

			INT gatesize = m_lstSubGates.Num();
			UBOOL bHasComposite = FALSE;
			for (INT i = 0; i < gatesize; ++i)
			{
				if (EBasicOperation::EBO_Composite == m_lstSubGates[i].m_eOp)
				{
					bHasComposite = TRUE;
					break;
				}
			}

			if (!bHasComposite)
			{
				PrintGateName(_T(""));
				for (BYTE i = 0; i < qubitCount; ++i)
				{
					if (m_lstQubits.FindItemIndex(i) >= 0)
					{
						appGeneral(_T(" ─│─ "));
					}
					else
					{
						appGeneral(_T("  │  "));
					}
				}
				appGeneral(_T("\n"));
			}
			else
			{
				PrintEmptyLine(qubitCount);
			}

			for (INT i = 0; i < gatesize; ++i)
			{
				INT gateIndex = m_bDagger ? (gatesize - i - 1) : i;
				m_lstSubGates[gateIndex].PrintOneGate(iNextDepth, qubitCount);
			}

			if (!bHasComposite)
			{
				PrintGateName(_T(""));
				for (BYTE i = 0; i < qubitCount; ++i)
				{
					if (m_lstQubits.FindItemIndex(i) >= 0)
					{
						appGeneral(_T(" ─│─ "));
					}
					else
					{
						appGeneral(_T("  │  "));
					}
				}
				appGeneral(_T("\n"));
			}
		}
	}
	else
	{
		INT opsize = m_lstOperations.Num();
		for (INT i = 0; i < opsize; ++i)
		{
			INT opIndex = m_bDagger ? (opsize - i - 1) : i;
			SBasicOperation newone;
			newone.m_eOperation = m_lstOperations[opIndex].m_eOperation;
			INT qubitsize = m_lstQubits.Num();
			for (INT j = 0; j < qubitsize; ++j)
			{
				newone.m_lstQubits.AddItem(m_lstQubits[j]);
			}
			newone.m_fClassicalParameter = m_bDagger ? -m_fClassicalParameter : m_fClassicalParameter;
			if (EBasicOperation::EBO_CC == newone.m_eOperation && m_bDagger)
			{
				printf("dagger of controlled collapse not support!\n");
				newone.m_fClassicalParameter = m_fClassicalParameter;
			}

			PrintOneOp(newone, qubitCount);
		}
	}
}

void QLGate::PrintOneCompositeGate(BYTE qubitCount) const
{
	PrintEmptyLine(qubitCount);
	BYTE minQubit = qubitCount;
	BYTE maxQubit = 0;
	for (INT i = 0; i < m_lstQubits.Num(); ++i)
	{
		if (m_lstQubits[i] > maxQubit)
		{
			maxQubit = m_lstQubits[i];
		}

		if (m_lstQubits[i] < minQubit)
		{
			minQubit = m_lstQubits[i];
		}
	}

	PrintGateName(_T(""));
	for (BYTE qubit = 0; qubit < qubitCount; ++qubit)
	{
		if (qubit < minQubit)
		{
			appGeneral(_T("  │  "));
		}
		else if (qubit <= maxQubit)
		{
			if (qubit == minQubit && qubit == maxQubit)
			{
				appGeneral(_T("┌─┴─┐"));
			}
			else if (qubit == minQubit)
			{
				if (m_lstQubits.FindItemIndex(qubit) < 0)
				{
					appGeneral(_T("┌────"));
				}
				else
				{
					appGeneral(_T("┌─┴──"));
				}
			}
			else if (qubit == maxQubit)
			{
				if (m_lstQubits.FindItemIndex(qubit) < 0)
				{
					appGeneral(_T("────┐"));
				}
				else
				{
					appGeneral(_T("──┴─┐"));
				}
			}
			else
			{
				if (m_lstQubits.FindItemIndex(qubit) < 0)
				{
					appGeneral(_T("─────"));
				}
				else
				{
					appGeneral(_T("──┴──"));
				}
			}
		}
		else
		{
			appGeneral(_T("  │  "));
		}
	}
	appGeneral(_T("\n"));

	PrintGateName(m_bDagger ? (_T("d") + m_sName) : m_sName);
	for (BYTE qubit = 0; qubit < qubitCount; ++qubit)
	{
		if (qubit < minQubit)
		{
			appGeneral(_T("  │  "));
		}
		else if (qubit <= maxQubit)
		{
			if (qubit == minQubit && qubit == maxQubit)
			{
				appGeneral(_T("│ 0 │"));
			}
			else if (qubit == minQubit)
			{
				INT innerQubit = m_lstQubits.FindItemIndex(qubit);
				if (innerQubit >= 0)
				{
					if (innerQubit < 10)
					{
						appGeneral(_T("│ %d  "), innerQubit);
					}
					else if (innerQubit < 100)
					{
						appGeneral(_T("│ %d "), innerQubit);
					}
					else
					{
						appGeneral(_T("│%d"), innerQubit);
					}
				}
				else
				{
					appGeneral(_T("│    "));
				}
			}
			else if (qubit == maxQubit)
			{
				INT innerQubit = m_lstQubits.FindItemIndex(qubit);
				if (innerQubit >= 0)
				{
					if (innerQubit < 10)
					{
						appGeneral(_T("  %d │"), innerQubit);
					}
					else if (innerQubit < 100)
					{
						appGeneral(_T(" %d │"), innerQubit);
					}
					else
					{
						appGeneral(_T("%d│"), innerQubit);
					}
				}
				else
				{
					appGeneral(_T("    │"));
				}
			}
			else
			{
				INT innerQubit = m_lstQubits.FindItemIndex(qubit);
				if (innerQubit < 0)
				{
					appGeneral(_T("     "));
				}
				else
				{
					if (innerQubit < 10)
					{
						appGeneral(_T("  %d  "), innerQubit);
					}
					else if (innerQubit < 100)
					{
						appGeneral(_T("  %d "), innerQubit);
					}
					else
					{
						appGeneral(_T(" %d "), innerQubit);
					}
				}
			}
		}
		else
		{
			appGeneral(_T("  │  "));
		}
	}
	appGeneral(_T("\n"));

	PrintGateName(_T(""));
	for (BYTE qubit = 0; qubit < qubitCount; ++qubit)
	{
		if (qubit < minQubit)
		{
			appGeneral(_T("  │  "));
		}
		else if (qubit <= maxQubit)
		{
			if (qubit == minQubit && qubit == maxQubit)
			{
				appGeneral(_T("└─┬─┘"));
			}
			else if (qubit == minQubit)
			{
				if (m_lstQubits.FindItemIndex(qubit) < 0)
				{
					appGeneral(_T("└────"));
				}
				else
				{
					appGeneral(_T("└─┬──"));
				}
			}
			else if (qubit == maxQubit)
			{
				if (m_lstQubits.FindItemIndex(qubit) < 0)
				{
					appGeneral(_T("────┘"));
				}
				else
				{
					appGeneral(_T("──┬─┘"));
				}
			}
			else
			{
				if (m_lstQubits.FindItemIndex(qubit) < 0)
				{
					appGeneral(_T("─────"));
				}
				else
				{
					appGeneral(_T("──┬──"));
				}
			}
		}
		else
		{ 
			appGeneral(_T("  │  "));
		}
	}
	appGeneral(_T("\n"));
}

void QLGate::PrintOneOp(const SBasicOperation& op, BYTE qubitCount)
{
	CCString sName;
	switch (op.m_eOperation)
	{
	case EBasicOperation::EBO_H:
		PrintBasicSingle(op.m_lstQubits[0], qubitCount, sName, _T("H"));
		break;
	case EBasicOperation::EBO_X:
		PrintBasicSingle(op.m_lstQubits[0], qubitCount, sName, _T("X"));
		break;
	case EBasicOperation::EBO_Y:
		PrintBasicSingle(op.m_lstQubits[0], qubitCount, sName, _T("Y"));
		break;
	case EBasicOperation::EBO_Z:
		PrintBasicSingle(op.m_lstQubits[0], qubitCount, sName, _T("Z"));
		break;
	case EBasicOperation::EBO_P:
		sName.Format(_T("%s:%.3f"), _T("p"), op.m_fClassicalParameter);
		PrintBasicSingle(op.m_lstQubits[0], qubitCount, sName, _T("p"));
		break;
	case EBasicOperation::EBO_Phase:
		sName.Format(_T("%s:%.3f"), _T("ph"), op.m_fClassicalParameter);
		PrintBasicSingle(op.m_lstQubits[0], qubitCount, sName, _T("Ph"));
		break;
	case EBasicOperation::EBO_RX:
		sName.Format(_T("%s:%.3f"), _T("rx"), op.m_fClassicalParameter);
		PrintBasicSingle(op.m_lstQubits[0], qubitCount, sName, _T("Rx"));
		break;
	case EBasicOperation::EBO_RY:
		sName.Format(_T("%s:%.3f"), _T("ry"), op.m_fClassicalParameter);
		PrintBasicSingle(op.m_lstQubits[0], qubitCount, sName, _T("Ry"));
		break;
	case EBasicOperation::EBO_RZ:
		sName.Format(_T("%s:%.3f"), _T("rz"), op.m_fClassicalParameter);
		PrintBasicSingle(op.m_lstQubits[0], qubitCount, sName, _T("Rz"));
		break;
	case EBasicOperation::EBO_CX:
		PrintBasicTwo(op.m_lstQubits[0], op.m_lstQubits[1], qubitCount, sName, _T("Not"));
		break;
	case EBasicOperation::EBO_CY:
		PrintBasicTwo(op.m_lstQubits[0], op.m_lstQubits[1], qubitCount, sName, _T("Y"));
		break;
	case EBasicOperation::EBO_CZ:
		PrintBasicTwo(op.m_lstQubits[0], op.m_lstQubits[1], qubitCount, sName, _T("Z"));
		break;
	case EBasicOperation::EBO_CP:
		sName.Format(_T("%s:%.3f"), _T("cp"), op.m_fClassicalParameter);
		PrintBasicTwo(op.m_lstQubits[0], op.m_lstQubits[1], qubitCount, sName, _T("p"));
		break;
	case EBasicOperation::EBO_CRX:
		sName.Format(_T("%s:%.3f"), _T("crx"), op.m_fClassicalParameter);
		PrintBasicTwo(op.m_lstQubits[0], op.m_lstQubits[1], qubitCount, sName, _T("Rx"));
		break;
	case EBasicOperation::EBO_CRY:
		sName.Format(_T("%s:%.3f"), _T("cry"), op.m_fClassicalParameter);
		PrintBasicTwo(op.m_lstQubits[0], op.m_lstQubits[1], qubitCount, sName, _T("Ry"));
		break;
	case EBasicOperation::EBO_CRZ:
		sName.Format(_T("%s:%.3f"), _T("crz"), op.m_fClassicalParameter);
		PrintBasicTwo(op.m_lstQubits[0], op.m_lstQubits[1], qubitCount, sName, _T("Rz"));
		break;
	case EBasicOperation::EBO_CC:
		PrintBasicSingle(op.m_lstQubits[0], qubitCount, sName, _T("m"));
		break;
	default:
		appCrucial(_T("not supported!\n"));
		break;
	}
}

void QLGate::PrintEmptyLine(BYTE qubitCount)
{
	PrintGateName(_T(""));
	for (INT i = 0; i < qubitCount; ++i)
	{
		appGeneral(_T("  │  "));
	}
	appGeneral(_T("\n"));
}

void QLGate::PrintBasicSingle(BYTE target, BYTE qubitCount, const CCString& sName, const CCString& sInner)
{
	PrintGateName(_T(""));
	for (BYTE qubit = 0; qubit < qubitCount; ++qubit)
	{
		if (qubit == target)
		{
			appGeneral(_T("┌─┴─┐"));
		}
		else
		{
			appGeneral(_T("  │  "));
		}
	}
	appGeneral(_T("\n"));

	PrintGateName(sName);
	for (BYTE qubit = 0; qubit < qubitCount; ++qubit)
	{
		if (qubit == target)
		{
			if (1 == sInner.GetLength())
			{
				appGeneral(_T("│ %s │"), sInner.c_str());
			}
			else if (2 == sInner.GetLength())
			{
				appGeneral(_T("│ %s│"), sInner.c_str());
			}
			else if (3 == sInner.GetLength())
			{
				appGeneral(_T("│%s│"), sInner.c_str());
			}
			else
			{
				appGeneral(_T("│   │"));
			}
		}
		else
		{
			appGeneral(_T("  │  "));
		}
	}
	appGeneral(_T("\n"));

	PrintGateName(_T(""));
	for (BYTE qubit = 0; qubit < qubitCount; ++qubit)
	{
		if (qubit == target)
		{
			appGeneral(_T("└─┬─┘"));
		}
		else
		{
			appGeneral(_T("  │  "));
		}
	}
	appGeneral(_T("\n"));
}

void QLGate::PrintBasicTwo(BYTE ctr, BYTE target, BYTE qubitCount, const CCString& sName, const CCString& sInner)
{
	PrintGateName(_T(""));
	for (BYTE qubit = 0; qubit < qubitCount; ++qubit)
	{
		if (qubit == target)
		{
			appGeneral(_T("┌─┴─┐"));
		}
		else
		{
			appGeneral(_T("  │  "));
		}
	}
	appGeneral(_T("\n"));

	PrintGateName(sName);
	if (ctr < target)
	{
		for (BYTE qubit = 0; qubit < qubitCount; ++qubit)
		{
			if (qubit == ctr)
			{
				appGeneral(_T("  ├──"));
			}
			else if (qubit > ctr && qubit < target)
			{
				appGeneral(_T("─────"));
			}
			else if (qubit == target)
			{
				if (1 == sInner.GetLength())
				{
					appGeneral(_T("┤ %s │"), sInner.c_str());
				}
				else if (2 == sInner.GetLength())
				{
					appGeneral(_T("┤ %s│"), sInner.c_str());
				}
				else if (3 == sInner.GetLength())
				{
					appGeneral(_T("┤%s│"), sInner.c_str());
				}
				else
				{
					appGeneral(_T("┤   │"));
				}
			}
			else
			{
				appGeneral(_T("  │  "));
			}
		}
	}
	else
	{
		for (BYTE qubit = 0; qubit < qubitCount; ++qubit)
		{
			if (qubit == ctr)
			{
				appGeneral(_T("──┤  "));
			}
			else if (qubit < ctr && qubit > target)
			{
				appGeneral(_T("─────"));
			}
			else if (qubit == target)
			{
				if (1 == sInner.GetLength())
				{
					appGeneral(_T("│ %s ├"), sInner.c_str());
				}
				else if (2 == sInner.GetLength())
				{
					appGeneral(_T("│ %s├"), sInner.c_str());
				}
				else if (3 == sInner.GetLength())
				{
					appGeneral(_T("│%s├"), sInner.c_str());
				}
				else
				{
					appGeneral(_T("│   ├"));
				}
			}
			else
			{
				appGeneral(_T("  │  "));
			}
		}
	}
	appGeneral(_T("\n"));

	PrintGateName(_T(""));
	for (BYTE qubit = 0; qubit < qubitCount; ++qubit)
	{
		if (qubit == target)
		{
			appGeneral(_T("└─┬─┘"));
		}
		else
		{
			appGeneral(_T("  │  "));
		}
	}
	appGeneral(_T("\n"));
}

void QLGate::PrintGateName(const CCString& sName)
{
	CCString sNameString = _T("%s");
	for (INT i = 0; i < (_kPrintName - sName.GetLength()); ++i)
	{
		sNameString = sNameString + _T(" ");
	}
	appGeneral(sNameString.c_str(), sName.c_str());
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================