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
		op_h.m_lstQubits.AddItem(0);
		m_lstOperations.AddItem(op_h);
		m_sName = "H";
		break;
	case EBasicOperation::EBO_X:
		m_lstQubits.AddItem(0);
		op_h.m_lstQubits.AddItem(0);
		m_lstOperations.AddItem(op_h);
		m_sName = "X";
		break;
	case EBasicOperation::EBO_Y:
		m_lstQubits.AddItem(0);
		op_h.m_lstQubits.AddItem(0);
		m_lstOperations.AddItem(op_h);
		m_sName = "Y";
		break;
	case EBasicOperation::EBO_Z:
		m_lstQubits.AddItem(0);
		op_h.m_lstQubits.AddItem(0);
		m_lstOperations.AddItem(op_h);
		m_sName = "Z";
		break;
	case EBasicOperation::EBO_P:
		m_lstQubits.AddItem(0);
		op_h.m_lstQubits.AddItem(0);
		m_lstOperations.AddItem(op_h);
		m_sName = "p";
		break;
	case EBasicOperation::EBO_Phase:
	{
		m_lstQubits.AddItem(0);
		SBasicOperationInGate op_1;
		op_1.m_eOperation = EBasicOperation::EBO_P;
		op_1.m_lstQubits.AddItem(0);
		SBasicOperationInGate op_2;
		op_2.m_eOperation = EBasicOperation::EBO_X;
		op_2.m_lstQubits.AddItem(0);
		SBasicOperationInGate op_3;
		op_3.m_eOperation = EBasicOperation::EBO_P;
		op_3.m_lstQubits.AddItem(0);
		SBasicOperationInGate op_4;
		op_4.m_eOperation = EBasicOperation::EBO_X;
		op_4.m_lstQubits.AddItem(0);
		m_lstOperations.AddItem(op_1);
		m_lstOperations.AddItem(op_2);
		m_lstOperations.AddItem(op_3);
		m_lstOperations.AddItem(op_4);
		m_sName = "phase";
	}
		break;

	case EBasicOperation::EBO_RX:
		m_lstQubits.AddItem(0);
		op_h.m_lstQubits.AddItem(0);
		m_lstOperations.AddItem(op_h);
		m_sName = "rx";
		break;
	case EBasicOperation::EBO_RY:
		m_lstQubits.AddItem(0);
		op_h.m_lstQubits.AddItem(0);
		m_lstOperations.AddItem(op_h);
		m_sName = "ry";
		break;
	case EBasicOperation::EBO_RZ:
		m_lstQubits.AddItem(0);
		op_h.m_lstQubits.AddItem(0);
		m_lstOperations.AddItem(op_h);
		m_sName = "rz";
		break;

	case EBasicOperation::EBO_CX:
		m_lstQubits.AddItem(0);
		m_lstQubits.AddItem(1);
		op_h.m_lstQubits.AddItem(0);
		op_h.m_lstQubits.AddItem(1);
		m_lstOperations.AddItem(op_h);
		m_sName = "CNOT";
		break;
	case EBasicOperation::EBO_CY:
		m_lstQubits.AddItem(0);
		m_lstQubits.AddItem(1);
		op_h.m_lstQubits.AddItem(0);
		op_h.m_lstQubits.AddItem(1);
		m_lstOperations.AddItem(op_h);
		m_sName = "CY";
		break;
	case EBasicOperation::EBO_CZ:
		m_lstQubits.AddItem(0);
		m_lstQubits.AddItem(1);
		op_h.m_lstQubits.AddItem(0);
		op_h.m_lstQubits.AddItem(1);
		m_lstOperations.AddItem(op_h);
		m_sName = "CZ";
		break;

	case EBasicOperation::EBO_CP:
		m_lstQubits.AddItem(0);
		m_lstQubits.AddItem(1);
		op_h.m_lstQubits.AddItem(0);
		op_h.m_lstQubits.AddItem(1);
		m_lstOperations.AddItem(op_h);
		m_sName = "cp";
		break;
	case EBasicOperation::EBO_CRX:
		m_lstQubits.AddItem(0);
		m_lstQubits.AddItem(1);
		op_h.m_lstQubits.AddItem(0);
		op_h.m_lstQubits.AddItem(1);
		m_lstOperations.AddItem(op_h);
		m_sName = "crx";
		break;
	case EBasicOperation::EBO_CRY:
		m_lstQubits.AddItem(0);
		m_lstQubits.AddItem(1);
		op_h.m_lstQubits.AddItem(0);
		op_h.m_lstQubits.AddItem(1);
		m_lstOperations.AddItem(op_h);
		m_sName = "crx";
		break;
	case EBasicOperation::EBO_CRZ:
		m_lstQubits.AddItem(0);
		m_lstQubits.AddItem(1);
		op_h.m_lstQubits.AddItem(0);
		op_h.m_lstQubits.AddItem(1);
		m_lstOperations.AddItem(op_h);
		m_sName = "crx";
		break;

	case EBasicOperation::EBO_CC:
		m_lstQubits.AddItem(0);
		op_h.m_lstQubits.AddItem(0);
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
	return Controlled(0, m_lstQubits);
}

QLGate QLGate::Controlled(BYTE controlledQubit, const TArray<BYTE>& lstMappingQubits) const
{
	if (EBasicOperation::EBO_Composite != m_eOp)
	{
		TArray<BYTE> combinedQubits;
		combinedQubits.AddItem(controlledQubit);
		combinedQubits.Append(lstMappingQubits);
		switch (m_eOp)
		{
		case EBasicOperation::EBO_H:
		{
			QLGate ret = CreateControlledZYZGate(_hadamard, FALSE);
			ret.ApplyOnQubits(combinedQubits);
			return ret;
		}
		case EBasicOperation::EBO_X:
			{
				QLGate ret(EBasicOperation::EBO_CX);
				ret.ApplyOnQubits(combinedQubits);
				return ret;
			}
		case EBasicOperation::EBO_Y:
		{
			QLGate ret(EBasicOperation::EBO_CY);
			ret.ApplyOnQubits(combinedQubits);
			return ret;
		}
		case EBasicOperation::EBO_Z:
		{
			QLGate ret(EBasicOperation::EBO_CZ);
			ret.ApplyOnQubits(combinedQubits);
			return ret;
		}
		case EBasicOperation::EBO_P:
		{
			QLGate ret(EBasicOperation::EBO_CP, m_fClassicalParameter);
			ret.ApplyOnQubits(combinedQubits);
			return ret;
		}
		case EBasicOperation::EBO_Phase:
		{
			TArray<BYTE> combinedQubits2;
			combinedQubits2.AddItem(controlledQubit);
			QLGate ret(EBasicOperation::EBO_P, m_fClassicalParameter);
			ret.ApplyOnQubits(combinedQubits2);
			return ret;
		}

		case EBasicOperation::EBO_RX:
		{
			QLGate ret(EBasicOperation::EBO_CRX, m_fClassicalParameter);
			ret.ApplyOnQubits(combinedQubits);
			return ret;
		}
		case EBasicOperation::EBO_RY:
		{
			QLGate ret(EBasicOperation::EBO_CRY, m_fClassicalParameter);
			ret.ApplyOnQubits(combinedQubits);
			return ret;
		}
		case EBasicOperation::EBO_RZ:
		{
			QLGate ret(EBasicOperation::EBO_CRZ, m_fClassicalParameter);
			ret.ApplyOnQubits(combinedQubits);
			return ret;
		}

		case EBasicOperation::EBO_CX:
			//ccnot
			break;
		case EBasicOperation::EBO_CY:
			//ccy
			break;
		case EBasicOperation::EBO_CZ:
			//ccz
			break;

		case EBasicOperation::EBO_CP:
			//ccp
			break;
		case EBasicOperation::EBO_CRX:
			//ccrx
			break;
		case EBasicOperation::EBO_CRY:
			//ccry
			break;
		case EBasicOperation::EBO_CRZ:
			//ccrz
			break;

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
	ret.m_lstQubits.AddItem(controlledQubit);
	for (INT i = 0; i < m_lstQubits.Num(); ++i)
	{
		ret.m_lstQubits.AddItem(lstMappingQubits[m_lstQubits[i]]);
	}

	for (INT i = 0; i < m_lstSubGates.Num(); ++i)
	{
		ret.m_lstSubGates.AddItem(m_lstSubGates[i].Controlled(controlledQubit, lstMappingQubits));
	}
	return ret;
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
			INT qubitsize = m_lstOperations[opIndex].m_lstQubits.Num();
			for (INT j = 0; j < qubitsize; ++j)
			{
				newone.m_lstQubits.AddItem(lstMappingQubits[m_lstQubits[m_lstOperations[opIndex].m_lstQubits[j]]]);
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

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================