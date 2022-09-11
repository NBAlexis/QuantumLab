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
    : m_bBasicOperation(TRUE)
	, m_bDagger(FALSE)
	, m_fClassicalParameter(fParam)
{
	SBasicOperationInGate op_h;
	op_h.m_eOperation = eOp;

	switch (eOp)
	{
	case EBasicOperation::EBO_H:
		m_lstQubits.push_back(0);
		op_h.m_lstQubits.push_back(0);
		m_lstOperations.push_back(op_h);
		m_sName = "H";
		break;
	case EBasicOperation::EBO_X:
		m_lstQubits.push_back(0);
		op_h.m_lstQubits.push_back(0);
		m_lstOperations.push_back(op_h);
		m_sName = "X";
		break;
	case EBasicOperation::EBO_Y:
		m_lstQubits.push_back(0);
		op_h.m_lstQubits.push_back(0);
		m_lstOperations.push_back(op_h);
		m_sName = "Y";
		break;
	case EBasicOperation::EBO_Z:
		m_lstQubits.push_back(0);
		op_h.m_lstQubits.push_back(0);
		m_lstOperations.push_back(op_h);
		m_sName = "Z";
		break;
	case EBasicOperation::EBO_P:
		m_lstQubits.push_back(0);
		op_h.m_lstQubits.push_back(0);
		m_lstOperations.push_back(op_h);
		m_sName = "p";
		break;
	case EBasicOperation::EBO_Phase:
	{
		m_lstQubits.push_back(0);
		SBasicOperationInGate op_1;
		op_1.m_eOperation = EBasicOperation::EBO_P;
		op_1.m_lstQubits.push_back(0);
		SBasicOperationInGate op_2;
		op_2.m_eOperation = EBasicOperation::EBO_X;
		op_2.m_lstQubits.push_back(0);
		SBasicOperationInGate op_3;
		op_3.m_eOperation = EBasicOperation::EBO_P;
		op_3.m_lstQubits.push_back(0);
		SBasicOperationInGate op_4;
		op_4.m_eOperation = EBasicOperation::EBO_X;
		op_4.m_lstQubits.push_back(0);
		m_lstOperations.push_back(op_1);
		m_lstOperations.push_back(op_2);
		m_lstOperations.push_back(op_3);
		m_lstOperations.push_back(op_4);
		m_sName = "phase";
	}
		break;

	case EBasicOperation::EBO_RX:
		m_lstQubits.push_back(0);
		op_h.m_lstQubits.push_back(0);
		m_lstOperations.push_back(op_h);
		m_sName = "rx";
		break;
	case EBasicOperation::EBO_RY:
		m_lstQubits.push_back(0);
		op_h.m_lstQubits.push_back(0);
		m_lstOperations.push_back(op_h);
		m_sName = "ry";
		break;
	case EBasicOperation::EBO_RZ:
		m_lstQubits.push_back(0);
		op_h.m_lstQubits.push_back(0);
		m_lstOperations.push_back(op_h);
		m_sName = "rz";
		break;

	case EBasicOperation::EBO_CX:
		m_lstQubits.push_back(0);
		m_lstQubits.push_back(1);
		op_h.m_lstQubits.push_back(0);
		op_h.m_lstQubits.push_back(1);
		m_lstOperations.push_back(op_h);
		m_sName = "CNOT";
		break;
	case EBasicOperation::EBO_CY:
		m_lstQubits.push_back(0);
		m_lstQubits.push_back(1);
		op_h.m_lstQubits.push_back(0);
		op_h.m_lstQubits.push_back(1);
		m_lstOperations.push_back(op_h);
		m_sName = "CY";
		break;
	case EBasicOperation::EBO_CZ:
		m_lstQubits.push_back(0);
		m_lstQubits.push_back(1);
		op_h.m_lstQubits.push_back(0);
		op_h.m_lstQubits.push_back(1);
		m_lstOperations.push_back(op_h);
		m_sName = "CZ";
		break;

	case EBasicOperation::EBO_CP:
		m_lstQubits.push_back(0);
		m_lstQubits.push_back(1);
		op_h.m_lstQubits.push_back(0);
		op_h.m_lstQubits.push_back(1);
		m_lstOperations.push_back(op_h);
		m_sName = "cp";
		break;
	case EBasicOperation::EBO_CRX:
		m_lstQubits.push_back(0);
		m_lstQubits.push_back(1);
		op_h.m_lstQubits.push_back(0);
		op_h.m_lstQubits.push_back(1);
		m_lstOperations.push_back(op_h);
		m_sName = "crx";
		break;
	case EBasicOperation::EBO_CRY:
		m_lstQubits.push_back(0);
		m_lstQubits.push_back(1);
		op_h.m_lstQubits.push_back(0);
		op_h.m_lstQubits.push_back(1);
		m_lstOperations.push_back(op_h);
		m_sName = "crx";
		break;
	case EBasicOperation::EBO_CRZ:
		m_lstQubits.push_back(0);
		m_lstQubits.push_back(1);
		op_h.m_lstQubits.push_back(0);
		op_h.m_lstQubits.push_back(1);
		m_lstOperations.push_back(op_h);
		m_sName = "crx";
		break;

	case EBasicOperation::EBO_CC:
		m_lstQubits.push_back(0);
		op_h.m_lstQubits.push_back(0);
		m_fClassicalParameter = m_fClassicalParameter + 0.00001; // to avoid 1.00000 been treated as 0.999999
		m_lstOperations.push_back(op_h);
		m_sName = "collapse";
		break;
	default:
		break;
	}
}

QLGate::~QLGate()
{
	
}

void QLGate::Dagger()
{
	m_bDagger = !m_bDagger;
	if (!m_bBasicOperation)
	{
		SIZE_T gateSize = m_lstSubGates.size();
		for (SIZE_T i = 0; i < gateSize; ++i)
		{
			m_lstSubGates[i].Dagger();
		}
	}
}

QLGate QLGate::Controlled(BYTE controlledQubitCount) const
{
	return *this;
}


void QLGate::AppendGate(QLGate toAppend, const std::vector<BYTE>& lstMappingQubits)
{
	if (m_bBasicOperation)
	{
		printf("Basic operation gate not allowed to append gates\n");
		return;
	}

	toAppend.ApplyOnQubits(lstMappingQubits);
	m_lstSubGates.push_back(toAppend);
}

std::vector<SBasicOperation> QLGate::GetOperation(const std::vector<BYTE>& lstMappingQubits) const
{
	std::vector<SBasicOperation> ret;
	if (m_bBasicOperation)
	{
		SIZE_T opsize = m_lstOperations.size();
		for (SIZE_T i = 0; i < opsize; ++i)
		{
			SIZE_T opIndex = m_bDagger ? (opsize - i - 1)  : i;
			SBasicOperation newone;
			newone.m_eOperation = m_lstOperations[opIndex].m_eOperation;
			SIZE_T qubitsize = m_lstOperations[opIndex].m_lstQubits.size();
			for (SIZE_T j = 0; j < qubitsize; ++j)
			{
				newone.m_lstQubits.push_back(lstMappingQubits[m_lstOperations[opIndex].m_lstQubits[j]]);
			}
			newone.m_fClassicalParameter = m_bDagger ? -m_fClassicalParameter : m_fClassicalParameter;
			if (EBasicOperation::EBO_CC == newone.m_eOperation && m_bDagger)
			{
				printf("dagger of controlled collapse not support!\n");
				newone.m_fClassicalParameter = m_fClassicalParameter;
			}

			ret.push_back(newone);
		}
	}
	else
	{
		SIZE_T gatesize = m_lstSubGates.size();
		for (SIZE_T i = 0; i < gatesize; ++i)
		{
			SIZE_T gateIndex = m_bDagger ? (gatesize - i - 1) : i;
			std::vector<BYTE> subgateQubits = ExchangeQubits(lstMappingQubits);
			std::vector<SBasicOperation> oplist = m_lstSubGates[gateIndex].GetOperation(subgateQubits);
			ret.insert(ret.end(), oplist.begin(), oplist.end());
		}
	}
	return ret;
}

void QLGate::ApplyOnQubits(const std::vector<BYTE>& lstMappingQubits)
{
	m_lstQubits = ExchangeQubits(lstMappingQubits);
	if (!m_bBasicOperation)
	{
		SIZE_T gatesize = m_lstSubGates.size();
		for (SIZE_T i = 0; i < gatesize; ++i)
		{
			m_lstSubGates[i].ApplyOnQubits(m_lstQubits);
		}
	}
}

std::vector<BYTE> QLGate::ExchangeQubits(const std::vector<BYTE>& lstMappingQubits) const
{
	std::vector<BYTE> ret;
	SIZE_T qubitSize = m_lstQubits.size();
	for (SIZE_T i = 0; i < qubitSize; ++i)
	{
		ret.push_back(lstMappingQubits[m_lstQubits[i]]);
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