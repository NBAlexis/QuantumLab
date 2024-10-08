﻿//=============================================================================
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
		m_fClassicalParameter = m_fClassicalParameter + F(0.00001); // to avoid 1.00000 been treated as 0.999999
		m_lstOperations.AddItem(op_h);
		m_sName = "collapse";
		break;

	case EBasicOperation::EBO_Noise_Damping:
		m_lstQubits.AddItem(0);
		m_fClassicalParameter = m_fClassicalParameter;
		m_lstOperations.AddItem(op_h);
		m_sName = "damping";
		break;

	case EBasicOperation::EBO_Noise_Dephaseing:
		m_lstQubits.AddItem(0);
		m_fClassicalParameter = m_fClassicalParameter;
		m_lstOperations.AddItem(op_h);
		m_sName = "dephase";
		break;

	case EBasicOperation::EBO_Noise_Depolarising:
		m_lstQubits.AddItem(0);
		m_fClassicalParameter = m_fClassicalParameter;
		m_lstOperations.AddItem(op_h);
		m_sName = "depolar";
		break;

	case EBasicOperation::EBO_Noise_MixPauliX:
		m_lstQubits.AddItem(0);
		m_fClassicalParameter = m_fClassicalParameter;
		m_lstOperations.AddItem(op_h);
		m_sName = "mixPauliX";
		break;

	case EBasicOperation::EBO_Noise_MixPauliY:
		m_lstQubits.AddItem(0);
		m_fClassicalParameter = m_fClassicalParameter;
		m_lstOperations.AddItem(op_h);
		m_sName = "mixPauliY";
		break;

	case EBasicOperation::EBO_Noise_MixPauliZ:
		m_lstQubits.AddItem(0);
		m_fClassicalParameter = m_fClassicalParameter;
		m_lstOperations.AddItem(op_h);
		m_sName = "mixPauliZ";
		break;

	case EBasicOperation::EBO_Noise_MixPauliAll:
		m_lstQubits.AddItem(0);
		m_fClassicalParameter = m_fClassicalParameter;
		m_lstOperations.AddItem(op_h);
		m_sName = "mixPauli";
		break;

	case EBasicOperation::EBO_Noise_TwoQubitDephaseing:
		m_lstQubits.AddItem(0);
		m_lstQubits.AddItem(1);
		m_fClassicalParameter = m_fClassicalParameter;
		m_lstOperations.AddItem(op_h);
		m_sName = "dephase";
		break;

	case EBasicOperation::EBO_Noise_TwoQubitDepolarising:
		m_lstQubits.AddItem(0);
		m_lstQubits.AddItem(1);
		m_fClassicalParameter = m_fClassicalParameter;
		m_lstOperations.AddItem(op_h);
		m_sName = "depolar";
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
			QLGate ret = CreateControlledZYZGate(_hadamard);
			ret.ApplyOnQubits(combinedQubits);
			ret.m_sName = _T("c") + m_sName;
			if (m_bDagger)
			{
				ret.Dagger();
			}
			return ret;
		}
		case EBasicOperation::EBO_X:
		{
			QLGate ret(EBasicOperation::EBO_CX);
			ret.ApplyOnQubits(combinedQubits);
			ret.m_sName = _T("c") + m_sName;
			if (m_bDagger)
			{
				ret.Dagger();
			}
			return ret;
		}
		case EBasicOperation::EBO_Y:
		{
			QLGate ret(EBasicOperation::EBO_CY);
			ret.ApplyOnQubits(combinedQubits);
			ret.m_sName = _T("c") + m_sName;
			if (m_bDagger)
			{
				ret.Dagger();
			}
			return ret;
		}
		case EBasicOperation::EBO_Z:
		{
			QLGate ret(EBasicOperation::EBO_CZ);
			ret.ApplyOnQubits(combinedQubits);
			ret.m_sName = _T("c") + m_sName;
			if (m_bDagger)
			{
				ret.Dagger();
			}
			return ret;
		}
		case EBasicOperation::EBO_P:
		{
			QLGate ret(EBasicOperation::EBO_CP, m_fClassicalParameter);
			ret.ApplyOnQubits(combinedQubits);
			ret.m_sName = _T("c") + m_sName;
			if (m_bDagger)
			{
				ret.Dagger();
			}
			return ret;
		}
		case EBasicOperation::EBO_Phase:
		{
			TArray<BYTE> combinedQubits2;
			combinedQubits2.AddItem(0);
			QLGate ret(EBasicOperation::EBO_P, m_fClassicalParameter);
			//ret.ApplyOnQubits(combinedQubits2);
			ret.m_sName = _T("c") + m_sName;
			if (m_bDagger)
			{
				ret.Dagger();
			}
			return ret;
		}

		case EBasicOperation::EBO_RX:
		{
			QLGate ret(EBasicOperation::EBO_CRX, m_fClassicalParameter);
			ret.ApplyOnQubits(combinedQubits);
			ret.m_sName = _T("c") + m_sName;
			if (m_bDagger)
			{
				ret.Dagger();
			}
			return ret;
		}
		case EBasicOperation::EBO_RY:
		{
			QLGate ret(EBasicOperation::EBO_CRY, m_fClassicalParameter);
			ret.ApplyOnQubits(combinedQubits);
			ret.m_sName = _T("c") + m_sName;
			if (m_bDagger)
			{
				ret.Dagger();
			}
			return ret;
		}
		case EBasicOperation::EBO_RZ:
		{
			QLGate ret(EBasicOperation::EBO_CRZ, m_fClassicalParameter);
			ret.ApplyOnQubits(combinedQubits);
			ret.m_sName = _T("c") + m_sName;
			if (m_bDagger)
			{
				ret.Dagger();
			}
			return ret;
		}

		case EBasicOperation::EBO_CX:
		{
			QLGate ret = CreateToffoliGate();
			ret.ApplyOnQubits(combinedQubits);
			if (m_bDagger)
			{
				ret.Dagger();
			}
			return ret;
		}
		case EBasicOperation::EBO_CY:
		{
			QLGate ret = CreateCnU(2, _PauliY);
			ret.m_sName = _T("CCY");
			ret.ApplyOnQubits(combinedQubits);
			if (m_bDagger)
			{
				ret.Dagger();
			}
			return ret;
		}
		case EBasicOperation::EBO_CZ:
		{
			QLGate ret = CreateCnU(2, _PauliZ);
			ret.m_sName = _T("CCZ");
			ret.ApplyOnQubits(combinedQubits);
			if (m_bDagger)
			{
				ret.Dagger();
			}
			return ret;
		}
		case EBasicOperation::EBO_CP:
		{
			QLGate ret = CreateCnP(2, m_fClassicalParameter);
			ret.ApplyOnQubits(combinedQubits);
			if (m_bDagger)
			{
				ret.Dagger();
			}
			return ret;
		}
		case EBasicOperation::EBO_CRX:
		{
			QLGate ret = CreateCnRX(2, m_fClassicalParameter);
			ret.ApplyOnQubits(combinedQubits);
			if (m_bDagger)
			{
				ret.Dagger();
			}
			return ret;
		}
		case EBasicOperation::EBO_CRY:
		{
			QLGate ret = CreateCnRY(2, m_fClassicalParameter);
			ret.ApplyOnQubits(combinedQubits);
			if (m_bDagger)
			{
				ret.Dagger();
			}
			return ret;
		}
		case EBasicOperation::EBO_CRZ:
		{
			QLGate ret = CreateCnRZ(2, m_fClassicalParameter);
			ret.ApplyOnQubits(combinedQubits);
			if (m_bDagger)
			{
				ret.Dagger();
			}
			return ret;
		}
		case EBasicOperation::EBO_CC:
		case EBasicOperation::EBO_Noise_Damping:
		case EBasicOperation::EBO_Noise_Dephaseing:
		case EBasicOperation::EBO_Noise_Depolarising:
		case EBasicOperation::EBO_Noise_MixPauliX:
		case EBasicOperation::EBO_Noise_MixPauliY:
		case EBasicOperation::EBO_Noise_MixPauliZ:
		case EBasicOperation::EBO_Noise_MixPauliAll:
		case EBasicOperation::EBO_Noise_TwoQubitDephaseing:
		case EBasicOperation::EBO_Noise_TwoQubitDepolarising:
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
	QLGate cg = CreateControlled();
	TArray<BYTE> qubits;
	qubits.AddItem(controlledQubit);
	qubits.Append(lstMappingQubits);
	cg.ApplyOnQubits(qubits);
	return cg;
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

void QLGate::AppendGate(QLGate toAppend, BYTE qubit, ...)
{
	va_list arg;
	{
		va_start(arg, qubit);
		TArray<BYTE> qubits;
		qubits.AddItem(qubit);
		for (BYTE i = 1; i < toAppend.GetQubitCount(); ++i)
		{
			qubits.AddItem(va_arg(arg, BYTE));
		}

		AppendGate(toAppend, qubits);
		va_end(arg);
	}
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
			if (EBasicOperation::EBO_Noise_Damping == newone.m_eOperation && m_bDagger
			 || EBasicOperation::EBO_Noise_Dephaseing == newone.m_eOperation && m_bDagger
			 || EBasicOperation::EBO_Noise_Depolarising == newone.m_eOperation && m_bDagger
			 || EBasicOperation::EBO_Noise_MixPauliX == newone.m_eOperation && m_bDagger
			 || EBasicOperation::EBO_Noise_MixPauliY == newone.m_eOperation && m_bDagger
			 || EBasicOperation::EBO_Noise_MixPauliZ == newone.m_eOperation && m_bDagger
			 || EBasicOperation::EBO_Noise_MixPauliAll == newone.m_eOperation && m_bDagger
			 || EBasicOperation::EBO_Noise_TwoQubitDephaseing == newone.m_eOperation && m_bDagger
			 || EBasicOperation::EBO_Noise_TwoQubitDepolarising == newone.m_eOperation && m_bDagger)
			{
				appParanoiac(_T("dagger of noise simulation not support!\n"));
				newone.m_fClassicalParameter = m_fClassicalParameter;
			}
			if (EBasicOperation::EBO_CC == newone.m_eOperation && m_bDagger)
			{
				appParanoiac(_T("dagger of controlled collapse not support!\n"));
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

TArray<UINT> QLGate::SummarizeGateCounts() const
{
	UINT uiSingle = 0;
	UINT uiTwo = 0;
	UINT uiNoise = 0;
	UINT uiMeasure = 0;
	UINT uiOther = 0;
	if (EBasicOperation::EBO_Composite != m_eOp)
	{
		INT opsize = m_lstOperations.Num();
		for (INT i = 0; i < opsize; ++i)
		{

			if (SBasicOperation::IsSingle(m_lstOperations[i].m_eOperation))
			{
				++uiSingle;
			}
			else if (SBasicOperation::IsBasicTwoQubit(m_lstOperations[i].m_eOperation))
			{
				++uiTwo;
			}
			else if (SBasicOperation::IsNoise(m_lstOperations[i].m_eOperation))
			{
				++uiNoise;
			}
			else if (SBasicOperation::IsMeasure(m_lstOperations[i].m_eOperation))
			{
				++uiMeasure;
			}
			else
			{
				++uiOther;
			}
		}
	}
	else
	{
		INT gatesize = m_lstSubGates.Num();
		for (INT i = 0; i < gatesize; ++i)
		{
			TArray<UINT> oplist = m_lstSubGates[i].SummarizeGateCounts();
			uiSingle += oplist[0];
			uiTwo += oplist[1];
			uiNoise += oplist[2];
			uiMeasure += oplist[3];
			uiOther += oplist[4];
		}
	}

	TArray<UINT> ret;
	ret.AddItem(uiSingle);
	ret.AddItem(uiTwo);
	ret.AddItem(uiNoise);
	ret.AddItem(uiMeasure);
	ret.AddItem(uiOther);
	return ret;
}

void QLGate::ApplyOnQubits(const TArray<BYTE>& lstMappingQubits)
{
	m_lstQubits = ExchangeQubits(lstMappingQubits); //thess changes are enough for basic gates
	TArray<BYTE> additionalQubitData;
	INT qubitSize = m_lstAdditionalQubitsData.Num();
	for (INT i = 0; i < qubitSize; ++i)
	{
		additionalQubitData.AddItem(lstMappingQubits[m_lstAdditionalQubitsData[i]]);
	}
	m_lstAdditionalQubitsData = additionalQubitData;
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

Real QLGate::PerformBasicOperation(const Qureg& reg, const SBasicOperation& op)
{
	Real ret = F(1.0);
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
		ret = static_cast<Real>(collapseToOutcome(reg, op.m_lstQubits[0], static_cast<INT>(op.m_fClassicalParameter)));
		break;

	case EBasicOperation::EBO_Noise_Damping:
		if (reg.isDensityMatrix)
		{
			mixDamping(reg, op.m_lstQubits[0], op.m_fClassicalParameter);
		}
		else
		{
			appDetailed(_T("Noise gate encountered, but the state is pure, ignored...\n"));
		}
		break;
	case EBasicOperation::EBO_Noise_Dephaseing:
		if (reg.isDensityMatrix)
		{
			mixDephasing(reg, op.m_lstQubits[0], op.m_fClassicalParameter);
		}
		else
		{
			appDetailed(_T("Noise gate encountered, but the state is pure, ignored...\n"));
		}
		break;
	case EBasicOperation::EBO_Noise_Depolarising:
		if (reg.isDensityMatrix)
		{
			mixDepolarising(reg, op.m_lstQubits[0], op.m_fClassicalParameter);
		}
		else
		{
			appDetailed(_T("Noise gate encountered, but the state is pure, ignored...\n"));
		}
		break;
	case EBasicOperation::EBO_Noise_MixPauliX:
		if (reg.isDensityMatrix)
		{
			mixPauli(reg, op.m_lstQubits[0], op.m_fClassicalParameter, F(0.0), F(0.0));
		}
		else
		{
			appDetailed(_T("Noise gate encountered, but the state is pure, ignored...\n"));
		}
		break;
	case EBasicOperation::EBO_Noise_MixPauliY:
		if (reg.isDensityMatrix)
		{
			mixPauli(reg, op.m_lstQubits[0], F(0.0), op.m_fClassicalParameter, F(0.0));
		}
		else
		{
			appDetailed(_T("Noise gate encountered, but the state is pure, ignored...\n"));
		}
		break;
	case EBasicOperation::EBO_Noise_MixPauliZ:
		if (reg.isDensityMatrix)
		{
			mixPauli(reg, op.m_lstQubits[0], F(0.0), F(0.0), op.m_fClassicalParameter);
		}
		else
		{
			appDetailed(_T("Noise gate encountered, but the state is pure, ignored...\n"));
		}
		break;
	case EBasicOperation::EBO_Noise_MixPauliAll:
		if (reg.isDensityMatrix)
		{
			mixPauli(reg, op.m_lstQubits[0], op.m_fClassicalParameter, op.m_fClassicalParameter, op.m_fClassicalParameter);
		}
		else
		{
			appDetailed(_T("Noise gate encountered, but the state is pure, ignored...\n"));
		}
		break;
	case EBasicOperation::EBO_Noise_TwoQubitDephaseing:
		if (reg.isDensityMatrix)
		{
			mixTwoQubitDephasing(reg, op.m_lstQubits[0], op.m_lstQubits[1], op.m_fClassicalParameter);
		}
		else
		{
			appDetailed(_T("Noise gate encountered, but the state is pure, ignored...\n"));
		}
		break;
	case EBasicOperation::EBO_Noise_TwoQubitDepolarising:
		if (reg.isDensityMatrix)
		{
			mixTwoQubitDepolarising(reg, op.m_lstQubits[0], op.m_lstQubits[1], op.m_fClassicalParameter);
		}
		else
		{
			appDetailed(_T("Noise gate encountered, but the state is pure, ignored...\n"));
		}
		break;

	default:
		printf("not supported type!\n");
		break;
	}
	return ret;
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
	case EBasicOperation::EBO_Noise_Damping:
	case EBasicOperation::EBO_Noise_Dephaseing:
	case EBasicOperation::EBO_Noise_Depolarising:
	case EBasicOperation::EBO_Noise_MixPauliX:
	case EBasicOperation::EBO_Noise_MixPauliY:
	case EBasicOperation::EBO_Noise_MixPauliZ:
	case EBasicOperation::EBO_Noise_MixPauliAll:
		PrintBasicSingle(op.m_lstQubits[0], qubitCount, sName, _T("N"));
		break;
	case EBasicOperation::EBO_Noise_TwoQubitDephaseing:
	case EBasicOperation::EBO_Noise_TwoQubitDepolarising:
		PrintBasicTwo(op.m_lstQubits[0], op.m_lstQubits[1], qubitCount, sName, _T("N"));
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


TArray<SBasicOperation> QLGate::ToBuiltIn(const TArray<SBasicOperation>& gates)
{
	TArray<SBasicOperation> ret;

	QLGate p1(EBasicOperation::EBO_P, -PI / 2);
	QLGate p2(EBasicOperation::EBO_P, PI / 2);
	QLGate h(EBasicOperation::EBO_H);
	QLGate cx(EBasicOperation::EBO_CX);

	for (INT i = 0; i < gates.Num(); ++i)
	{
		switch (gates[i].m_eOperation)
		{
		case EBasicOperation::EBO_CC:
		case EBasicOperation::EBO_H:
		case EBasicOperation::EBO_X:
		case EBasicOperation::EBO_Y:
		case EBasicOperation::EBO_Z:
		case EBasicOperation::EBO_P:
		case EBasicOperation::EBO_RX:
		case EBasicOperation::EBO_RY:
		case EBasicOperation::EBO_RZ:
		case EBasicOperation::EBO_CX:
			{
				ret.AddItem(gates[i]);
			}
			break;
		case EBasicOperation::EBO_CY:
			{
				SBasicOperation cy_h1;
				cy_h1.m_eOperation = EBasicOperation::EBO_H;
				cy_h1.m_lstQubits.AddItem(gates[i].m_lstQubits[1]);
				ret.AddItem(cy_h1);

				SBasicOperation cy_p1;
				cy_p1.m_eOperation = EBasicOperation::EBO_P;
				cy_p1.m_lstQubits.AddItem(gates[i].m_lstQubits[1]);
				cy_p1.m_fClassicalParameter = PI / 2;
				ret.AddItem(cy_p1);

				SBasicOperation cy_cx;
				cy_cx.m_eOperation = EBasicOperation::EBO_CX;
				cy_cx.m_lstQubits.AddItem(gates[i].m_lstQubits[0]);
				cy_cx.m_lstQubits.AddItem(gates[i].m_lstQubits[1]);
				ret.AddItem(cy_cx);

				SBasicOperation cy_p2;
				cy_p2.m_eOperation = EBasicOperation::EBO_P;
				cy_p2.m_lstQubits.AddItem(gates[i].m_lstQubits[1]);
				cy_p2.m_fClassicalParameter = -PI / 2;
				ret.AddItem(cy_p2);

				SBasicOperation cy_h2;
				cy_h2.m_eOperation = EBasicOperation::EBO_H;
				cy_h2.m_lstQubits.AddItem(gates[i].m_lstQubits[1]);
				ret.AddItem(cy_h2);
			}
			break;
		case EBasicOperation::EBO_CZ:
			{
				SBasicOperation cz_h1;
				cz_h1.m_eOperation = EBasicOperation::EBO_H;
				cz_h1.m_lstQubits.AddItem(gates[i].m_lstQubits[1]);
				ret.AddItem(cz_h1);

				SBasicOperation cz_cx;
				cz_cx.m_eOperation = EBasicOperation::EBO_CX;
				cz_cx.m_lstQubits.AddItem(gates[i].m_lstQubits[0]);
				cz_cx.m_lstQubits.AddItem(gates[i].m_lstQubits[1]);
				ret.AddItem(cz_cx);

				SBasicOperation cz_h2;
				cz_h2.m_eOperation = EBasicOperation::EBO_H;
				cz_h2.m_lstQubits.AddItem(gates[i].m_lstQubits[1]);
				ret.AddItem(cz_h2);
			}
			break;
		case EBasicOperation::EBO_CP:
			{
				SBasicOperation cp_cx1;
				cp_cx1.m_eOperation = EBasicOperation::EBO_CX;
				cp_cx1.m_lstQubits.AddItem(gates[i].m_lstQubits[0]);
				cp_cx1.m_lstQubits.AddItem(gates[i].m_lstQubits[1]);
				ret.AddItem(cp_cx1);

				SBasicOperation cp_p1;
				cp_p1.m_eOperation = EBasicOperation::EBO_P;
				cp_p1.m_lstQubits.AddItem(gates[i].m_lstQubits[1]);
				cp_p1.m_fClassicalParameter = -gates[i].m_fClassicalParameter / 2;
				ret.AddItem(cp_p1);

				SBasicOperation cp_cx2;
				cp_cx2.m_eOperation = EBasicOperation::EBO_CX;
				cp_cx2.m_lstQubits.AddItem(gates[i].m_lstQubits[0]);
				cp_cx2.m_lstQubits.AddItem(gates[i].m_lstQubits[1]);
				ret.AddItem(cp_cx2);

				SBasicOperation cp_p2;
				cp_p2.m_eOperation = EBasicOperation::EBO_P;
				cp_p2.m_lstQubits.AddItem(gates[i].m_lstQubits[1]);
				cp_p2.m_fClassicalParameter = gates[i].m_fClassicalParameter / 2;
				ret.AddItem(cp_p2);

				SBasicOperation cp_p3;
				cp_p3.m_eOperation = EBasicOperation::EBO_P;
				cp_p3.m_lstQubits.AddItem(gates[i].m_lstQubits[0]);
				cp_p3.m_fClassicalParameter = gates[i].m_fClassicalParameter / 2;
				ret.AddItem(cp_p3);
			}
			break;
		case EBasicOperation::EBO_CRX:
			{
				SBasicOperation crx_rz1;
				crx_rz1.m_eOperation = EBasicOperation::EBO_RZ;
				crx_rz1.m_lstQubits.AddItem(gates[i].m_lstQubits[1]);
				crx_rz1.m_fClassicalParameter = -PI / 2;
				ret.AddItem(crx_rz1);

				SBasicOperation crx_cx1;
				crx_cx1.m_eOperation = EBasicOperation::EBO_CX;
				crx_cx1.m_lstQubits.AddItem(gates[i].m_lstQubits[0]);
				crx_cx1.m_lstQubits.AddItem(gates[i].m_lstQubits[1]);
				ret.AddItem(crx_cx1);

				SBasicOperation crx_ry1;
				crx_ry1.m_eOperation = EBasicOperation::EBO_RY;
				crx_ry1.m_lstQubits.AddItem(gates[i].m_lstQubits[1]);
				crx_ry1.m_fClassicalParameter = gates[i].m_fClassicalParameter / 2;
				ret.AddItem(crx_ry1);

				SBasicOperation crx_cx2;
				crx_cx2.m_eOperation = EBasicOperation::EBO_CX;
				crx_cx2.m_lstQubits.AddItem(gates[i].m_lstQubits[0]);
				crx_cx2.m_lstQubits.AddItem(gates[i].m_lstQubits[1]);
				ret.AddItem(crx_cx2);

				SBasicOperation crx_ry2;
				crx_ry2.m_eOperation = EBasicOperation::EBO_RY;
				crx_ry2.m_lstQubits.AddItem(gates[i].m_lstQubits[1]);
				crx_ry2.m_fClassicalParameter = -gates[i].m_fClassicalParameter / 2;
				ret.AddItem(crx_ry2);

				SBasicOperation crx_rz2;
				crx_rz2.m_eOperation = EBasicOperation::EBO_RZ;
				crx_rz2.m_lstQubits.AddItem(gates[i].m_lstQubits[1]);
				crx_rz2.m_fClassicalParameter = PI / 2;
				ret.AddItem(crx_rz2);
			}
			break;
		case EBasicOperation::EBO_CRY:
			{
				SBasicOperation cry_cx1;
				cry_cx1.m_eOperation = EBasicOperation::EBO_CX;
				cry_cx1.m_lstQubits.AddItem(gates[i].m_lstQubits[0]);
				cry_cx1.m_lstQubits.AddItem(gates[i].m_lstQubits[1]);
				ret.AddItem(cry_cx1);

				SBasicOperation cry_ry1;
				cry_ry1.m_eOperation = EBasicOperation::EBO_RY;
				cry_ry1.m_lstQubits.AddItem(gates[i].m_lstQubits[1]);
				cry_ry1.m_fClassicalParameter = -gates[i].m_fClassicalParameter / 2;
				ret.AddItem(cry_ry1);

				SBasicOperation cry_cx2;
				cry_cx2.m_eOperation = EBasicOperation::EBO_CX;
				cry_cx2.m_lstQubits.AddItem(gates[i].m_lstQubits[0]);
				cry_cx2.m_lstQubits.AddItem(gates[i].m_lstQubits[1]);
				ret.AddItem(cry_cx2);

				SBasicOperation cry_ry2;
				cry_ry2.m_eOperation = EBasicOperation::EBO_RY;
				cry_ry2.m_lstQubits.AddItem(gates[i].m_lstQubits[1]);
				cry_ry2.m_fClassicalParameter = gates[i].m_fClassicalParameter / 2;
				ret.AddItem(cry_ry2);
			}
			break;
		case EBasicOperation::EBO_CRZ:
			{
				SBasicOperation crz_cx1;
				crz_cx1.m_eOperation = EBasicOperation::EBO_CX;
				crz_cx1.m_lstQubits.AddItem(gates[i].m_lstQubits[0]);
				crz_cx1.m_lstQubits.AddItem(gates[i].m_lstQubits[1]);
				ret.AddItem(crz_cx1);

				SBasicOperation crz_rz1;
				crz_rz1.m_eOperation = EBasicOperation::EBO_RZ;
				crz_rz1.m_lstQubits.AddItem(gates[i].m_lstQubits[1]);
				crz_rz1.m_fClassicalParameter = -gates[i].m_fClassicalParameter / 2;
				ret.AddItem(crz_rz1);

				SBasicOperation crz_cx2;
				crz_cx2.m_eOperation = EBasicOperation::EBO_CX;
				crz_cx2.m_lstQubits.AddItem(gates[i].m_lstQubits[0]);
				crz_cx2.m_lstQubits.AddItem(gates[i].m_lstQubits[1]);
				ret.AddItem(crz_cx2);

				SBasicOperation crz_rz2;
				crz_rz2.m_eOperation = EBasicOperation::EBO_RZ;
				crz_rz2.m_lstQubits.AddItem(gates[i].m_lstQubits[1]);
				crz_rz2.m_fClassicalParameter = gates[i].m_fClassicalParameter / 2;
				ret.AddItem(crz_rz2);
			}
			break;
		default:
			appParanoiac(_T("not supported operator : %d\n"), gates[i].m_eOperation);
			break;
		}
	}
	return ret;
}

/**
* About phase:
* in normal openQASM,
* U(a,b,c)=Rz(b)Ry(a)Rz(c)|psi>
*         =  [ exp(i (alpha+beta)/2) cos(theta/2)  exp(i (alpha-beta)/2) sin(theta/2)]
*            [-exp(-i(alpha-beta)/2) sin(theta/2)  exp(-i(alpha+beta)/2) cos(theta/2)]
* 
* in qiskit
* U(a,b,c) = [cos(a/2)           -exp(ic)    sin(a/2)]
*            [exp(ib) sin(a/2)   exp(i(b+c)) cos(a/2)]
* 
* 
*/
CCString QLGate::ToOpenOASM(const TArray<SBasicOperation>& gates, BYTE byQubit, const TArray<BYTE>& measureAtLast)
{
	CCString ret = _T("");

	ret = ret + _T("qreg q[") + appToString(byQubit) + _T("];\n");
	BYTE byMeasurebit = 0;
	THashMap<BYTE, BYTE> measuremap;

	QLMatrix single = _I2;
	CCString operations = _T("");
	CCString sOpAdd;
	TArray<QLMatrix> singlematrices;
	TArray<UINT> qubitLevel;
	Real gloablePhaseOpenQASM = F(0.0);
	Real gloablePhaseOpenQiskit = F(0.0);

	for (BYTE i = 0; i < byQubit; ++i)
	{
		singlematrices.AddItem(single);
		qubitLevel.AddItem(0);
	}

	for (INT i = 0; i < gates.Num(); ++i)
	{
		//UBOOL bIsTwoQubits = FALSE;
		//BYTE dirtyQubit1 = 0;
		//BYTE dirtyQubit2 = 0;

		switch (gates[i].m_eOperation)
		{
		case EBasicOperation::EBO_CC:
			{
				TArray<Real> degrees = GetZYZDecompose(singlematrices[gates[i].m_lstQubits[0]]);
				if (abs(degrees[0]) + abs(degrees[1]) + abs(degrees[2]) > F(0.00001))
				{
					sOpAdd.Format(_T("U(%.8f, %.8f, %.8f) q[%d];\n"), -degrees[0], -degrees[2], -degrees[1], gates[i].m_lstQubits[0]);
					gloablePhaseOpenQASM += degrees[3];
					gloablePhaseOpenQiskit += degrees[3] + F(0.5) * (degrees[1] + degrees[2]);
					operations = operations + sOpAdd;
				}
				else
				{
					gloablePhaseOpenQASM += degrees[3];
					gloablePhaseOpenQiskit += degrees[3] + F(0.5) * (degrees[1] + degrees[2]);
				}

				singlematrices[gates[i].m_lstQubits[0]] = _I2;
				qubitLevel[gates[i].m_lstQubits[0]] = qubitLevel[gates[i].m_lstQubits[0]] + 1;

				//find a classical bit
				BYTE classicalByte = 0;
				if (measuremap.Exist(gates[i].m_lstQubits[0]))
				{
					classicalByte = measuremap[gates[i].m_lstQubits[0]];
				}
				else
				{
					classicalByte = byMeasurebit;
					measuremap.SetAt(gates[i].m_lstQubits[0], classicalByte);
					byMeasurebit = byMeasurebit + 1;
				}

				//apply measure
				sOpAdd.Format(_T("measure q[%d] -> c[%d];\n"), gates[i].m_lstQubits[0], classicalByte);
				operations = operations + sOpAdd;
			}
			break;
		case EBasicOperation::EBO_H:
		case EBasicOperation::EBO_X:
		case EBasicOperation::EBO_Y:
		case EBasicOperation::EBO_Z:
		case EBasicOperation::EBO_P:
		case EBasicOperation::EBO_RX:
		case EBasicOperation::EBO_RY:
		case EBasicOperation::EBO_RZ:
			{
				singlematrices[gates[i].m_lstQubits[0]] = CreateSingleQubitMatrix(gates[i].m_eOperation, gates[i].m_fClassicalParameter) * singlematrices[gates[i].m_lstQubits[0]];
			}
			break;
		case EBasicOperation::EBO_CX:
		case EBasicOperation::EBO_CY:
		case EBasicOperation::EBO_CZ:
		case EBasicOperation::EBO_CP:
		case EBasicOperation::EBO_CRX:
		case EBasicOperation::EBO_CRY:
		case EBasicOperation::EBO_CRZ:
			{
				TArray<Real> degrees = GetZYZDecompose(singlematrices[gates[i].m_lstQubits[0]]);
				if (abs(degrees[0]) + abs(degrees[1]) + abs(degrees[2]) > F(0.00001))
				{
					sOpAdd.Format(_T("U(%.8f, %.8f, %.8f) q[%d];\n"), -degrees[0], -degrees[2], -degrees[1], gates[i].m_lstQubits[0]);
					gloablePhaseOpenQASM += degrees[3];
					gloablePhaseOpenQiskit += degrees[3] + F(0.5) * (degrees[1] + degrees[2]);
					operations = operations + sOpAdd;
				}
				else
				{
					gloablePhaseOpenQASM += degrees[3];
					gloablePhaseOpenQiskit += degrees[3] + F(0.5) * (degrees[1] + degrees[2]);
				}
				degrees = GetZYZDecompose(singlematrices[gates[i].m_lstQubits[1]]);
				if (abs(degrees[0]) + abs(degrees[1]) + abs(degrees[2]) > F(0.00001))
				{
					sOpAdd.Format(_T("U(%.8f, %.8f, %.8f) q[%d];\n"), -degrees[0], -degrees[2], -degrees[1], gates[i].m_lstQubits[1]);
					gloablePhaseOpenQASM += degrees[3];
					gloablePhaseOpenQiskit += degrees[3] + F(0.5) * (degrees[1] + degrees[2]);
					operations = operations + sOpAdd;
				}
				else
				{
					gloablePhaseOpenQASM += degrees[3];
					gloablePhaseOpenQiskit += degrees[3] + F(0.5) * (degrees[1] + degrees[2]);
				}
				singlematrices[gates[i].m_lstQubits[0]] = _I2;
				singlematrices[gates[i].m_lstQubits[1]] = _I2;
				UINT iLevel1 = qubitLevel[gates[i].m_lstQubits[0]] + 1;
				UINT iLevel2 = qubitLevel[gates[i].m_lstQubits[1]] + 1;
				qubitLevel[gates[i].m_lstQubits[0]] = iLevel1 > iLevel2 ? iLevel1 : iLevel2;
				qubitLevel[gates[i].m_lstQubits[1]] = iLevel1 > iLevel2 ? iLevel1 : iLevel2;

				switch (gates[i].m_eOperation)
					{
						case EBasicOperation::EBO_CX:
						{
							sOpAdd.Format(_T("CX q[%d], q[%d];\n"), gates[i].m_lstQubits[0], gates[i].m_lstQubits[1]);
							operations = operations + sOpAdd;
						}
						break;
						case EBasicOperation::EBO_CY:
						{
							sOpAdd.Format(_T("cy q[%d], q[%d];\n"), gates[i].m_lstQubits[0], gates[i].m_lstQubits[1]);
							operations = operations + sOpAdd;
						}
						break;
						case EBasicOperation::EBO_CZ:
						{
							sOpAdd.Format(_T("cz q[%d], q[%d];\n"), gates[i].m_lstQubits[0], gates[i].m_lstQubits[1]);
							operations = operations + sOpAdd;
						}
						break;
						case EBasicOperation::EBO_CP:
						{
							sOpAdd.Format(_T("cp(%.8f) q[%d], q[%d];\n"), gates[i].m_fClassicalParameter, gates[i].m_lstQubits[0], gates[i].m_lstQubits[1]);
							operations = operations + sOpAdd;
						}
						break;
						case EBasicOperation::EBO_CRX:
						{
							sOpAdd.Format(_T("crx(%.8f) q[%d], q[%d];\n"), gates[i].m_fClassicalParameter, gates[i].m_lstQubits[0], gates[i].m_lstQubits[1]);
							operations = operations + sOpAdd;
						}
						break;
						case EBasicOperation::EBO_CRY:
						{
							sOpAdd.Format(_T("cry(%.8f) q[%d], q[%d];\n"), gates[i].m_fClassicalParameter, gates[i].m_lstQubits[0], gates[i].m_lstQubits[1]);
							operations = operations + sOpAdd;
						}
						break;
						case EBasicOperation::EBO_CRZ:
						{
							sOpAdd.Format(_T("crz(%.8f) q[%d], q[%d];\n"), gates[i].m_fClassicalParameter, gates[i].m_lstQubits[0], gates[i].m_lstQubits[1]);
							operations = operations + sOpAdd;
						}
						break;
					}

				}
			break;
		}
	}

	UINT uiMaxLevel = 0;
	for (BYTE i = 0; i < byQubit; ++i)
	{
		TArray<Real> degrees = GetZYZDecompose(singlematrices[i]);
		if (abs(degrees[0]) + abs(degrees[1]) + abs(degrees[2]) > F(0.00001))
		{
			sOpAdd.Format(_T("U(%.8f, %.8f, %.8f) q[%d];\n"), -degrees[0], -degrees[2], -degrees[1], i);
			gloablePhaseOpenQASM += degrees[3];
			gloablePhaseOpenQiskit += degrees[3] + F(0.5) * (degrees[1] + degrees[2]);
			operations = operations + sOpAdd;
		}
		else
		{
			gloablePhaseOpenQASM += degrees[3];
			gloablePhaseOpenQiskit += degrees[3] + F(0.5) * (degrees[1] + degrees[2]);
		}

		if (uiMaxLevel < qubitLevel[i])
		{
			uiMaxLevel = qubitLevel[i];
		}
	}

	for (INT i = 0; i < measureAtLast.Num(); ++i)
	{
		BYTE classicalByte = 0;
		if (measuremap.Exist(gates[i].m_lstQubits[0]))
		{
			classicalByte = measuremap[gates[i].m_lstQubits[0]];
		}
		else
		{
			classicalByte = byMeasurebit;
			measuremap.SetAt(gates[i].m_lstQubits[0], classicalByte);
			byMeasurebit = byMeasurebit + 1;
		}

		//apply measure
		sOpAdd.Format(_T("measure q[%d] -> c[%d];\n"), gates[i].m_lstQubits[0], classicalByte);
		operations = operations + sOpAdd;
	}

	if (byMeasurebit > 0)
	{
		ret = ret + _T("creg c[") + appToString(byMeasurebit) + _T("];\n");
	}
	while (gloablePhaseOpenQASM > PI)
	{
		gloablePhaseOpenQASM -= PI2;
	}
	while (gloablePhaseOpenQASM < -PI)
	{
		gloablePhaseOpenQASM += PI2;
	}
	while (gloablePhaseOpenQiskit > PI)
	{
		gloablePhaseOpenQiskit -= PI2;
	}
	while (gloablePhaseOpenQiskit < -PI)
	{
		gloablePhaseOpenQiskit += PI2;
	}

	CCString sPhase = _T("");
	if (abs(gloablePhaseOpenQASM) > 0.000001f)
	{
		sPhase = _T("// gloable phase (QASM): ") + appToString(gloablePhaseOpenQASM) + _T("\n");
	}
	if (abs(gloablePhaseOpenQiskit) > 0.000001f)
	{
		sPhase = sPhase + _T("// gloable phase (Qiskit): ") + appToString(gloablePhaseOpenQiskit) + _T("\n");
	}
	ret = _T("\nOPENQASM 2.0;\ninclude \"qelib1.inc\";\n// depth: ") + appToString(uiMaxLevel) + _T("\n") + sPhase + ret;

	return ret + operations + _T("\n");
}

QLMatrix QLGate::CreateSingleQubitMatrix(EBasicOperation eop, Real fParam)
{
	QLComplex matrixdata[4];

	switch (eop)
	{
	case EBasicOperation::EBO_H:
		return _hadamard;
	case EBasicOperation::EBO_X:
		return _PauliX;
	case EBasicOperation::EBO_Y:
		return _PauliY;
	case EBasicOperation::EBO_Z:
		return _PauliZ;
	case EBasicOperation::EBO_P:
		matrixdata[0] = _onec;
		matrixdata[1] = _zeroc;
		matrixdata[2] = _zeroc;
		matrixdata[3] = _make_cuComplex(cos(fParam), sin(fParam));
		return QLMatrix::CopyCreate(2U, 2U, matrixdata);
	case EBasicOperation::EBO_RX:
		matrixdata[0] = _make_cuComplex(cos(fParam * F(0.5)), F(0.0));
		matrixdata[1] = _make_cuComplex(F(0.0), sin(fParam * F(-0.5)));
		matrixdata[2] = _make_cuComplex(F(0.0), sin(fParam * F(-0.5)));
		matrixdata[3] = _make_cuComplex(cos(fParam * F(0.5)), F(0.0));
		return QLMatrix::CopyCreate(2U, 2U, matrixdata);
	case EBasicOperation::EBO_RY:
		matrixdata[0] = _make_cuComplex(cos(fParam * F(0.5)), F(0.0));
		matrixdata[1] = _make_cuComplex(sin(fParam * F(0.5)), F(0.0));
		matrixdata[2] = _make_cuComplex(sin(fParam * F(-0.5)), F(0.0));
		matrixdata[3] = _make_cuComplex(cos(fParam * F(0.5)), F(0.0));
		return QLMatrix::CopyCreate(2U, 2U, matrixdata);
	case EBasicOperation::EBO_RZ:
		matrixdata[0] = _make_cuComplex(cos(fParam * F(0.5)), sin(fParam * F(-0.5)));
		matrixdata[1] = _zeroc;
		matrixdata[2] = _zeroc;
		matrixdata[3] = _make_cuComplex(cos(fParam * F(0.5)), sin(fParam * F(0.5)));
		return QLMatrix::CopyCreate(2U, 2U, matrixdata);
	}

	appWarning(_T("CreateSingleQubitMatrix not supported: %d\n"), eop);
	return _I2;
}

CCString QLGate::ToQLISP(const TArray<SBasicOperation>& gates, BYTE byQubit, const TArray<BYTE>& measureAtLast)
{
	BYTE byMeasurebit = 0;
	THashMap<BYTE, BYTE> measuremap;

	QLMatrix single = _I2;
	CCString operations = _T("\n[\n");
	CCString sOpAdd;
	TArray<QLMatrix> singlematrices;
	for (BYTE i = 0; i < byQubit; ++i)
	{
		singlematrices.AddItem(single);
	}

	for (INT i = 0; i < gates.Num(); ++i)
	{
		//UBOOL bIsTwoQubits = FALSE;
		//BYTE dirtyQubit1 = 0;
		//BYTE dirtyQubit2 = 0;

		switch (gates[i].m_eOperation)
		{
		case EBasicOperation::EBO_CC:
			{
				TArray<Real> degrees = GetZYZDecompose(singlematrices[gates[i].m_lstQubits[0]]);
				if (abs(degrees[0]) + abs(degrees[1]) + abs(degrees[2]) > F(0.00001))
				{
					sOpAdd.Format(_T("    ((\'U\', %.5f, %.5f, %.5f), \'Q%d\'),\n"), -degrees[0], -degrees[2], -degrees[1], gates[i].m_lstQubits[0]);
					operations = operations + sOpAdd;
				}
				singlematrices[gates[i].m_lstQubits[0]] = _I2;

				//find a classical bit
				BYTE classicalByte = 0;
				if (measuremap.Exist(gates[i].m_lstQubits[0]))
				{
					classicalByte = measuremap[gates[i].m_lstQubits[0]];
				}
				else
				{
					classicalByte = byMeasurebit;
					measuremap.SetAt(gates[i].m_lstQubits[0], classicalByte);
					byMeasurebit = byMeasurebit + 1;
				}

				//apply measure
				sOpAdd.Format(_T("    ((\'Measure\', %d), \'Q%d\'),\n"), classicalByte, gates[i].m_lstQubits[0]);
				operations = operations + sOpAdd;
			}
			break;
		case EBasicOperation::EBO_H:
		case EBasicOperation::EBO_X:
		case EBasicOperation::EBO_Y:
		case EBasicOperation::EBO_Z:
		case EBasicOperation::EBO_P:
		case EBasicOperation::EBO_RX:
		case EBasicOperation::EBO_RY:
		case EBasicOperation::EBO_RZ:
			{
				singlematrices[gates[i].m_lstQubits[0]] = CreateSingleQubitMatrix(gates[i].m_eOperation, gates[i].m_fClassicalParameter) * singlematrices[gates[i].m_lstQubits[0]];
			}
			break;
		case EBasicOperation::EBO_CX:
			{
				TArray<Real> degrees = GetZYZDecompose(singlematrices[gates[i].m_lstQubits[0]]);
				if (abs(degrees[0]) + abs(degrees[1]) + abs(degrees[2]) > F(0.00001))
				{
					sOpAdd.Format(_T("    ((\'U\', %.5f, %.5f, %.5f), \'Q%d\'),\n"), -degrees[0], -degrees[2], -degrees[1], gates[i].m_lstQubits[0]);
					operations = operations + sOpAdd;
				}
				degrees = GetZYZDecompose(singlematrices[gates[i].m_lstQubits[1]]);
				if (abs(degrees[0]) + abs(degrees[1]) + abs(degrees[2]) > F(0.00001))
				{
					sOpAdd.Format(_T("    ((\'U\', %.5f, %.5f, %.5f), \'Q%d\'),\n"), -degrees[0], -degrees[2], -degrees[1], gates[i].m_lstQubits[1]);
					operations = operations + sOpAdd;
				}
				singlematrices[gates[i].m_lstQubits[0]] = _I2;
				singlematrices[gates[i].m_lstQubits[1]] = _I2;

				sOpAdd.Format(_T("    (\'Cnot\', (\'Q%d\', \'Q%d\')),\n"), gates[i].m_lstQubits[0], gates[i].m_lstQubits[1]);
				operations = operations + sOpAdd;

			}
			break;
		}
	}

	for (BYTE i = 0; i < byQubit; ++i)
	{
		TArray<Real> degrees = GetZYZDecompose(singlematrices[i]);
		if (abs(degrees[0]) + abs(degrees[1]) + abs(degrees[2]) > F(0.00001))
		{
			sOpAdd.Format(_T("    ((\'U\', %.5f, %.5f, %.5f), \'Q%d\'),\n"), -degrees[0], -degrees[2], -degrees[1], i);
			operations = operations + sOpAdd;
		}
	}

	for (INT i = 0; i < measureAtLast.Num(); ++i)
	{
		BYTE classicalByte = 0;
		if (measuremap.Exist(gates[i].m_lstQubits[0]))
		{
			classicalByte = measuremap[gates[i].m_lstQubits[0]];
		}
		else
		{
			classicalByte = byMeasurebit;
			measuremap.SetAt(gates[i].m_lstQubits[0], classicalByte);
			byMeasurebit = byMeasurebit + 1;
		}

		//apply measure
		sOpAdd.Format(_T("    ((\'Measure\', %d), \'Q%d\'),\n"), classicalByte, gates[i].m_lstQubits[0]);
		operations = operations + sOpAdd;
	}

	return operations + _T("]\n");
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================