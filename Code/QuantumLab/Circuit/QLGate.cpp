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
{
	SBasicOperation op_h;
	op_h.m_eOperation = eOp;

	switch (eOp)
	{
	case EBasicOperation::EBO_H:
		m_lstQubits.push_back(0);
		op_h.m_lstQubits.push_back(0);
		m_lstOperations.push_back(op_h);
		m_sName = "H";
		break;
	default:
		break;
	}
}

QLGate::~QLGate()
{

}

QLGate* QLGate::Controlled(BYTE controlledQubitCount) const
{
	return NULL;
}

QLGate* QLGate::Dagger() const
{
	return NULL;
}

void QLGate::AppendGate(const QLGate& toAppend, const std::vector<BYTE>& lstMappingQubits, const std::vector<SClassicalParamters>& addtoinalParameters)
{

}

std::vector<SBasicOperation> QLGate::GetOperation(const std::vector<BYTE>& lstMappingQubits) const
{
	std::vector<SBasicOperation> ret;
	if (m_bBasicOperation)
	{
		SIZE_T opsize = m_lstOperations.size();
		for (INT i = 0; i < opsize; ++i)
		{
			SBasicOperation newone;
			newone.m_eOperation = m_lstOperations[i].m_eOperation;
			SIZE_T qubitsize = m_lstOperations[i].m_lstClassicalParameters.size();
			for (INT j = 0; j < opsize; ++j)
			{
				newone.m_lstQubits.push_back(lstMappingQubits[m_lstOperations[i].m_lstClassicalParameters[j]]);
			}
			SIZE_T paramsize = m_lstOperations[i].m_lstQubits.size();
			for (INT j = 0; j < paramsize; ++j)
			{
				newone.m_lstClassicalParameters.push_back(m_lstOperations[i].m_lstQubits[j]);
			}
			ret.push_back(newone);
		}
	}
	else
	{
		SIZE_T gatesize = m_lstSubGates.size();
		for (INT i = 0; i < gatesize; ++i)
		{
			std::vector<SBasicOperation> oplist = m_lstSubGates[i]->GetOperation(lstMappingQubits);
			ret.insert(ret.end(), oplist.begin(), oplist.end());
		}
	}
	return ret;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================