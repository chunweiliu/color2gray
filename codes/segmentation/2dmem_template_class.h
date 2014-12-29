#ifndef __2DMEM_TEMPLATE_CLASS_H__
#define __2DMEM_TEMPLATE_CLASS_H__

/**********************************************
* Add m_bFree for directly assigning a memory *
* created from other part.                    *
* The memory should be freed from other part. *
*                                             *
*                 Written by Kai-Yueh Chang   *
*                 Released 2008/11/25 11:53   *
**********************************************/

/**********************************************
* Add the operator[] for directly accessing   *
*                                             *
*                 Written by Kai-Yueh Chang   *
*                 Released 2006/04/19 07:40   *
**********************************************/

/**********************************************
* For the problem of the function reveal in   *
* VC++ 6.0, the functions are moved into      *
* class bady.                                 *
*                                             *
* For acceleration, some functions are        *
* changed to be inline function.              *
*                                             *
*                 Written by Kai-Yueh Chang   *
*                 Released 2004/04/06 20:02   *
**********************************************/

#include <stdlib.h>

template<class T>
class TC2DMem{
public:
	int m_row;
	int m_col;
	T** m_ppMem;
private:
	bool m_bFree;
public:
	TC2DMem(){
		m_row = 0;
		m_col = 0;
		m_bFree = 0;
		m_ppMem = NULL;
	}
	~TC2DMem(){
		Free2DMem();
	}

	T** Alloc2DMem(int nRow, int nCol){
		T** ppMem;
		T* pMem;

		T** ppCur;
		T* pCur;

		int tmp;

		Free2DMem();
		if(nRow<1)
			return NULL;
		if(nCol<1)
			return NULL;
		m_row = nRow;
		m_col = nCol;

		ppMem = (T**)malloc(m_row*sizeof(T*));
		pMem = new T[m_row*m_col];
		ppCur = ppMem;
		pCur = pMem;

		for(tmp=0; tmp<m_row; tmp++){
			*ppCur = pCur;
			pCur += m_col;
			ppCur++;
		}

		m_ppMem = ppMem;

		m_bFree = 1;
		return m_ppMem;
	}
	T** AssignMemory(int nRow, int nCol, T* pMem){
		T** ppMem;

		T** ppCur;
		T* pCur;

		int tmp;

		Free2DMem();
		if(nRow<1)
			return NULL;
		if(nCol<1)
			return NULL;
		m_row = nRow;
		m_col = nCol;

		ppMem = (T**)malloc(m_row*sizeof(T*));
		ppCur = ppMem;
		pCur = pMem;

		for(tmp=0; tmp<m_row; tmp++){
			*ppCur = pCur;
			pCur += m_col;
			ppCur++;
		}

		m_ppMem = ppMem;

		return m_ppMem;
	}
	void Free2DMem(){
		if(NULL!=(m_ppMem)){
			if(m_bFree){
				delete [](*m_ppMem);
				m_bFree = 0;
			}
			free(m_ppMem);
			m_ppMem = NULL;
		}
		m_row = 0;
		m_col = 0;
	}

	inline void GetSize(int* pRow, int* pCol){
		*pRow = m_row;
		*pCol = m_col;
	}
	inline T* Get1DPointer(){
		if(NULL==m_ppMem)
			return NULL;
		return *m_ppMem;
	}
	inline T** Get2DPointer(){
		return m_ppMem;
	}
public:
	inline T* operator[](int nRow) const{
		return m_ppMem[nRow];
	}
};
/*
template<class T>
TC2DMem<T>::TC2DMem(){
	m_row = 0;
	m_col = 0;
	m_ppMem = NULL;
}


template<class T>	
TC2DMem<T>::~TC2DMem(){
	Free2DMem();
}

template<class T>
T** TC2DMem<T>::Alloc2DMem(int nRow, int nCol){
	T** ppMem;
	T* pMem;

	T** ppCur;
	T* pCur;

	int tmp;

	if(nRow<1)
		return NULL;
	if(nCol<1)
		return NULL;
	m_row = nRow;
	m_col = nCol;
	Free2DMem();

	ppMem = (T**)malloc(m_row*sizeof(T*));
	pMem = new T[m_row*m_col];
	ppCur = ppMem;
	pCur = pMem;

	for(tmp=0; tmp<m_row; tmp++){
		*ppCur = pCur;
		pCur += m_col;
		ppCur++;
	}

	m_ppMem = ppMem;

	return m_ppMem;
}

template<class T>
void TC2DMem<T>::Free2DMem(){
	if(NULL!=(m_ppMem)){
		delete [](*m_ppMem);
		free(m_ppMem);
		m_ppMem = NULL;
	}
}

template<class T>
void TC2DMem<T>::GetSize(int* pRow, int* pCol){
	*pRow = m_row;
	*pCol = m_col;
}

template<class T>
T* TC2DMem<T>::Get1DPointer(){
	if(NULL==m_ppMem)
		return NULL;
	return *m_ppMem;
}

template<class T>
T** TC2DMem<T>::Get2DPointer(){
	return m_ppMem;
}
*/
#endif