#ifndef _STIFFMATRIX2_H_
#define _STIFFMATRIX2_H_
#include <BilinearOperator.h>
#include <DGFEMSpace.h>
#include <lac/vector.h>

class StiffMatrix2 : public StiffMatrix<2,double>
{
    public:
    StiffMatrix2() {};
    StiffMatrix2(DGFEMSpace<double,2>& dgsp):StiffMatrix<2,double>(dgsp),dg_fem_space(&dgsp) {};
    virtual ~StiffMatrix2() {};
    private:
    DGFEMSpace<double,2>* dg_fem_space;
    public:
    virtual void getElementMatrix(const Element<double,2>&,
                                  const Element<double,2>&,
                                  const ActiveElementPairIterator<2>::State state);
    virtual void getDGElementMatrix(const DGElement<double,2>&,
                                    const Element<double,2>&); 
    DGFEMSpace<double,2> * DGFEMSpace0() {return dg_fem_space;};
    const DGFEMSpace<double,2> * DGFEMSpace0() const {return dg_fem_space;};
};

void L2DGDiscretize(double (*f0)(const double *), 
		            const DGFEMSpace<double,2>& dgfem_space,
                    Vector<double>& f1, 
                    int algebric_accuracy);

#endif
