#include <Matrix.h>
#include <FEMSpace.h>
#include <DGFEMSpace.h>

void StiffMatrix2::getElementMatrix(const Element<double,2>& element0,
                                    const Element<double,2>& element1,
                                    const ActiveElementPairIterator<2>::State state)
{
  const std::vector<int>& ele_dof0 = element0.dof();
  const std::vector<int>& ele_dof1 = element1.dof();
  int n_element_dof0 = ele_dof0.size();
  int n_element_dof1 = ele_dof1.size();
  double volume = element0.templateElement().volume();
  const int& alg_acc = this->algebricAccuracy();
  const QuadratureInfo<2>& quad_info = element0.findQuadratureInfo(alg_acc);
  std::vector<double> jacobian = element0.local_to_global_jacobian(quad_info.quadraturePoint());
  int n_quadrature_point = quad_info.n_quadraturePoint();
  std::vector<AFEPack::Point<2> > q_point = element0.local_to_global(quad_info.quadraturePoint());
  std::vector<std::vector<std::vector<double> > > basis_gradient = element0.basis_function_gradient(q_point);
  for (int l = 0;l < n_quadrature_point;l ++) {
    double Jxw = quad_info.weight(l)*jacobian[l]*volume;
    for (int j = 0;j < n_element_dof0;j ++) {
      for (int k = 0;k < n_element_dof1;k ++) {
	      this->elementMatrix(j,k) += Jxw*innerProduct(basis_gradient[j][l],basis_gradient[k][l]);
      }
    }
  }
  GeometryBM& geo = element0.geometry();
  for(int i = 0; i<geo.n_boundary(); i++)
  {
      int index = geo.boundary(i);
      const DGElement<double,2>& the_dgele = DGFEMSpace0()->dgElement(index);
      if(the_dgele.p_neighbourElement(0)==NULL || the_dgele.p_neighbourElement(1)==NULL)
          getDGElementMatrix(the_dgele,element0);
  }
}


void StiffMatrix2::getDGElementMatrix(const DGElement<double,2>& dgele0, const Element<double,2>& element0)
{
    int n_dof_ele0 = element0.n_dof();
    double volume = dgele0.templateElement().volume();
    const int& alg_acc = this->algebricAccuracy();
    const QuadratureInfo<1>& quad_info = dgele0.findQuadratureInfo(alg_acc);
    std::vector<double> jacobian = dgele0.local_to_global_jacobian(quad_info.quadraturePoint());
    int n_quadrature_point = quad_info.n_quadraturePoint();
    std::vector<AFEPack::Point<2> > q_point = dgele0.local_to_global(quad_info.quadraturePoint());
    std::vector<std::vector<double>> basis_value0 = element0.basis_function_value(q_point);
    //std::vector<std::vector<std::vector<double>>> basis_gradient0 = element0.basis_function_gradient(q_point);
    //std::vector<std::vector<double>> out_normal = unitOutNormal(q_point,element0,dgele0);
    for(int l = 0;l<n_quadrature_point;l++)
    {
        double Jxw = quad_info.weight(l)*jacobian[l]*volume;
        double h = Jxw;
        for(int j = 0;j<n_dof_ele0;j++)
            for(int k = 0;k<n_dof_ele0;k++)
            {
                this->elementMatrix(j,k) += Jxw*100*(basis_value0[j][l]*basis_value0[k][l])/h;
                //this->elementMatrix(j,k) -= Jxw*innerProduct(basis_gradient0[j][l],out_normal[l])*basis_value0[k][l];
                //this->elementMatrix(j,k) -= Jxw*innerProduct(basis_gradient0[k][l],out_normal[l])*basis_value0[j][l];
            }
    }
}


/*void StiffMatrix2::buildDofInfo2()
{
  nDof0() = fem_space0->n_dof();
  nDof1() = fem_space1->n_dof();
  std::vector<int> n_coupling_dof(nDof0(), 0);
  if ((void *)fem_space0 == (void *)fem_space1) {
    DGElementIterator the_dgelement0 = fem_space0->beginDGElement();
    DGElementIterator end_dgelement0 = fem_space0->endDGElement();
    //DGElementIterator the_dgelement1 = fem_space1->beginElement();
    for (;the_dgelement0 != end_dgelement0;the_dgelement0 ++) {
      const Element<double,2>& nei0 = (*the_dgelement0).neighbourElement(0);
      const Element<double,2>& nei1 = (*the_dgelement0).neighbourElement(1);
      getElementPattern(nei0,nei1);
      int n_element_dof0 = elementDof0().size();
      int n_element_dof1 = elementDof1().size();
      for (int i = 0;i < n_element_dof0;i ++)
	      n_coupling_dof[elementDof0(i)] += (n_element_dof0+n_element_dof1);
      for (int i = 0;i < n_element_dof1;i ++)
          n_coupling_dof[elementDof1(i)] += (n_element_dof0+n_element_dof1);
    }
  }
  else
      std::cout<<"The fem_space must be same!"<<std::endl;
  nMaxCouplingDof() = *max_element(n_coupling_dof.begin(), n_coupling_dof.end());
  if (nMaxCouplingDof() > nDof1()) nMaxCouplingDof() = nDof1();
}*/
/*
void StiffMatrix2::buildSparsityPattern()
{
    buildDofInfo();
    sparsity_pattern.reinit(nDof0(), nDof1(), nMaxCouplingDof());
    if ((void *)fem_space0 == (void *)fem_space1) {
    //单元内部
    typename std::vector<Element<value_type0,DIM,DOW,TDIM0> >::iterator the_element0 = fem_space0->beginElement();
    typename std::vector<Element<value_type0,DIM,DOW,TDIM0> >::iterator end_element0 = fem_space0->endElement();
    typename std::vector<Element<value_type1,DIM,DOW,TDIM1> >::iterator the_element1 = fem_space1->beginElement();
    for (;the_element0 != end_element0;the_element0 ++, the_element1 ++) {
      getElementPattern(*the_element0, *the_element1);
      addElementPattern();
    }
    //边界
    DGElementIterator the_dgelement0 = fem_space0->beginDGElement();
    DGElementIterator end_dgelement0 = fem_space0->endDGElement();
    for(;the_dgelement0 != end_dgelement0; the_dgelement0++){
        const Element<double,2>& nei0 = (*the_dgelement0).neighbourElement(0);
        const Element<double,2>& nei1 = (*the_dgelement0).neighbourElement(1);
        getElementPattern(nei0,nei1);
        addElementPattern();
    }
  }
  else
    std::cout<<"The fem_space must be same!"<<std::endl;
  sparsity_pattern.compress();
}
*/
/*
void StiffMatrix2::buildSparseMatrix()
{
    SparseMatrix<double>::reinit(sparsity_pattern);
    if ((void *)fem_space0 == (void *)fem_space1) {
    //单元内部
    typename std::vector<Element<value_type0,DIM,DOW,TDIM0> >::iterator the_element0 = fem_space0->beginElement();
    typename std::vector<Element<value_type0,DIM,DOW,TDIM0> >::iterator end_element0 = fem_space0->endElement();
    typename std::vector<Element<value_type1,DIM,DOW,TDIM1> >::iterator the_element1 = fem_space1->beginElement();
    for (;the_element0 != end_element0;the_element0 ++, the_element1 ++) {
         getElementPattern(*the_element0, *the_element1);
         elementMatrix().reinit(elementDof0().size(), elementDof1().size());
         getElementMatrix(*the_element0, *the_element1);
         addElementMatrix();
        }
    //边界
    DGElementIterator the_dgelement0 = fem_space0->beginDGElement();
    DGElementIterator end_dgelement0 = fem_space0->endDGElement();
    DGElementIterator the_dgelement1 = fem_space1->beginDGElement();
    for (;the_dgelement0 != end_dgelement0; the_dgelement0++,the_dgelement1++)
    {
        const Element<double,2>& nei0 = (*the_dgelement0).neighbourElement(0);
        const Element<double,2>& nei1 = (*the_dgelement0).neighbourElement(1);
        getElementPattern(nei0,nei1);
        elementMatrix().reinit(elementDof0().size(),elementDof1().size());
        getDGElementMatrix(*the_dgelement0, *the_dgelement1);
        addElementMatrix();
    }
    }
}
*/
