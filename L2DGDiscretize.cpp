#include "Matrix.h"

void L2DGDiscretize(double (*f0)(const double *), 
		            const DGFEMSpace<double,2>& dgfem_space,
                    Vector<double>& f1, 
                    int algebric_accuracy)
{
   typename DGFEMSpace<double,2>::ConstDGElementIterator the_dgele = dgfem_space.beginDGElement(),
   end_dgele = dgfem_space.endDGElement();
   for(;the_dgele != end_dgele; the_dgele++)
    {
        if((the_dgele->p_neighbourElement(1))==NULL) 
        {
        const Element<double,2>& nei = the_dgele->neighbourElement(0);
        double volume = the_dgele->templateElement().volume();
        const QuadratureInfo<1>& quad_info = the_dgele->findQuadratureInfo(algebric_accuracy);
        std::vector<double> jacobian = the_dgele->local_to_global_jacobian(quad_info.quadraturePoint());
        int n_quadrature_point = quad_info.n_quadraturePoint();
        std::vector<AFEPack::Point<2>> q_point = the_dgele->local_to_global(quad_info.quadraturePoint());
        std::vector<std::vector<double>> basis_value = nei.basis_function_value(q_point);
        const std::vector<int>& element_dof = nei.dof();
        int n_element_dof = element_dof.size();
        for (int l = 0;l < n_quadrature_point;l ++) {
            double f0_value = (*f0)(q_point[l]);
            double Jxw = quad_info.weight(l)*jacobian[l]*volume;
            double h = Jxw;
            for (int j = 0;j < n_element_dof;j ++) {
	            f1(element_dof[j]) += Jxw*100*f0_value*basis_value[j][l]/h;
            }
        }
        }
        else continue;
    }
}
