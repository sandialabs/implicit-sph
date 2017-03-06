#include "NOX.H"
#include "NOX_Epetra.H"

namespace LAMMPS_NS {

  class ProblemScalar : public NOX::Epetra::Interface::Required,
                        public NOX::Epetra::Interface::Jacobian,
                        public NOX::Epetra::Interface::Preconditioner {
  public:
    ProblemScalar() { }
    virtual~ProblemScalar() { }
    
    bool computeF(const Epetra_Vector &x, Epetra_Vector &f,
                  const FillType flag = Residual); 
    bool computeJacobian(const Epetra_Vector &x, Epetra_Operator &J);
    bool computePreconditioner(const Epetra_Vector &x, Epetra_Operator &M,
                               Teuchos::ParameterList *param = NULL);
  };

  // impl.
  inline bool 
  ProblemScalar::computeF(const Epetra_Vector &x, Epetra_Vector &f,
                          const FillType flag) {
    f[0] = x[0]*x[0] - 4.0;
    return true;
  }
  inline bool 
  ProblemScalar::computeJacobian(const Epetra_Vector &x, 
                                 Epetra_Operator &J) {
    double val = 2.0;
    int idx = 0;
    auto *jacobian = dynamic_cast<Epetra_CrsMatrix*>(&J);
    jacobian->InsertMyValues(0, 1, &val, &idx); 

    return true;
  }
  inline bool 
  ProblemScalar::computePreconditioner(const Epetra_Vector &x, Epetra_Operator &M,
                                       Teuchos::ParameterList *param) {
    return false;
  }

}
