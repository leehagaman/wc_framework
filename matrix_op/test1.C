#include <TMatrixD.h>
#include <TMatrixDSymEigen.h>

void test1(){

  TFile *file = new TFile("cov_17.root");
  TMatrixD* mtest = (TMatrixD*)(file->Get("frac_cov_xf_mat_17"));
  
  mtest->Draw("colz");
  
  cout << "After modifying Xs matrix:\n";
  
  // copy un-collapsed Xs matrix, make sure it's symmetrical
  TMatrixD myMatrix = *mtest;
  /*if (!myMatrix.IsSymmetric()) {
    std::cout << "    Matrix is not symmetric!\n";
    myMatrix = 0.5 * (myMatrix + myMatrix.T()); // Symmetrization
    }*/
  
  // put it into a symmetric matrix type
  TMatrixDSym symMatrix(myMatrix.GetNrows());
  for (Int_t i = 0; i < symMatrix.GetNrows(); ++i) {
    for (Int_t j = 0; j <= i; ++j) {
      symMatrix(i, j) = myMatrix(i, j);
    }
  }

  // get the eigenvalues and eigenvectors
  TMatrixDSymEigen eig(symMatrix);
  TVectorD eigenvalues = eig.GetEigenValues();
  TMatrixD eigenvectors = eig.GetEigenVectors();
  
  // swap eigenvalues with their absolute values, print them
  for (Int_t i = 0; i < eigenvalues.GetNrows(); ++i) {
    Double_t eigenvalue = eigenvalues[i];
    //if (eigenvalue < 0) cout << "    negative eigenvalue: " << eigenvalue << "\n";        
    cout << "    " << i << " eigenvalue: " << eigenvalue << "\n";
    //eigenvalues[i] = TMath::Abs(eigenvalue);
  }
  
  // making diagonal eigenvalues matrix
  TMatrixD diagEigenvalues(eigenvalues.GetNrows(), eigenvalues.GetNrows());
  diagEigenvalues.Zero();
  for (Int_t i = 0; i < eigenvalues.GetNrows(); ++i) {
    diagEigenvalues(i, i) = eigenvalues[i];
  }
  
  // getting the new full matrix
  TMatrixD modifiedMatrix = eigenvectors * diagEigenvalues * eigenvectors.T();
  //  (*map_matrix_flux_Xs_frac[idx]) = modifiedMatrix;
  
  cout << "After putting it back:\n";
  
  // put it into a symmetric matrix type
  TMatrixDSym symMatrix2(modifiedMatrix.GetNrows());
  for (Int_t i = 0; i < modifiedMatrix.GetNrows(); ++i) {
    for (Int_t j = 0; j <= i; ++j) {
      symMatrix2(i, j) = modifiedMatrix(i, j);
    }
  }
  
  // get the eigenvalues and eigenvectors
  TMatrixDSymEigen eig2(symMatrix2);
  TVectorD eigenvalues2 = eig2.GetEigenValues();
  TMatrixD eigenvectors2 = eig2.GetEigenVectors();
  
  // print negatives
  for (Int_t i = 0; i < eigenvalues2.GetNrows(); ++i) {
    Double_t eigenvalue2 = eigenvalues2[i];
    //if (eigenvalue < 0) cout << "    negative eigenvalue: " << eigenvalue << "\n";        
    cout << "    " << i << " eigenvalue: " << eigenvalue2 << " " << eigenvalues[i] << "\n";
  }
}
