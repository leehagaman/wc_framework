#include <TMatrixD.h>
#include <TMatrixDSymEigen.h>

int test() {
  TFile *file = new TFile("cov_17.root");
  TMatrixD* mtest = (TMatrixD*)(file->Get("frac_cov_xf_mat_17"));

  mtest->Draw("colz");
  
  TMatrixD myMatrix(*mtest);
    //myMatrix(0, 0) = 0;     myMatrix(0, 1) = 1;     myMatrix(0, 2) = 2;
    //myMatrix(1, 0) = 1;     myMatrix(1, 1) = 4;     myMatrix(1, 2) = -0.3;
    //myMatrix(2, 0) = 2;     myMatrix(2, 1) = -0.3;  myMatrix(2, 2) = 6;

    // for (Int_t i=0;i!=100;i++){
    //   for (Int_t j=0;j!=100;j++){
    // 	myMatrix(i,j) = (*mtest)(i,j);
    //   }
    // }
    
    
    if (!myMatrix.IsSymmetric()) {
        std::cout << "Matrix is not symmetric!\n";
        myMatrix = 0.5 * (myMatrix + myMatrix.T()); // Symmetrization
    }


    //    myMatrix.Print();

    TMatrixDSym symMatrix(myMatrix.GetNrows());
    for (Int_t i = 0; i < symMatrix.GetNrows(); ++i) {
        for (Int_t j = 0; j <= i; ++j) {
            symMatrix(i, j) = myMatrix(i, j);
        }
    }

    // Step 2: Get the eigenvalues and eigenvectors
    TMatrixDSymEigen eig(symMatrix);
    TVectorD eigenvalues = eig.GetEigenValues();

    std::cout << "eigenvalues:\n";
    eigenvalues.Print();

    
    // Step 3: Swap eigenvalues with their absolute values
    for (Int_t i = 0; i < eigenvalues.GetNrows(); ++i) {
        Double_t eigenvalue = eigenvalues[i];
        //eigenvalues[i] = TMath::Abs(eigenvalue);
    }
    std::cout << "absoluted eigenvalues:\n";
    eigenvalues.Print();
    

    // Step 4: Reconstruct the modified matrix
    TMatrixD diagEigenvalues(eigenvalues.GetNrows(), eigenvalues.GetNrows());
    diagEigenvalues.Zero();
    for (Int_t i = 0; i < eigenvalues.GetNrows(); ++i) {
        diagEigenvalues(i, i) = eigenvalues[i];
    }

    TMatrixD eigenvectorsCopy = eig.GetEigenVectors();
    TMatrixD modifiedMatrix = eigenvectorsCopy * diagEigenvalues * eigenvectorsCopy.T();

    // modifiedMatrix.Print();

    cout << "After building from eigenvalues:\n";

    
    TMatrixDSym symMatrix2(modifiedMatrix.GetNrows());
    for (Int_t i = 0; i < modifiedMatrix.GetNrows(); ++i) {
            for (Int_t j = 0; j <= i; ++j) {
                    symMatrix2(i, j) = myMatrix(i, j);
            }
    }

    // get the eigenvalues and eigenvectors
    TMatrixDSymEigen eig2(symMatrix2);
    TVectorD eigenvalues2 = eig2.GetEigenValues();
    TMatrixD eigenvectors2 = eig2.GetEigenVectors();

    // print negatives
    for (Int_t i = 0; i < eigenvalues2.GetNrows(); ++i) {
            Double_t eigenvalue2 = eigenvalues2[i];
            cout << i << "  eigenvalue: " << eigenvalue2 << " " << eigenvalues[i] << "\n";
    }


    return 0;
}
