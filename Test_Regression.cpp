#include "FHE-SI.h"
#include "Ciphertext.h"
#include "SingleCRT.h"
#include "Regression.h"
#include <time.h>
#include <ctime>
#include <fstream>
#include "Util.h"

static int RunRegressionTest(Matrix<ZZ> &rawData, vector<ZZ> &labels, 
                             const ZZ &p, const FHEcontext &context) {
                              
  vector<ZZ> theta;
  ZZ det;
  RegressPT(theta, det, rawData, labels);

  cout << "Expected values: " << endl;
  for (unsigned i = 0; i < theta.size(); i++) {
      cout << "  theta[" << i << "] = " << theta[i] % p << endl;
  }
  cout << "  Determinant: " << det % p << endl;
  cout << endl << endl;
  
  double start(clock());
  Regression regress(context);
  cout << "Setup time: " << (clock()-start)/CLOCKS_PER_SEC << endl;
  
  vector<vector<Plaintext>> ptxtData;
  vector<Plaintext> ptxtLabels;
  
  double batchTime = BatchData(ptxtData, ptxtLabels, rawData, labels, context);
  
  cout << "Batch time: " << batchTime << endl;
  
  double encStart(clock());
  regress.AddData(ptxtData, ptxtLabels);
  cout << "Encryption time: " << (clock()-encStart)/CLOCKS_PER_SEC << endl;
  
  vector<Ciphertext> encTheta;
  Ciphertext encDet(regress.GetPublicKey());
  
  double regressionStart(clock());
  regress.Regress(encTheta, encDet);
  cout << "Regression time: " << (clock()-regressionStart)/CLOCKS_PER_SEC << endl;
  
  FHESISecKey secretKey = regress.GetSecretKey();
  Plaintext tmp(context);
  vector<ZZ_pX> msgs;
  ZZ_pX msg;
  
  double decStart(clock());
  cout << endl << "Computed values: " << endl;
  for (unsigned i = 0; i < encTheta.size(); i++) {
    secretKey.Decrypt(tmp, encTheta[i]);
    tmp.DecodeSlots(msgs);
    cout << "  theta[" << i << "] = " << msgs[0] << endl;
  }

  secretKey.Decrypt(tmp, encDet);
  tmp.DecodeSlots(msgs);

  // Output keys
  cout << "Output secret key\n";
  std::ofstream ofs_secret_key;
  ofs_secret_key.open("secret_key.txt", std::ofstream::out | std::ofstream::app);
  secretKey.Export(ofs_secret_key);
  ofs_secret_key.close();

  cout << "Output public key\n";
  FHESIPubKey publicKey = regress.GetPublicKey();
  std::ofstream ofs_public_key;
  ofs_public_key.open("public_key.txt", std::ofstream::out | std::ofstream::app);
  publicKey.Export(ofs_public_key);
  ofs_public_key.close();

  cout << "Output encypted data\n";
  std::ofstream ofs_enc_data;
  ofs_enc_data.open("encrypted_data.txt", std::ofstream::out | std::ofstream::app);
  ofs_enc_data << regress.data;
  ofs_enc_data.close();




  
  cout << "  Determinant: " << msgs[0] << endl << endl;
  cout << "Decryption time: " << (clock() - decStart)/CLOCKS_PER_SEC << endl;
  cout << "Total time: " << (clock()-start)/CLOCKS_PER_SEC << endl;
  
  return 0;
}

int main(int argc, char *argv[]) {
  srand48(time(NULL));
  SetSeed(to_ZZ(time(NULL)));

  unsigned p, g, logQ = 0;
  char *datafile;
  
  if (argc >= 4) {
    datafile = argv[1];
    p = atoi(argv[2]);
    g = atoi(argv[3]);
  } else {
    cout << "usage: Test_Regression_x datafile p generator" << endl;
    return 1;
  }
  
  unsigned blockSize = 1;
  unsigned val = (p-1)/2-1;
  while (val > 1) {
    blockSize <<= 1;
    val >>= 1;
  }
  
  Matrix<ZZ> rawData;
  vector<ZZ> labels;
  unsigned dim;
  
  if (!LoadData(rawData, labels, dim, datafile)) {
    return 1;
  }

  // Write out matrix
  cout << "Output raw matrix\n";
  std::ofstream ofs_raw_data;
  ofs_raw_data.open("raw_data.txt", std::ofstream::out | std::ofstream::app);
  ofs_raw_data << rawData;
  ofs_raw_data.close();
  
  unsigned n = (p-1)/2-1;
  unsigned nBlocks = labels.size() / blockSize;
  if (labels.size() % blockSize != 0) {
    nBlocks++;
  }
  unsigned xi = max(nBlocks, dim);

  double lgQ = 4.5*log(n)+max(1,(int) dim-1)*(log(1280)+2*log(n)+log(xi));
  logQ = (unsigned) ceil(lgQ / log(2) + 24.7);

  cout << "================================================" << endl
       << "Running regression tests using Brakerski system." << endl
       << "================================================" << endl;
  
  cout << "Parameters: " << endl
   << "  data file: " << datafile << endl
   << "  logQ: " << logQ << endl
   << "  p: " << p << endl
   << "  generator: " << g << endl
   << "  block size: " << blockSize << endl
   << "  num blocks: " << nBlocks << endl;
   
  FHEcontext context(p-1, logQ, p, g, 3);
  activeContext = &context;
  
  cout << "Running " << dim << "-dimensional regression on " << rawData.NumRows()
       << " datapoints in " << (rawData.NumRows()+blockSize-1)/blockSize
       << " blocks, modulo prime " << p << endl;
  
  context.SetUpSIContext(xi);
  return RunRegressionTest(rawData, labels, to_ZZ(p), context);
}
