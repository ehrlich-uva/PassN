void CrvPedestal(const std::string &inputFileName, const std::string &outputFileName)
{
    TFile *inputFile = TFile::Open(inputFileName.c_str(),"update");
    inputFile->cd("CrvPedestalFinder");
    TList *keys = gDirectory->GetListOfKeys();

    TF1 funcPedestal("f0", "gaus");

    std::ofstream outputFile;
    outputFile.open(outputFileName);
    outputFile<<"TABLE CRVSiPM"<<std::endl;
    outputFile<<"#channel, pedestal, calibPulseHeight, calibPulseArea"<<std::endl;

    for(int iKey=0; iKey<keys->GetSize(); ++iKey)
    {
      TKey *key=(TKey*)keys->At(iKey);
      if(strncmp(key->GetName(),"crvPedestalHist_",16)!=0) continue;
      TH1F *hist=(TH1F*)key->ReadObj();
      int channel=atoi(key->GetName()+16);

      if(hist->GetEntries()<100) //not enough data
      {
        outputFile<<channel<<","<<0<<",-1,-1"<<std::endl;
        continue;
      }

/*
      int n=hist->GetNbinsX();
      double overflow=hist->GetBinContent(0)+hist->GetBinContent(n+1);
      if(overflow/((double)hist->GetEntries())>0.1) //too much underflow/overflow. something may be wrong.
      {
        outputFile<<channel<<","<<0<<",-1,-1"<<std::endl;
        continue;
      }
*/

      int maxbinPedestal = hist->GetMaximumBin();
      double peakPedestal = hist->GetBinCenter(maxbinPedestal);
      funcPedestal.SetRange(peakPedestal-4,peakPedestal+4);
      funcPedestal.SetParameter(1,peakPedestal);
      hist->Fit(&funcPedestal, "0QR");
      outputFile<<channel<<","<<funcPedestal.GetParameter(1)<<",-1,-1"<<std::endl;  //only print out pedestal values
    }

    outputFile.close();

    inputFile->Write(0,TFile::kWriteDelete);
    inputFile->Close();
}

