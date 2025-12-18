const double fitRangeStart=0.8;
const double fitRangeEnd=1.2;
const int    minHistEntries=100;
const int    spectrumNPeaks=6;
const double spectrumPeakSigma=4.0;
const double spectrumPeakThreshold=0.01;
const double peakRatioTolerance=0.2;

bool FindSPEpeak(TH1F *hist, TSpectrum &spectrum, double &SPEpeak);

void CrvCalibration(const std::string &inputFileName, const std::string &outputFileName)
{
    TFile *inputFile = TFile::Open(inputFileName.c_str(),"update");
    inputFile->cd("CrvCalibration");
    TTree *treePedestals = (TTree*)gDirectory->FindObjectAny("crvPedestals");
    size_t channel;
    double pedestal;
    treePedestals->SetBranchAddress("channel", &channel);
    treePedestals->SetBranchAddress("pedestal", &pedestal);
    std::map<size_t,double> pedestals;  //first need to fill this map to filter out multiple entries of the same channel due to $ROOTSYS/bin/hadd
    for(int i=0; i<treePedestals->GetEntries(); ++i)
    {
      treePedestals->GetEntry(i);
      pedestals[channel]=pedestal;
    }

    TF1 funcCalib("SPEpeak", "gaus");
    TSpectrum spectrum(spectrumNPeaks); //any value of 3 or less results in a "Peak buffer full" warning.

    std::ofstream outputFile;
    outputFile.open(outputFileName);
    outputFile<<"TABLE CRVSiPM "<<std::endl;
    outputFile<<"#channel, pedestal, calibPulseHeight, calibPulseArea"<<std::endl;

    for(auto iter=pedestals.begin(); iter!=pedestals.end(); ++iter)
    {
      channel=iter->first;
      pedestal=iter->second;

      TH1F *hist;
      double calibValue[2];
      for(int i=0; i<2; ++i) //loop over hisograms with pulse areas and pulse heights
      {
        if(i==1) hist=(TH1F*)gDirectory->FindObjectAny(Form("crvCalibrationHistPulseArea_%zu",channel));
        else hist=(TH1F*)gDirectory->FindObjectAny(Form("crvCalibrationHistPulseHeight_%zu",channel));

        double peakCalib=0;
        if(!FindSPEpeak(hist, spectrum, peakCalib))
        {
          calibValue[i]=-1;
          continue;
        }

        funcCalib.SetRange(peakCalib*fitRangeStart,peakCalib*fitRangeEnd);
        if(hist->FindBin(peakCalib*fitRangeStart)==hist->FindBin(peakCalib*fitRangeEnd)) //fit range start/end are in the same bin
        {
          calibValue[i]=-1;
          continue;
        }
        funcCalib.SetParameter(1,peakCalib);
        hist->Fit(&funcCalib, "QR");
        calibValue[i]=funcCalib.GetParameter(1);
      }

      outputFile<<channel<<","<<pedestal<<","<<calibValue[0]<<","<<calibValue[1]<<std::endl;
    }

    outputFile<<std::endl;
    inputFile->Write(0,TFile::kWriteDelete);

    //time offsets
    TTree *treeTimeOffsets = (TTree*)gDirectory->FindObjectAny("crvTimeOffsets");
    double offset;
    treeTimeOffsets->SetBranchAddress("channel", &channel);
    treeTimeOffsets->SetBranchAddress("timeOffset", &offset);
    std::map<size_t,double> timeOffsets;  //first need to fill this map to filter out multiple entries of the same channel due to $ROOTSYS/bin/hadd
    for(int i=0; i<treeTimeOffsets->GetEntries(); ++i)
    {
      treeTimeOffsets->GetEntry(i);
      timeOffsets[channel]=offset;
    }

    outputFile<<"TABLE CRVTime"<<std::endl;
    outputFile<<"#channel, timeOffset"<<std::endl;
    for(auto iter=timeOffsets.begin(); iter!=timeOffsets.end(); ++iter)
    {
      outputFile<<iter->first<<","<<iter->second<<std::endl;
    }

    outputFile.close();
    inputFile->Close();
}

bool FindSPEpeak(TH1F *hist, TSpectrum &spectrum, double &SPEpeak)
{
    if(hist->GetEntries()<minHistEntries) return false; //not enough data

    int nPeaks = spectrum.Search(hist,spectrumPeakSigma,"nodraw",spectrumPeakThreshold);
    if(nPeaks==0) return false;

    //peaks are not returned sorted
    double *peaksX = spectrum.GetPositionX();
    double *peaksY = spectrum.GetPositionY();
    std::vector<std::pair<double,double> > peaks;
    for(int iPeak=0; iPeak<nPeaks; ++iPeak) peaks.emplace_back(peaksX[iPeak],peaksY[iPeak]);
    std::sort(peaks.begin(),peaks.end(), [](const std::pair<double,double> &a, const std::pair<double,double> &b) {return a.first<b.first;});

    int peakToUse=0;
    if(nPeaks>1 && peaks[0].first>0)   //if more than one peak is found, the first peak could be due to baseline fluctuations
    {
      if(fabs(peaks[1].first/peaks[0].first-2.0)>peakRatioTolerance) peakToUse=1; //2nd peak is not twice the 1st peak, so the 1st peak is not the SPE peak
                                                                                  //assume that the 2nd peak is the SPE peak
                                                                                  //we have never seen that the 3rd peak was the SPE peak - no need to test it
    }
    SPEpeak = peaks[peakToUse].first;
    return true;
}
