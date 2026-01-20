const double fitRangeStart=0.8;
const double fitRangeEnd=1.2;
const int    minHistEntries=100;
const double minPeakPulseHeight=10.0;
const double minPeakPulseArea=250.0;
const int    spectrumNPeaks=100;
const double spectrumPeakSigma=4.0;
//const double spectrumPeakThreshold=0.001;
const double spectrumPeakThreshold=0.01;
//const double peakRatioTolerance=0.3;

//template<size_t N>
//bool FindSPEpeak(TH1F *hist, TSpectrum &spectrum, std::array<TF1,N> &functions, double &SPEpeak, double minPeak);
bool FindSPEpeak(TH1F *hist, TSpectrum &spectrum, TF1 &function, double &SPEpeak, double minPeak);

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

//    std::array functions={TF1("calibPeak1","gaus"), TF1("calibPeak2","gaus"), TF1("calibPeak3","gaus")}; //only need to fit three peaks
    TF1 function("calibPeak","gaus");
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
      for(int i=0; i<2; ++i) //loop over histograms with pulse areas and pulse heights
      {
        if(i==1) hist=(TH1F*)gDirectory->FindObjectAny(Form("crvCalibrationHistPulseArea_%zu",channel));
        else hist=(TH1F*)gDirectory->FindObjectAny(Form("crvCalibrationHistPulseHeight_%zu",channel));
        hist->GetListOfFunctions()->Delete();

        double SPEpeak=-1;
//        if(!FindSPEpeak(hist, spectrum, functions, SPEpeak, (i==0?minPeakPulseHeight:minPeakPulseArea)))
        if(!FindSPEpeak(hist, spectrum, function, SPEpeak, (i==0?minPeakPulseHeight:minPeakPulseArea)))
        {
          calibValue[i]=-1;
          continue;
        }
        calibValue[i]=SPEpeak;
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

bool FindSPEpeak(TH1F *hist, TSpectrum &spectrum, TF1 &function, double &SPEpeak, double minPeak)
{
    if(hist->GetEntries()<minHistEntries) return false; //not enough data

    int nPeaks = spectrum.Search(hist,spectrumPeakSigma,"nodraw",spectrumPeakThreshold);
    if(nPeaks<=0) return false;

    //peaks are returned sorted by Y
    //from our long-time experience:
    //-SPE peak is either the highest peak or second highest peak
    //-if the SPE peak is the second highest, then the highest peak comes from the baseline
    //-the peak from the baseline is always below the minPeak threshold, while the SPE peak is not
    //-the minPeak threshold may have to be adjusted for non-standard bias voltages
    double *peaksX = spectrum.GetPositionX();
    double x=peaksX[0];
    if(x<minPeak)
    {
      if(nPeaks==1) return false;
      x=peaksX[1];
      if(x<minPeak) return false;
    }

    if(hist->FindBin(x*fitRangeStart)==hist->FindBin(x*fitRangeEnd)) return false; //fit range start/end are in the same bin
    function.SetRange(x*fitRangeStart,x*fitRangeEnd);
    function.SetParameter(1,x);
    hist->Fit(&function, "QR");
    SPEpeak = function.GetParameter(1);

    return true;
}

/*
template<size_t N>
bool FindSPEpeak(TH1F *hist, TSpectrum &spectrum, std::array<TF1,N> &functions, double &SPEpeak, double minPeak)
{
    if(hist->GetEntries()<minHistEntries) return false; //not enough data

    size_t nPeaks = spectrum.Search(hist,spectrumPeakSigma,"nodraw",spectrumPeakThreshold);
    if(nPeaks==0) return false;

    //peaks are returned sorted by Y
    double *peaksX = spectrum.GetPositionX();
    std::vector<double> fittedPeaks;
    for(size_t iPeak=0; iPeak<nPeaks && iPeak<functions.size(); ++iPeak)
    {
      double x=peaksX[iPeak];
      if(hist->FindBin(x*fitRangeStart)==hist->FindBin(x*fitRangeEnd)) continue; //fit range start/end are in the same bin
      functions[iPeak].SetRange(x*fitRangeStart,x*fitRangeEnd);
      functions[iPeak].SetParameter(1,x);
      hist->Fit(&functions[iPeak], "QR+");
      fittedPeaks.emplace_back(functions[iPeak].GetParameter(1));
    }
    if(fittedPeaks.size()==0) return false;

    int peakToUse=-1;
    //only need to test two highest peaks (=first two entries in vector)
    //one of the peaks could be due to baseline fluctuations
    for(size_t iPeak=0; iPeak<fittedPeaks.size() && iPeak<2; ++iPeak)
    for(size_t jPeak=iPeak+1; jPeak<fittedPeaks.size(); ++jPeak)
    {
      if(fabs(fittedPeaks.at(jPeak)/fittedPeaks.at(iPeak)-2.0)<peakRatioTolerance) peakToUse=iPeak;
    }
    if(peakToUse==-1) return false;
    if(fittedPeaks.at(peakToUse)<minPeak) return false;

    SPEpeak = fittedPeaks.at(peakToUse);

    return true;
}
*/
