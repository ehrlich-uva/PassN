void Plot(TFile *file, TVirtualPad *c, int pad, const std::string &pdfFileName, int feb1, int fpga1, int feb2, int fpga2, std::map<std::pair<int,int>,float> &measuredTimeDiffs)
{
  if(feb1==1) feb1=28;
  if(feb2==1) feb2=28;
  TH1F *plot = (TH1F*)file->FindObjectAny(Form("fpgaTimeDiff_%i_%i",feb1*4+fpga1,feb2*4+fpga2));
  if(plot==NULL) return;

  c->cd(pad);
  plot->SetTitle(Form("FEB/FPGA %i/%i, FEB/FPGA %i/%i",feb1+1,fpga1,feb2+1,fpga2));
  plot->Draw();

  float mean=plot->GetMean();
  measuredTimeDiffs[std::pair<int,int>(feb1*4+fpga1,feb2*4+fpga2)]=mean;
}

void CrvTimeOffsets_extracted(const std::string &rootFileName, const std::string &calibFileName, const std::string &pdfFileName)
{
  std::string pdfFileName="timeDifferences.pdf";

  std::map<std::pair<int,int>,float> measuredTimeDiffs;
  std::map<int,float> timeOffsets;
  timeOffsets[(1-1)*4]=0;      //FEB1 is used as reference
  timeOffsets[(17-1)*4]=-112;  //FEB17 has 72ft longer cable (for testing purpose)
  timeOffsets[(25-1)*4]=-7.8;  //FEB25 has 5ft longer cable

  TCanvas c0;
  c0.Print(Form("%s[", pdfFileName.c_str()), "pdf");

  TFile *file = new TFile(rootFileName.c_str());
  file->cd("CrvTimingStudies");

  //compare all FPGAs of an FEB
  for(int feb=0; feb<28; ++feb)
  {
    TCanvas *c = new TCanvas(Form("FEB %i",feb+1),Form("FEB %i",feb+1),800,800);
    TText *t = new TText(.1,.8,Form("FEB %i",(feb!=1?feb+1:29)));

    c->Divide(1,2);
    c->cd(1);
    t->SetTextSize(0.15);
    t->Draw();

    c->cd(2);
    TVirtualPad *pad=gPad;
    pad->Divide(3,1);

    Plot(file,pad,1,pdfFileName, feb,0, feb,2, measuredTimeDiffs);
    Plot(file,pad,2,pdfFileName, feb,1, feb,3, measuredTimeDiffs);
    Plot(file,pad,3,pdfFileName, feb,0, feb,3, measuredTimeDiffs);

    c->Print(pdfFileName.c_str(),"pdf");
  }

  //compare both FEBs at a module's readout side
  for(int crvmodule=0; crvmodule<8; ++crvmodule)
  {
    TCanvas *c = new TCanvas(Form("Module %i",crvmodule+1),Form("Module %i",crvmodule+1),800,800);
    TText *t = new TText(.1,.8,Form("Module %i (FEBs on same side)",crvmodule+1));

    c->Divide(1,2);
    c->cd(1);
    t->SetTextSize(0.1);
    t->Draw();

    c->cd(2);
    TVirtualPad *pad=gPad;
    pad->Divide(3,1);

    if(crvmodule<6)
    {
      Plot(file,pad,1,pdfFileName, crvmodule*4+0,3, crvmodule*4+1,1, measuredTimeDiffs);
      Plot(file,pad,2,pdfFileName, crvmodule*4+2,3, crvmodule*4+3,1, measuredTimeDiffs);
    }
    else if(crvmodule==6) Plot(file,pad,1,pdfFileName, 24,3, 25,1, measuredTimeDiffs);
    else if(crvmodule==7) Plot(file,pad,1,pdfFileName, 26,3, 27,1, measuredTimeDiffs);

    c->Print(pdfFileName.c_str(),"pdf");
  }

  //compare neighboring modules
  for(int crvmodule=0; crvmodule<8; ++crvmodule)
  {
    if(crvmodule==3) continue; //last Tmodule
    if(crvmodule==5) continue; //last LEmodule
    if(crvmodule==7) continue; //last DSmodule
    TCanvas *c = new TCanvas(Form("Module %i / Module %i",crvmodule+1,crvmodule+2),Form("Module %i / Module %i",crvmodule+1,crvmodule+2),800,800);
    TText *t = new TText(.1,.8,Form("Module %i / Module %i",crvmodule+1,crvmodule+2));

    c->Divide(1,2);
    c->cd(1);
    t->SetTextSize(0.15);
    t->Draw();

    c->cd(2);
    TVirtualPad *pad=gPad;
    pad->Divide(3,1);

    if(crvmodule<6)
    {
      Plot(file,pad,1,pdfFileName, crvmodule*4+0,1, crvmodule*4+5,2, measuredTimeDiffs);
      Plot(file,pad,2,pdfFileName, crvmodule*4+2,1, crvmodule*4+7,2, measuredTimeDiffs);
    }
    else if(crvmodule==6) Plot(file,pad,1,pdfFileName, 24,1, 27,2, measuredTimeDiffs);

    c->Print(pdfFileName.c_str(),"pdf");
  }

  //compare FEBs of opposite readout sides
  {
    TCanvas *c = new TCanvas("Modules opposite sides","Modules opposite sides",800,800);
    TText *t1 = new TText(.1,.85,"Module 1 (FEBs on opposide sides)");
    TText *t2 = new TText(.1,.35,"Module 5 (FEBs on opposide sides)");

    c->Divide(1,4);
    c->cd(1);
    t1->SetTextSize(0.15);
    t1->Draw();
    c->cd(3);
    t2->SetTextSize(0.15);
    t2->Draw();

    c->cd(2);
    TVirtualPad *pad1=gPad;
    pad1->Divide(3,1);

    c->cd(4);
    TVirtualPad *pad2=gPad;
    pad2->Divide(3,1);

    Plot(file,pad1,1,pdfFileName, 0*4+0,0, 0*4+2,0, measuredTimeDiffs);
    Plot(file,pad2,1,pdfFileName, 4*4+0,1, 4*4+2,2, measuredTimeDiffs);

    c->Print(pdfFileName.c_str(),"pdf");
  }

  //missing: between different sectors

  c0.Print(Form("%s]", pdfFileName.c_str()), "pdf");

  /************************************************************************/

  size_t prevSize=measuredTimeDiffs.size();
  for(auto measuredTimeDiff=measuredTimeDiffs.cbegin(); ; )
  {
    if(measuredTimeDiff==measuredTimeDiffs.cend())
    {
      if(measuredTimeDiffs.size()==0) break;  //all measured time differences used.
      if(prevSize==measuredTimeDiffs.size()) break;  //loops don't seem to make a difference anymore.

      prevSize=measuredTimeDiffs.size();
      measuredTimeDiff=measuredTimeDiffs.cbegin();  //there are still some measured time differences left. start the loop again.
    }

    //check if this particular measured time difference can be connected to a point where the time offset is known already.
    bool erased=false;
    for(const auto& timeOffset: timeOffsets)
    {
      //Longer cables result in an earlier reco time, because the start t0 arrives later at the FEBs, so that the time difference between t0 and tHit is shorter (=tReco).
      //Therefore, the dt plots result in a negative mean value.
      //The time offset needs to be positive to counter act it. That's why a negative sign is used below.
      if(measuredTimeDiff->first.first==timeOffset.first)
      {
        timeOffsets[measuredTimeDiff->first.second]=timeOffset.second-measuredTimeDiff->second;
        measuredTimeDiff=measuredTimeDiffs.erase(measuredTimeDiff);
	erased=true;
	break;
      }
      if(measuredTimeDiff->first.second==timeOffset.first)
      {
        timeOffsets[measuredTimeDiff->first.first]=timeOffset.second+measuredTimeDiff->second;
        measuredTimeDiff=measuredTimeDiffs.erase(measuredTimeDiff);
	erased=true;
	break;
      }
    }
    if(erased) continue;

    ++measuredTimeDiff;
  }

  if(measuredTimeDiffs.size()>0) std::cout<<"There are still some unused measured time diffs!"<<std::endl;
  for(const auto& measuredTimeDiff: measuredTimeDiffs)
  {
    std::cout<<measuredTimeDiff.first.first<<"/"<<measuredTimeDiff.first.second<<" "<<measuredTimeDiff.second<<std::endl;
  }
  for(const auto& timeOffset: timeOffsets)
  {
    int fpga=timeOffset.first%4;
    int globalfeb=timeOffset.first/4;
    int feb=globalfeb%24+1;
    int roc=globalfeb/24+1;
    std::cout<<"ROC "<<roc<<"    FEB "<<feb<<"   fpga "<<fpga<<"     timeOffset "<<timeOffset.second<<std::endl;
  }

  /************************************************************************/

  TTree *treeChannelMap = (TTree*)file->FindObjectAny("channelMap");
  TTree *treeCalib = (TTree*)file->FindObjectAny("crvCalib");
  if(treeChannelMap==NULL || treeCalib==NULL)
  {
    std::cout<<"Couldn't find channel map or calibration constants!"<<std::endl;
    return;
  }

  std::ofstream calibFile(calibFileName);

  size_t channel;
  double pedestal, calibPulseHeight, calibPulseArea;
  treeCalib->SetBranchAddress("channel", &channel);
  treeCalib->SetBranchAddress("pedestal", &pedestal);
  treeCalib->SetBranchAddress("calibPulseHeight", &calibPulseHeight);
  treeCalib->SetBranchAddress("calibPulseArea", &calibPulseArea);

  calibFile<<"TABLE CRVSiPM"<<std::endl;
  calibFile<<"#channel,pedestal,calibPulseHeight,calibPulseArea"<<std::endl;
  for(int i=0; i<treeCalib->GetEntries(); ++i)
  {
    treeCalib->GetEntry(i);
    calibFile<<channel<<","<<pedestal<<","<<calibPulseHeight<<","<<calibPulseArea<<std::endl;
  }

  int    roc, feb, febChannel;
  std::map<int,float> timeOffsetChannelMap;
  treeChannelMap->SetBranchAddress("channel", &channel);
  treeChannelMap->SetBranchAddress("roc", &roc);
  treeChannelMap->SetBranchAddress("feb", &feb);
  treeChannelMap->SetBranchAddress("febChannel", &febChannel);
  for(int i=0; i<treeChannelMap->GetEntries(); ++i)
  {
    treeChannelMap->GetEntry(i);
    int globalfeb = (roc-1)*24 + (feb-1);
    int globalfpga = globalfeb*4 + febChannel/16;
    timeOffsetChannelMap[channel]=timeOffsets.at(globalfpga);
  }

  calibFile<<"TABLE CRVTime"<<std::endl;
  calibFile<<"#channel,timeOffset"<<std::endl;
  for(int i=0; i<treeCalib->GetEntries(); ++i) //require the same channels for the SiPM table and Time table
                                               //incl. for non-existing channels to get a continuously running index
  {
    treeCalib->GetEntry(i); //just to get the channel number
    float timeOffset=0;
    if(timeOffsetChannelMap.find(channel)!=timeOffsetChannelMap.end()) timeOffset=timeOffsetChannelMap.at(channel);
    calibFile<<channel<<","<<timeOffset<<std::endl;
  }

  calibFile.close();
}
