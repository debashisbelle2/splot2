//Make the splot
//This is a simple example showing how to make an sPlot using RooStats and RooFit.
//The method is described in the paper paper by Pvik and De Liberder, 
// Nucl. Instrum. Meth. A555 356 (2005) [arXiv:physics/0402083]
//Here the data consists of signal events and background events. The method describes how to separate signal events from background.
//So two variables are used , one is discriminating variable and another one is the control varible. Here discriminating variable is jpsi mass and the control variable is the muon liklihood ratio. At first, We generate a two component data set from two varibles stored in the ntuple. Then fit the jpsi mass with two pdfs. Then calculate the signal and background weights to plot MuId splots. Two splots are found with signal contributions and backgrounds.

//Author: Debashis Sahoo, Utkal University, Bhubaneswar & TIFR, Mumbai & KEK,Japan


void splot(TString fname="allmc_stream0.root") {
using namespace RooFit;
using namespace RooStats;
gSystem->Load("libRooFit");
gSystem->Load("libRooStats.so");
RooRealVar mass("mass", "GeV",2.9,3.3);
RooRealVar muid("muid", "MuId",0,1);
RooDataSet* data = new RooDataSet("data", "data",RooArgSet(mass,muid));

// fit variable
Float_t jpsimass,mumid_kl;

// Read data file
TFile* input=new TFile(fname);
TTree* t1 = (TTree*) input->Get("h1");

// no of entries in the datafile

Int_t n_tot = (int)t1->GetEntries();
cout << " n_tot = " << n_tot << endl;
t1->SetBranchAddress("jpsimass",&jpsimass);
t1->SetBranchAddress("mum_l",&mumid_kl);
// Import the data points
  for(int i=0; i<n_tot; i++) {
      t1->GetEntry(i); 
  if(jpsimass >=2.9 && jpsimass <=3.3) {   
     mass.setVal(jpsimass); 
     muid.setVal(mumid_kl);
     data->add(RooArgSet(mass,muid));       
     }                                                     
  }

  
 cout << "make mass model" << endl;
 //make mass model for signal(jpsi) .............add jp infront of each var
 RooRealVar jpMmean("jpM#mu_{x}","jpMmean",3.09,3.05,3.12);
 RooRealVar jpMsigma1("jpM#sigma_{1}","jpmssigma1",0.006,0.0,0.0075);
 RooGaussian jpMsignal1("jpMsignal1","jpsignal1",mass,jpMmean,jpMsigma1);
 RooRealVar jpMsigma2("jpM#sigma_{2}","jpsigma2",0.03,0.0,0.08);
 RooGaussian jpMsignal2("jpMsignal2","jpsignal2",mass,jpMmean,jpMsigma2);
 RooRealVar jpMsfrac("jpMsfrac","Area fraction",0.5,0.0 ,1.0);
 RooAddPdf jpMsignald("jpMsignald","jpMsignald",RooArgList(jpMsignal1,jpMsignal2),jpMsfrac);

//make mass model for background ............  
 RooRealVar bkMa("bkMa","bkMa",-2.0,-4.0,0);
 RooExponential bkMbkg("bkMbkg","bkMexp",mass,bkMa);

 RooRealVar nsig("nsig", "nsig",100000,50000,900000);
 RooRealVar nbkg("nbkg", "nbkg",100000,50000,950000);

 RooAddPdf massPdf("massPdf", "", RooArgList(jpMsignald, bkMbkg), RooArgList(nsig, nbkg));

  // Fit the data for jpsi mass in order to make splots of muid

  RooFitResult * r = (RooFitResult*)massPdf.fitTo(*data,"etrm", RooFit::Save());
  r->Print("v");

  
  // Create the sPlot data from the fit to mass.  The sweights computed here are valid for plotting pt
  RooStats::SPlot* sDataX = new RooStats::SPlot("sData","An SPlot", *data, &massPdf, RooArgList(nsig, nbkg) );
  cout << "Yield of signal is " << nsig.getVal() << " From sWeights it is " << sDataX->GetYieldFromSWeight("nsig") <<endl;
  cout << "Yield of background is " << nbkg.getVal() << "  From sWeights it is " << sDataX->GetYieldFromSWeight("nbkg") <<endl;

 // cout << "signal Weight   " << sDataX->GetSWeight("nsig") 
 //      << " bkg Weight   " <<  sDataX->GetSWeight("nbkg") <<endl;
   
  for(Int_t i=0; i < 100; i++)
    {
      cout << "signal Weight   " << sDataX->GetSWeight(i,"nsig") 
		<< " bkg Weight   " << sDataX->GetSWeight(i,"nbkg") 
		<< "  Total Weight   " << sDataX->GetSumOfEventSWeight(i) 
		<< endl;
    }
  
  // create weighted data sets
  RooDataSet * dataw_sig = new RooDataSet(data->GetName(),data->GetTitle(),data,*data->get(),0,"nsig_sw") ;
  RooDataSet * dataw_bg  = new RooDataSet(data->GetName(),data->GetTitle(),data,*data->get(),0,"nbkg_sw") ;

  cout << "Making splots of the signal and background" <<endl;
  RooPlot* frame_sig_muid = muid.frame(100) ;
  frame_sig_muid->SetTitle("sPlot for the signal MuId distribution");
  dataw_sig->plotOn(frame_sig_muid, RooFit::DataError(RooAbsData::SumW2) ) ;

  RooPlot* frame_bg_muid = muid.frame(100) ;
  frame_bg_muid->SetTitle("sPlot for the background MuId distribtuion");
  dataw_bg->plotOn(frame_bg_muid, RooFit::DataError(RooAbsData::SumW2) ) ;

  // plot PDFs onto the sPlots; if the sPlot correctly reweights the data then the signal / background only shapes should
  // provide an adequate description of the data.
 // signaldjpt.plotOn(frame_sig_pt);
 // signal2bpt.plotOn(frame_bg_pt);

  TCanvas *c1 = new TCanvas("c1" ,"c1");
  c1->cd();
  frame_sig_muid->Draw();
 
  TCanvas *c2 = new TCanvas("c2","c2");
  c2->cd();
  frame_bg_muid->Draw();
  
 // c1.Print("splot_pt.pdf");
 c1->SaveAs("splot_muid_sig.C");
 c2->SaveAs("splot_muid_bkg.C");
  }
















