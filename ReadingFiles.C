TString promptdir = "/u/group/halla/parity/software/japan_offline/prompt/prex-prompt";
TString rrwork = "/w/hallc-scifs17exp/qweak/rradloff/crex-runlist/prex-runlist";

TString pcrex = "PREX";

vector<Double_t> OpenRun(Int_t runnum, TString ihwps, TString wiens);

vector<Double_t> GetSlugVals(Int_t slugnum, vector<Double_t> runs, vector<Double_t> slugs, vector<Double_t> arms, vector<Double_t> means, vector<Double_t> errors);

void ReadingFiles(){

  vector<Double_t> Runs;
  vector<Double_t> RunMeans;
  vector<Double_t> RunErrors;
  vector<Double_t> Slugs;
  vector<TString> ALL_PREXs;
  vector<TString> goodbads;
  vector<TString> IHWPs;
  vector<TString> Wiens;
  vector<Double_t> Arms;
  vector<Double_t> Zeros;

  ifstream infile;
  if(pcrex == "PREX"){
    infile.open(Form("%s/prex-runlist/all_production.list", promptdir.Data()));
  }else{
    infile.open(Form("%s/all_production_crex.list", rrwork.Data()));
  }
  
  string ins;
  vector<vector<string>> result;
  vector<string> subresult;
  vector<Double_t> Vals;
  Int_t i = 0;
  Int_t First = 0;
  Int_t Last = 5000;
  while(infile >> ins){
    stringstream line(ins);
    subresult.clear();
    while(line.good()){
      string substr;
      getline(line, substr, ',');
      subresult.push_back(substr);
    }
    
    if(stoi(subresult[0]) >= First && stoi(subresult[0]) <= Last){
      Vals = OpenRun(stoi(subresult[0]), subresult[4], subresult[5]);
      if(Vals[0] == -1){
        continue;
      }

      result.push_back(subresult);
      Runs.push_back(stoi(result[i][0]));
      Slugs.push_back(stoi(result[i][1]));
      ALL_PREXs.push_back(result[i][2]);
      goodbads.push_back(result[i][3]);
      IHWPs.push_back(result[i][4]);
      Wiens.push_back(result[i][5]);
      Arms.push_back(stoi(result[i][6]));

      RunMeans.push_back(Vals[0]);
      RunErrors.push_back(Vals[1]);
      Zeros.push_back(0);

      i++;
    }
  }
  
  infile.close();
  
  Double_t RunsArr[Runs.size()];
  Double_t RunMeansArr[Runs.size()];
  Double_t ZerosArr[Runs.size()];
  Double_t RunErrorsArr[Runs.size()];
  for(Int_t j=0; j<i; j++){
    RunsArr[j] = Runs[j];
    ZerosArr[j] = 0;
    RunMeansArr[j] = RunMeans[j];
    RunErrorsArr[j] = RunErrors[j];
  }
  
  TCanvas *c1 = new TCanvas();
  TGraphErrors *plotallruns = new TGraphErrors(Runs.size(), RunsArr, RunMeansArr, ZerosArr, RunErrorsArr);
  plotallruns->SetTitle("PREX Runs Upstream DD");
  plotallruns->SetMarkerStyle(7);
  plotallruns->Draw("AP");

  vector<Double_t> SlugNumbers;
  vector<Double_t> SlugMeans;
  vector<Double_t> SlugErrors;
  vector<Double_t> SlugZeros;
  vector<Double_t> slugouts;
  for(Int_t p = 0; p<200; p++){
    slugouts = GetSlugVals(p, Runs, Slugs, Arms, RunMeans, RunErrors);
    if(slugouts[0] == -1){
      continue;
    }
    SlugNumbers.push_back(p);
    SlugMeans.push_back(slugouts[0]);
    SlugErrors.push_back(slugouts[1]);
    SlugZeros.push_back(0);
  }

  Double_t SlugNumbersArr[SlugNumbers.size()];
  Double_t SlugMeansArr[SlugNumbers.size()];
  Double_t SlugErrorsArr[SlugNumbers.size()];
  Double_t SlugZerosArr[SlugNumbers.size()];
  for(Int_t j=0; j<SlugNumbers.size(); j++){
    SlugNumbersArr[j] = SlugNumbers[j];
    SlugZerosArr[j] = 0;
    SlugMeansArr[j] = SlugMeans[j];
    SlugErrorsArr[j] = SlugErrors[j];
  }

  TCanvas *c2 = new TCanvas();
  TGraphErrors *plotallslugs = new TGraphErrors(SlugNumbers.size(), SlugNumbersArr, SlugMeansArr, SlugZerosArr, SlugErrorsArr);
  plotallslugs->SetTitle("PREX Slugs Upstream DD");
  plotallslugs->SetMarkerStyle(7);
  plotallslugs->Draw("AP");
 
  return;
}

vector<Double_t> OpenRun(Int_t runnum, TString ihwps, TString wiens){

  vector<Double_t> vals;
  TFile *runfile = TFile::Open(Form("%s/japanOutput/prexPrompt_pass2_%d.000.root", promptdir.Data(), runnum), "READ");
  if(runfile==NULL){
    vals = {-1, -1};
    return vals;
  }

  //TTree *evt_tree = (TTree*)gROOT->FindObject("evt");
  TTree *mul_tree = (TTree*)gROOT->FindObject("mul");
  TTree *mulc_lrb_alldet_tree = (TTree*)gROOT->FindObject("mulc_lrb_alldet");

  mulc_lrb_alldet_tree->AddFriend("mul");

  cout << "inside " << runnum << endl;
  Int_t test = mulc_lrb_alldet_tree->Draw("cor_asym_us_dd:Entry$>>htemp()","mul.ErrorFlag==0","goff");
  if(test==-1){
    vals = {-1, -1};
    return vals;
  }
  
  Int_t ihwp = 1;
  Int_t wien = 1;
  if(ihwps == "IN"){
    ihwp = -1;
  }
  if(wiens == "FLIP-RIGHT"){
    wien = -1;
  }

  TH2F *hout = (TH2F*)gROOT->FindObject("htemp");
  vals = {ihwp*wien*hout->GetMean(2), hout->GetRMS(2)};

  runfile->Close();

  return vals;
}

vector<Double_t> GetSlugVals(Int_t slugnum, vector<Double_t> runs, vector<Double_t> slugs, vector<Double_t> arms, vector<Double_t> means, vector<Double_t> errors){
  vector<Double_t> slugmeanerr = {-1,-1};

  vector<Double_t> slugruns;
  vector<Double_t> slugmeans;
  vector<Double_t> slugerrors;
  vector<Double_t> slugzeros;

  for(Int_t i = 0; i < runs.size(); i++){
    if(slugs[i] == slugnum){
      if(arms[i]!=0){
        return slugmeanerr;
      }
      slugruns.push_back(runs[i]);
      slugmeans.push_back(means[i]);
      slugerrors.push_back(errors[i]);
      slugzeros.push_back(0);
    }
  }

  if(slugruns.size()==0){
    return slugmeanerr; 
  }

  Double_t slugrunsarr[slugruns.size()];
  Double_t slugerrorsarr[slugruns.size()];
  Double_t slugmeansarr[slugruns.size()];
  Double_t slugzerosarr[slugruns.size()];

  for(Int_t j = 0; j < slugruns.size(); j++){
    slugrunsarr[j] = slugruns[j];
    slugerrorsarr[j] = slugerrors[j];
    slugmeansarr[j] = slugmeans[j];
    slugzerosarr[j] = 0;
  }

  TGraphErrors *plot = new TGraphErrors(slugruns.size(), slugrunsarr, slugmeansarr, slugzerosarr, slugerrorsarr);

  slugmeanerr = {plot->GetMean(2), plot->GetRMS(2)};

  return slugmeanerr;  
}

