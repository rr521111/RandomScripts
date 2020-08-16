TString promptdir = "/u/group/halla/parity/software/japan_offline/prompt/prex-prompt";
TString rrwork = "/w/hallc-scifs17exp/qweak/rradloff/crex-runlist/prex-runlist";

TString pcrex = "PREX";
vector<TString> types = {"arm_flag", "beam_current", "beam_energy", "bmw", "components", "component_stats", "event_count", "event_rate", "experiment", "feedback", "FFB", "flip_state", "good_charge", " helicity_frequency", "helicity_pattern", "horizontal_wien", "ihwp", "is_valid_run_end", "prompt_analysis", "respin_comment", "rhwp", "rtvs", "run_config", "run_end_time", "run_flag", "run_length", "run_prestart_time", "run_start_epoch", "run_start_time", "run_type", "session", "slug", "target_45encoder", "target_90encoder", "target_encoder", "target_type", "total_charge", "user_comment", "vertical_wien", "wac_comment", "run_number"};
vector<TString> types_line = {" ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " "};

vector<Double_t> OpenRun(Int_t runnum, TString ihwps, TString wiens);

vector<Double_t> GetSlugVals(Int_t slugnum, vector<Double_t> runs, vector<Double_t> slugs, vector<Double_t> arms, vector<Double_t> means, vector<Double_t> errors);

void ReadingFiles(){
  
  vector<vector<TString>> type_table;

  ifstream infile;
  infile.open(Form("./pcrex_run_data.list", promptdir.Data()));
  
  string ins;
  vector<Double_t> Vals;
  vector<Double_t> RunMeans;
  vector<Double_t> RunErrors;
  vector<Double_t> Zeros;

  Int_t i = 0;
  Int_t First = 5000;
  Int_t Last = 5010;
  while(infile >> ins){
    stringstream line(ins);
    types_line.clear();
    while(line.good()){
      string substr;
      getline(line, substr, ',');
      types_line.push_back(substr);
    }
    
    TString arm = types_line[0];
    TString wien = types_line[11];
    TString ihwp = types_line[16];
    TString goodbad = types_line[24];
    TString production = types_line[29];
    Int_t slug = stoi(types_line[31]);
    TString target = types_line[35];
    Int_t run = stoi(types_line[40]);

    if(run >= First && run <= Last && slug >= 4000 && slug <= 4500 && production == "\'Production\'" && goodbad == "\'Good\'"){
      Vals = OpenRun(run, ihwp, wien);
      if(Vals[0] == -1){
        continue;
      }

      type_table.push_back(types_line);

      RunMeans.push_back(Vals[0]);
      RunErrors.push_back(Vals[1]);
      Zeros.push_back(0);

      i++;
    }
  }
  
  infile.close();
  
  Double_t RunsArr[type_table.size()];
  Double_t RunMeansArr[type_table.size()];
  Double_t ZerosArr[type_table.size()];
  Double_t RunErrorsArr[type_table.size()];
  for(Int_t j=0; j<i; j++){
    RunsArr[j] = stoi(type_table[j][40]);
    ZerosArr[j] = 0;
    RunMeansArr[j] = RunMeans[j];
    RunErrorsArr[j] = RunErrors[j];
  }
  
  TCanvas *c1 = new TCanvas();
  TGraphErrors *plotallruns = new TGraphErrors(type_table.size(), RunsArr, RunMeansArr, ZerosArr, RunErrorsArr);
  plotallruns->SetTitle("PREX Runs Upstream DD");
  plotallruns->SetMarkerStyle(7);
  plotallruns->Draw("AP");

  vector<Double_t> SlugNumbers;
  vector<Double_t> SlugMeans;
  vector<Double_t> SlugErrors;
  vector<Double_t> SlugZeros;
  vector<Double_t> slugouts;
  for(Int_t p = 0; p<5000; p++){
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

