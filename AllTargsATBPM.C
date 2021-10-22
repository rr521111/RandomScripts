//##This script is deigned to loop through a runlist and calculate means and errors of values in those runs. Currently setup to find the asymmetry across all runs.

TString promptdir = "/u/group/halla/parity/software/japan_offline/prompt/prex-prompt";
TString rootfiles = "$QW_ROOTFILES";
TString rootfiles2 = "/volatile/halla/parity/crex-respin2/japanOutput";

//goal settings
TString Target = "Carbon";
//TString Tree = "burst_mulc_lrb_alldet";
TString Branch = "cor_asym_us_dd";
TString ValueLeaf = "hw_sum";
TString ErrorLeaf = "hw_sum_err";
TString VarCombo = "diff_bpm4aX";
Int_t xpos = 0;

//mostly used for reference; all of the rcdb data types in order, with an added entry for run number at the end.
vector<TString> types = { "arm_flag", "beam_current", "beam_energy", "bmw", "components", "component_stats", "event_count", "event_rate", "experiment", "feedback", "FFB", "flip_state", "good_charge", " helicity_frequency", "helicity_pattern", "horizontal_wien", "ihwp", "is_valid_run_end", "prompt_analysis", "respin_comment", "rhwp", "rtvs", "run_config", "run_end_time", "run_flag", "run_length", "run_prestart_time", "run_start_epoch", "run_start_time", "run_type", "session", "slug", "target_45encoder", "target_90encoder", "target_encoder", "target_type", "total_charge", "user_comment", "vertical_wien", "wac_comment", "run_number" };
vector<TString> types_line = { " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " " };

void oldmain(TString targ, TString pcrex, TString rungroup);

//opens the rootfile for a a given run, sign corrects, then returns the mean and error.
vector<vector<Double_t>> OpenRun(TString directory, Int_t runnum, Int_t slugnum, TString ihwps, TString wiens, TString arm);

//averages the means for miniruns in a given slug, returns the slug mean and error.
vector<Double_t> GetSlugVals(Int_t slugnum, vector<Double_t> miniruns, vector<Int_t> slugs, vector<TString> arms, vector<Double_t> means, vector<Double_t> errors);

//takes a vector of vectors and returns a tgrapherrors. (lol why is this not built in?)
TGraphErrors* PlotFromMatrix(vector<vector<Double_t>> data, TString title);

//takes a vector of vectors and returns a tgrapherrors of the pull.
TH1D* PullFromMatrix(vector<vector<Double_t>> data, TString title);
  
//takes the matrix format mentioned above and plots the second column(means)
TH1D* HistFromMatrix(vector<vector<Double_t>> data, TString title);

void AllTargsATBPM(){

    //cout << "run, mini, ihwp, wien, targ, 0vdd, D0vdd, 1vdd, D1vdd, 2vdd, D2vdd, 3vdd, D3vdd, 4vdd, D4vdd, asym, Dasym, ex, Dex, ey, Dey, ax, Dax, ay, Day, E, DE, 1x, D1x, 12x, D12x" << endl;

    vector<vector<TString>> allinputs = {
        {"Pb",     "PREX", "1"},
        {"Pb",     "PREX", "2"},
        {"40",     "PREX", "1"},
        {"40",     "PREX", "2"},
        {"Carbon", "PREX", "1"},
        {"Carbon", "PREX", "2"}/*,
        {"Pb",     "CREX", "1"},
        {"40",     "CREX", "1"},
        {"40",     "CREX", "2"},
        {"48",     "CREX", "1"},
        {"48",     "CREX", "2"},
        {"48",     "CREX", "3"},
        {"48",     "CREX", "4"},
        {"Carbon", "CREX", "1"},
        {"Carbon", "CREX", "2"}*/
    };

    system("rm ./Outputs/temp/*");
    //system("rm ./Outputs/AllTargsAT.pdf");

    for(int i = 0; i<allinputs.size(); i++){
        //cout << allinputs[i][0] << ",  " << allinputs[i][1] << endl;
        oldmain(allinputs[i][0], allinputs[i][1], allinputs[i][2]);
    }

    //system("pdfunite ./Outputs/temp/* ./Outputs/AllTargsATexvdd.pdf");
    system("pdfunite ./Outputs/temp/* ./Outputs/AllTargsATBPM.pdf");

    return;
}

void oldmain(TString targ = "Pb", TString pcrex = "CREX", TString rungroup = "1") {
    
    Target = targ;
    Int_t slugmin;
    Int_t slugmax;
    TString prodstring;

    vector<vector<TString>> type_table;

    ifstream infile;
    infile.open("./another_runlist/pcrex_run_data.list");

    string ins;
    vector<Int_t> Runs;
    vector<Double_t> Miniruns;
    vector<Int_t> Slugs;
    vector<TString> Arms;
    vector<TString> Targets;
    vector<vector<Double_t>> Vals_table;
    vector<vector<Double_t>> Vals;
    vector<vector<Double_t>> Vals2;
    vector<Double_t> RunMeans;
    vector<Double_t> RunErrors;
    vector<vector<Double_t>> textmatrix;

    Int_t total = 0;
    Int_t i = xpos;
    Int_t First = 4106; // 3305
    //Int_t Last = 5420;
    Int_t Last =6408; // 4980

    if(pcrex.Contains("PREX")){
        slugmax = 600;
        slugmin = 400;
        prodstring = "A_T";
        rootfiles = "~/rrvolatile/rootfiles";
        if(targ=="Pb" && rungroup=="1"){
            First = 4110;
            Last = 4119;
        }else if(targ=="Pb" && rungroup=="2"){
            First = 4128;
            Last = 4130;
        }else if(targ=="Carbon" && rungroup=="1"){
            First = 4106;
            Last = 4109;
        }else if(targ=="Carbon" && rungroup=="2"){
            First = 4131;
            Last = 4133;
        }else if(targ=="40" && rungroup=="1"){
            First = 4120;
            Last = 4121;
        }else if(targ=="40" && rungroup=="2"){
            First = 4122;
            Last = 4127;
        }
    }else if(pcrex.Contains("CREX")){
        slugmax = 6000;
        slugmin = 4000;
        prodstring = "Production";
        rootfiles = "$QW_ROOTFILES";
        if(targ=="Pb" && rungroup=="1"){
            First = 6367;
            Last = 6379;
        }else if(targ=="Carbon" && rungroup=="1"){
            First = 6359;
            Last = 6366;
        }else if(targ=="Carbon" && rungroup=="2"){
            First = 6386;
            Last = 6393;
        }else if(targ=="40" && rungroup=="1"){
            First = 6349;
            Last = 6353;
        }else if(targ=="40" && rungroup=="2"){
            First = 6394;
            Last = 6404;
        }else if(targ=="48" && rungroup=="1"){
            First = 6344;
            Last = 6348;
        }else if(targ=="48" && rungroup=="2"){
            First = 6354;
            Last = 6358;
        }else if(targ=="48" && rungroup=="3"){
            First = 6380;
            Last = 6385;
        }else if(targ=="48" && rungroup=="4"){
            First = 6405;
            Last = 6408;
        }
    }

    while (!infile.eof()) {
        //splits the file at line breaks
        getline(infile, ins, '\n');
        stringstream line(ins);
        types_line.clear();
        //splits the line at commas
        while (line.good()) {
            string substr;
            getline(line, substr, ',');
            types_line.push_back(substr);

        }

        //avoids slug and run numbers that somehow aren't integer. I dont think any more of these exist after fixing the runlist.
        if (types_line[31] == " " || types_line[31] == "BEGONE_COMMAS" || types_line[40] == " " || types_line[40] == "BEGONE_COMMAS") {
            continue;
        }

        TString arm = types_line[0];
        TString wien = types_line[11];
        TString ihwp = types_line[16];
        TString goodbad = types_line[24];
        TString production = types_line[29];
        Int_t slug = stoi(types_line[31].Data());
        TString target = types_line[35];
        Int_t run = stoi(types_line[40].Data());

        //first filter to decide which runs to keep. usually production and good are safe bets.
        if (run >= First && run <= Last && run != 3140 && slug >= slugmin && slug <= slugmax && production.Contains(prodstring) && goodbad.Contains("Good") && target.Contains(Target) && (wien.Contains("UP") || wien.Contains("DOWN"))){
            Vals = OpenRun(rootfiles, run, slug, ihwp, wien, arm);
            //Vals2 = OpenRun(rootfiles2, run, ihwp, wien, arm);
            if (Vals[0][0] == -1) {
                continue;
            }

            //Vals[0][1] = Vals[0][1] - Vals2[0][1];
            total += Vals[0][1];
            //cout << run << " " << Vals[0][1] << endl;

            type_table.push_back(types_line);
            int minirun = 0;
            for(minirun = 0; minirun<Vals.size(); minirun++){
                Vals_table.push_back(Vals[minirun]);
                Miniruns.push_back(Vals[minirun][0]);
                RunMeans.push_back(Vals[minirun][1]);
                RunErrors.push_back(Vals[minirun][3]);
                
                Runs.push_back(run);
                Slugs.push_back(slug);
                Arms.push_back(arm);
                Targets.push_back(target);
            }

            /*/wiggly status message
            for (int j = 0; j < round(70 * (sin((i)*M_PI / 20)) + 70); j++) {
                cout << " ";
            }
            cout << "Found run " << run << " with " << minirun << " miniruns..." << endl;

            i++;
            xpos = i;*/
            
        }
    }

    infile.close();

    vector<TString> allrowvars;
    if(pcrex.Contains("PREX")){
        allrowvars = {
            "BPM_4aX",
            "BPM_4aY",
            "BPM_4eX",
            "BPM_4eY",
            "BPM_E"
        };
    }else if(pcrex.Contains("CREX")){
        allrowvars = {
            "BPM_1X",
            "BPM_4aY",
            "BPM_4eX",
            "BPM_4eY",
            "BPM_E"
        };
    }

    vector<TString> allcolvars = {
        "usl",
        "usr",
        "avg",
        "dd"
    };
    /*
    vector<vector<Double_t>> plotmatrix;
    TCanvas* c1 = new TCanvas();
    c1->Divide(2,1);
    Int_t ivar = 0;
    
    for(int rows = 0; rows < allrowvars.size(); rows++){
        for(int cols = 0; cols < allcolvars.size(); cols++){
            ivar = cols + rows*allcolvars.size();
            VarCombo = Form("%sVs%s", allrowvars[rows].Data(), allcolvars[cols].Data());
            c1->cd(1);
            for(int miniindex = 0; miniindex < Vals_table.size(); miniindex++){
                plotmatrix.push_back({Vals_table[miniindex][0], Vals_table[miniindex][(ivar*2)+1], 0, Vals_table[miniindex][(ivar*2)+2]});
            }
            TGraphErrors* RunPlot = PlotFromMatrix(plotmatrix, Form("%s, %s Target, Part %s, Miniruns %s", pcrex.Data(), Target.Data(), rungroup.Data(), VarCombo.Data()));
            RunPlot->Fit("pol0","q");
            gStyle->SetOptFit(1);
            RunPlot->Draw("AP");
            c1->cd(2);
            TH1D* RunHist = HistFromMatrix(plotmatrix, Form("%s, %s Target, Part %s, Miniuns %s", pcrex.Data(), Target.Data(), rungroup.Data(), VarCombo.Data()));
            RunHist->Fit("gaus","q");
            gStyle->SetOptFit(1);
            RunHist->Draw();
            plotmatrix.clear();

            c1->SaveAs(Form("./Outputs/temp/AllTargsAT%s%s%s%s.pdf", pcrex.Data(), Target.Data(), rungroup.Data(), VarCombo.Data()));
        }
    }

    vector<Int_t> Avusl;
    vector<Int_t> Avusr;
    vector<Int_t> Avavg;
    vector<Int_t> Avdd;
    if(pcrex.Contains("PREX")){
        Avusl = {3, 63, 123, 183, 243};
        Avusr = {5, 65, 125, 185, 245};
        Avavg = {9, 69, 129, 189, 249};
        Avdd = {11, 71, 131, 191, 251};
    }else if(pcrex.Contains("CREX")){
        Avusl = {1, 9, 17, 25, 33};
        Avusr = {3, 11, 19, 27, 35};
        Avavg = {5, 13, 21, 29, 37};
        Avdd = {7, 15, 23, 31, 39};
    }
   
    TCanvas* c2 = new TCanvas();
    c2->Divide(5,5);
    for(int rows1 = 0; rows1 < allrowvars.size(); rows1++){
        for(int rows2 = 0; rows2 < allrowvars.size(); rows2++){
            ivar = 1 + rows2 + rows1*allrowvars.size();
            VarCombo = Form("%sVsuslVs%sVsusl", allrowvars[rows1].Data(), allrowvars[rows2].Data());
            for(int miniindex = 0; miniindex < Vals_table.size(); miniindex++){
                c2->cd(ivar);
                plotmatrix.push_back({Vals_table[miniindex][Avusl[rows1]], Vals_table[miniindex][Avusl[rows2]], Vals_table[miniindex][Avusl[rows1]+1], Vals_table[miniindex][Avusl[rows2]+1]});
            }
            TGraphErrors* RunPlot = PlotFromMatrix(plotmatrix, Form("%s, %s Target, Part %s, Miniruns %s", pcrex.Data(), Target.Data(), rungroup.Data(), VarCombo.Data()));
            TF1* f1 = new TF1("f1", "[0] + [1]*x");
            f1->SetParameter(0, RunPlot->GetMean(2));
            RunPlot->Fit("f1","q");
            gStyle->SetOptFit(1);
            RunPlot->Draw("AP");
            plotmatrix.clear();
        }
    }
    c2->SaveAs(Form("./Outputs/temp/AllTargsAT%s%s%sAlluslSlopes.pdf", pcrex.Data(), Target.Data(), rungroup.Data()));

    TCanvas* c3 = new TCanvas();
    c3->Divide(5,5);
    for(int rows1 = 0; rows1 < allrowvars.size(); rows1++){
        for(int rows2 = 0; rows2 < allrowvars.size(); rows2++){
            ivar = 1 + rows2 + rows1*allrowvars.size();
            VarCombo = Form("%sVsusrVs%sVsusr", allrowvars[rows1].Data(), allrowvars[rows2].Data());
            for(int miniindex = 0; miniindex < Vals_table.size(); miniindex++){
                c3->cd(ivar);
                plotmatrix.push_back({Vals_table[miniindex][Avusr[rows1]], Vals_table[miniindex][Avusr[rows2]], Vals_table[miniindex][Avusr[rows1]+1], Vals_table[miniindex][Avusr[rows2]+1]});
            }
            TGraphErrors* RunPlot = PlotFromMatrix(plotmatrix, Form("%s, %s Target, Part %s, Miniruns %s", pcrex.Data(), Target.Data(), rungroup.Data(), VarCombo.Data()));
            TF1* f1 = new TF1("f1", "[0] + [1]*x");
            f1->SetParameter(0, RunPlot->GetMean(2));
            RunPlot->Fit("f1","q");
            gStyle->SetOptFit(1);
            RunPlot->Draw("AP");
            plotmatrix.clear();
        }
    }
    c3->SaveAs(Form("./Outputs/temp/AllTargsAT%s%s%sAllusrSlopes.pdf", pcrex.Data(), Target.Data(), rungroup.Data()));

    TCanvas* c4 = new TCanvas();
    c4->Divide(5,5);
    for(int rows1 = 0; rows1 < allrowvars.size(); rows1++){
        for(int rows2 = 0; rows2 < allrowvars.size(); rows2++){
            ivar = 1 + rows2 + rows1*allrowvars.size();
            VarCombo = Form("%sVsavgVs%sVsavg", allrowvars[rows1].Data(), allrowvars[rows2].Data());
            for(int miniindex = 0; miniindex < Vals_table.size(); miniindex++){
                c4->cd(ivar);
                plotmatrix.push_back({Vals_table[miniindex][Avavg[rows1]], Vals_table[miniindex][Avavg[rows2]], Vals_table[miniindex][Avavg[rows1]+1], Vals_table[miniindex][Avavg[rows2]+1]});
            }
            TGraphErrors* RunPlot = PlotFromMatrix(plotmatrix, Form("%s, %s Target, Part %s, Miniruns %s", pcrex.Data(), Target.Data(), rungroup.Data(), VarCombo.Data()));
            TF1* f1 = new TF1("f1", "[0] + [1]*x");
            f1->SetParameter(0, RunPlot->GetMean(2));
            RunPlot->Fit("f1","q");
            gStyle->SetOptFit(1);
            RunPlot->Draw("AP");
            plotmatrix.clear();
        }
    }
    c4->SaveAs(Form("./Outputs/temp/AllTargsAT%s%s%sAllavgSlopes.pdf", pcrex.Data(), Target.Data(), rungroup.Data()));

    TCanvas* c5 = new TCanvas();
    c5->Divide(5,5);
    for(int rows1 = 0; rows1 < allrowvars.size(); rows1++){
        for(int rows2 = 0; rows2 < allrowvars.size(); rows2++){
            ivar = 1 + rows2 + rows1*allrowvars.size();
            VarCombo = Form("%sVsddVs%sVsdd", allrowvars[rows1].Data(), allrowvars[rows2].Data());
            for(int miniindex = 0; miniindex < Vals_table.size(); miniindex++){
                c5->cd(ivar);
                plotmatrix.push_back({Vals_table[miniindex][Avdd[rows1]], Vals_table[miniindex][Avdd[rows2]], Vals_table[miniindex][Avdd[rows1]+1], Vals_table[miniindex][Avdd[rows2]+1]});
            }
            TGraphErrors* RunPlot = PlotFromMatrix(plotmatrix, Form("%s, %s Target, Part %s, Miniruns %s", pcrex.Data(), Target.Data(), rungroup.Data(), VarCombo.Data()));
            TF1* f1 = new TF1("f1", "[0] + [1]*x");
            f1->SetParameter(0, RunPlot->GetMean(2));
            RunPlot->Fit("f1","q");
            gStyle->SetOptFit(1);
            RunPlot->Draw("AP");
            plotmatrix.clear();

        }
    }
    c5->SaveAs(Form("./Outputs/temp/AllTargsAT%s%s%sAllddSlopes.pdf", pcrex.Data(), Target.Data(), rungroup.Data()));
    */

    if(pcrex.Contains("PREX")){
        for(int miniindex = 0; miniindex < Vals_table.size(); miniindex++){
            //cout << Vals_table[miniindex][11] << ", " << Vals_table[miniindex][12] << ", " << Vals_table[miniindex][71] << ", " << Vals_table[miniindex][72] << ", " << Vals_table[miniindex][191] << ", " << Vals_table[miniindex][192] << ", " << Vals_table[miniindex][251] << ", " << Vals_table[miniindex][252] << endl;
            //11, 71, 131, 191, 251 and errors
            //plotmatrix.push_back({Vals_table[miniindex][Avavg[rows1]], Vals_table[miniindex][Avavg[rows2]], Vals_table[miniindex][Avavg[rows1]+1], Vals_table[miniindex][Avavg[rows2]+1]});
        }
    }else if(pcrex.Contains("CREX")){
        //cout << Vals_table[miniindex][Avavg[rows1]] << endl;
        for(int miniindex = 0; miniindex < Vals_table.size(); miniindex++){
            //cout << Vals_table[miniindex][Avavg[rows1]] << endl;
            continue;
            //plotmatrix.push_back({Vals_table[miniindex][Avavg[rows1]], Vals_table[miniindex][Avavg[rows2]], Vals_table[miniindex][Avavg[rows1]+1], Vals_table[miniindex][Avavg[rows2]+1]});
        }
    }

    return;
}

vector<vector<Double_t>> OpenRun(TString directory, Int_t runnum, Int_t slugnum, TString ihwps, TString wiens, TString arm) {

    vector<vector<Double_t>> vals;
    //open the rootfile, crash if not available.
    TFile* runfile = TFile::Open(Form("%s/prexPrompt_pass2_%d.000.root", directory.Data(), runnum), "READ");
    if (runfile == NULL || arm != "0") {
        vals = {{ -1, -1 }};
        return vals;
    }

    //find the trees we need.
    TTree* mul = (TTree*)gROOT->FindObject("mul");
    TTree* mulc = (TTree*)gROOT->FindObject("mulc");
    TTree* mulc_lrb_alldet = (TTree*)gROOT->FindObject("mulc_lrb_alldet");
    TTree* mulc_dit_combo = (TTree*)gROOT->FindObject("mulc_dit_combo");
    TTree* mulc_dit = (TTree*)gROOT->FindObject("mulc_dit");
    TTree* burst_mulc_lrb_burst = (TTree*)gROOT->FindObject("burst_mulc_lrb_burst");
    TTree* burst_lrb_std = (TTree*)gROOT->FindObject("burst_lrb_std");

    mul->AddFriend("mulc");
    mul->AddFriend("mulc_lrb_alldet");
    mul->AddFriend("mulc_dit_combo");
    
    TBranch* b2 = mul->GetBranch("diff_bpm4eX");
    TBranch* b5 = mul->GetBranch("diff_bpm4eY");
    TBranch* b6 = mul->GetBranch("diff_bpm4aX");
    TBranch* b7 = mul->GetBranch("diff_bpm4aY");
    TBranch* b8 = mulc->GetBranch("diff_bpmE");
    TBranch* b9 = mulc->GetBranch("diff_bpm1X");
    TBranch* b10 = mulc->GetBranch("diff_bpm12X");
    TBranch* b11 = burst_lrb_std->GetBranch("|statA");
    TBranch* b12 = burst_lrb_std->GetBranch("|statdA");

    TBranch* b1 = mulc_lrb_alldet->GetBranch("cor_asym_us_dd");

    TBranch* b3 = mul->GetBranch("BurstCounter");
    TBranch* b4 = burst_mulc_lrb_burst->GetBranch("cor_usr");

    if(b1 == NULL || b2 == NULL || b3 == NULL){
        //cout << "Cant find branches in tree for run " << runnum << endl;
        vals = {{ -1, -1 }};
        return vals;
    }
    
    TLeaf* asymvalues = b1->GetLeaf(ValueLeaf);

    TLeaf* varex = b2->GetLeaf(ValueLeaf);

    TLeaf* varey = b5->GetLeaf(ValueLeaf);

    TLeaf* varax = b6->GetLeaf(ValueLeaf);

    TLeaf* varay = b7->GetLeaf(ValueLeaf);

    TLeaf* vare = b8->GetLeaf(ValueLeaf);

    TLeaf* slopeA = b11->GetLeaf("A");
    TLeaf* slopedA = b12->GetLeaf("dA");

    TLeaf* bcounter = b3->GetLeaf("BurstCounter");
    Int_t bursts = b4->GetEntries();
    Int_t entries = b1->GetEntries();
    

    //selects the sign.
    Int_t ihwp = 1;
    Int_t wien = 1;
    if (ihwps.Contains("IN")) {
        ihwp = -1;
    }
    if (wiens.Contains("RIGHT")) {
        wien = -1;
    }
    
    vector<vector<TH1F*>> temp;
    for(int j = 0; j<bursts; j++){
        temp.push_back({new TH1F, new TH1F, new TH1F, new TH1F, new TH1F});
        if(runnum==4117 && j==0){
            continue;
        }

        burst_lrb_std->GetEntry(j);

        if(runnum>=5337){
            vector<Double_t> TempRunVector;
            TempRunVector.push_back(runnum+j/static_cast<double>(bursts));
            for(int ii = 0; ii < 20; ii++){
                TempRunVector.push_back(slopeA->GetValue(ii));
                TempRunVector.push_back(slopedA->GetValue(ii));
            }
            vals.push_back(TempRunVector);
            TempRunVector.clear();
        }else{
            vector<Double_t> TempRunVector;
            TempRunVector.push_back(runnum+j/static_cast<double>(bursts));
            for(int ii = 0; ii < 150; ii++){
                TempRunVector.push_back(slopeA->GetValue(ii));
                TempRunVector.push_back(slopedA->GetValue(ii));
            }
            vals.push_back(TempRunVector);
            cout << runnum << ", " << j << ", " << TempRunVector[11] << ", " << TempRunVector[12] << ", " << TempRunVector[71] << ", " << TempRunVector[72] << ", " << TempRunVector[131] << ", " << TempRunVector[132] << ", " << TempRunVector[191] << ", " << TempRunVector[192] << ", " << TempRunVector[251] << ", " << TempRunVector[252] << endl;
            TempRunVector.clear();
        }
    }

    runfile->Close();

    return vals;
}

vector<Double_t> GetSlugVals(Int_t slugnum, vector<Double_t> miniruns, vector<Int_t> slugs, vector<TString> arms, vector<Double_t> means, vector<Double_t> errors) {
    vector<Double_t> slugmeanerr = { -1,-1 };

    vector<vector<Double_t>> slugouts;
    //TH2D* plot = new TH2D("", "", 500, 1, 0, 500, 1, 0);

    Double_t numerator = 0;
    Double_t denominator = 0;

    //loops through miniruns in list for ones that match slugnum
    for (Int_t i = 0; i < miniruns.size(); i++) {
        if (slugs[i] == slugnum) {
            //if(!arms[i].Contains("0")){
            //  continue;
            //}
            numerator += means[i]/(errors[i]*errors[i]);
            denominator += 1.0/(errors[i]*errors[i]);

            //plot->Fill(miniruns[i], means[i]);
            slugouts.push_back({ (Double_t)miniruns[i], means[i], 0, errors[i] });
        }
    }

    if (slugouts.size() < 1) {
        return slugmeanerr;
    }

    slugmeanerr = { (Double_t)slugnum, numerator/denominator, 0, sqrt(1.0/denominator) };

    return slugmeanerr;
}

TGraphErrors* PlotFromMatrix(vector<vector<Double_t>> data, TString title) {

    Double_t x[data.size()];
    Double_t y[data.size()];
    Double_t ex[data.size()];
    Double_t ey[data.size()];
    for (Int_t i = 0; i < data.size(); i++) {
        x[i] = data[i][0];
        y[i] = data[i][1];
        ex[i] = data[i][2];
        ey[i] = data[i][3];
    }

    TGraphErrors* output = new TGraphErrors(data.size(), x, y, ex, ey);
    output->SetTitle(title);
    output->SetMarkerStyle(7);

    return output;
}

TH1D* PullFromMatrix(vector<vector<Double_t>> data, TString title) {

    Double_t x[data.size()];
    Double_t y[data.size()];
    Double_t ex[data.size()];
    Double_t ey[data.size()];
    for (Int_t i = 0; i < data.size(); i++) {
        x[i] = data[i][0];
        y[i] = data[i][1];
        ex[i] = data[i][2];
        ey[i] = data[i][3];
    }

    TGraphErrors* meanhist = new TGraphErrors(data.size(), x, y, ex, ey);

    TH1D* output = new TH1D(title, title, 50, 1, 0);
    for (Int_t i = 0; i < data.size(); i++) {
        // fix me use get mean on y axis output->Fill((y[i]-meanhist->GetMean())/ey[i]);
        output->Fill((y[i]-meanhist->GetMean(2))/ey[i]);
    }

    output->SetTitle(title);

    return output;
}

TH1D* HistFromMatrix(vector<vector<Double_t>> data, TString title) {

    TH1D* output = new TH1D(title, title, 50, 1, 0);

    //Double_t x[data.size()];
    for (Int_t i = 0; i < data.size(); i++) {
        //x[i] = data[i][1];
        output->Fill(data[i][1]);
    }

    output->SetTitle(title);

    return output;
}
