//##This script is deigned to loop through a runlist and calculate means and errors of values in those runs. Currently setup to find the asymmetry across all runs.

TString promptdir = "/u/group/halla/parity/software/japan_offline/prompt/prex-prompt";
TString rootfiles = "$QW_ROOTFILES";
TString rootfiles2 = "/volatile/halla/parity/crex-respin1/japanOutput";

//goal settings
TString Tree = "burst_mulc_lrb_alldet";
TString Branch = "cor_asym_us_avg";
TString ValueLeaf = "hw_sum";
TString ErrorLeaf = "hw_sum_err";
TString outfilePath = "./ComparisonOutputs/outputca48old.txt";

/*//postpan identical
TString Tree = "burst_mulc_lrb_burst";
TString Branch = "cor_asym_us_dd";
TString ValueLeaf = "hw_sum";
TString ErrorLeaf = "hw_sum_err";*/

/*
TString Tree = "burst_mulc_lrb";
TString Branch = "cor_asym_us_dd";
TString ValueLeaf = "hw_sum";
TString ErrorLeaf = "hw_sum_err";*/

/*
TString Tree = "burst_mulc";
TString Branch = "asym_us_dd";
TString ValueLeaf = "hw_sum";
TString ErrorLeaf = "hw_sum_err";*/

//mostly used for reference; all of the rcdb data types in order, with an added entry for run number at the end.
vector<TString> types = { "arm_flag", "beam_current", "beam_energy", "bmw", "components", "component_stats", "event_count", "event_rate", "experiment", "feedback", "FFB", "flip_state", "good_charge", " helicity_frequency", "helicity_pattern", "horizontal_wien", "ihwp", "is_valid_run_end", "prompt_analysis", "respin_comment", "rhwp", "rtvs", "run_config", "run_end_time", "run_flag", "run_length", "run_prestart_time", "run_start_epoch", "run_start_time", "run_type", "session", "slug", "target_45encoder", "target_90encoder", "target_encoder", "target_type", "total_charge", "user_comment", "vertical_wien", "wac_comment", "run_number" };
vector<TString> types_line = { " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " " };

//opens the rootfile for a a given run, sign corrects, then returns the mean and error.
vector<vector<Double_t>> OpenRun(TString directory, Int_t runnum, Int_t slugnum, TString ihwps, TString wiens, TString arm, vector<vector<Double_t>> textmatrix, Int_t WindowSize);

//averages the means for miniruns in a given slug, returns the slug mean and error.
vector<Double_t> GetSlugVals(Int_t slugnum, vector<Double_t> miniruns, vector<Int_t> slugs, vector<TString> arms, vector<Double_t> means, vector<Double_t> errors);

//takes a vector of vectors and returns a tgrapherrors. (lol why is this not built in?)
TGraphErrors* PlotFromMatrix(vector<vector<Double_t>> data, TString title);
  
//takes the matrix format mentioned above and plots the second column(means)
TH1D* HistFromMatrix(vector<vector<Double_t>> data, TString title);

void StatsCheckanus() {

    vector<vector<TString>> type_table;
    remove(outfilePath.Data());

    vector<Int_t> SizeList = {1000, 3000, 6000, 9000};
    vector<TGraphErrors*> PlotList;
    vector<TH1D*> GausList;
    TCanvas* c1 = new TCanvas();
    c1->Divide(2,SizeList.size());
    for(int isize; isize < SizeList.size(); isize++){
 
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
    Int_t i = 0;
    Int_t First = 5400; // 3305
    //Int_t Last = 5420;
    Int_t Last =5500;//8559; // 4980
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
        if (run >= First && run <= Last && run != 3140 && production.Contains("Production") && goodbad.Contains("Good") && target.Contains("48") && (wien.Contains("RIGHT") || wien.Contains("LEFT"))){
            Vals = OpenRun(rootfiles, run, slug, ihwp, wien, arm, textmatrix, SizeList[isize]);
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

            //wiggly status message
            for (int j = 0; j < round(70 * (sin((i)*M_PI / 20)) + 70); j++) {
                cout << " ";
            }
            cout << "Found run " << run << " with " << minirun << " miniruns..." << endl;

            i++;
            
        }
    }

    infile.close();

    PlotList.push_back(PlotFromMatrix(Vals_table, Form("PREX an_us, %d Window", SizeList[isize])));
    PlotList[isize]->Fit("pol0");

    c1->cd(1+(isize*2));
    gStyle->SetOptFit(1);
    PlotList[isize]->Draw("AP");
    c1->cd(2+(isize*2));
    GausList.push_back(HistFromMatrix(Vals_table, Form("PREX an_us, %d Window", SizeList[isize])));
    GausList[isize]->Fit("gaus");
    gStyle->SetOptFit(1);
    GausList[isize]->Draw();

    }

    //cout << endl << endl << "Total difference in ErrorFlag==0 events: " << total << endl;
    /*
    //loops through all possible slug numbers included in your runs. this is fast enough to not matter.
    for (Int_t p = 0; p < 5000; p++) {
        slugouts = GetSlugVals(p, Miniruns, Slugs, Arms, RunMeans, RunErrors);
        if (slugouts[0] == -1) {
            continue;
        }
        slugouts_table.push_back(slugouts);
    }

    TGraphErrors* SlugPlot = PlotFromMatrix(slugouts_table, Form("PREX Slugs %s", Branch.Data()));
    SlugPlot->Fit("pol0");

    TCanvas* c2 = new TCanvas();
    c2->Divide(2,1);
    c2->cd(1);
    gStyle->SetOptFit(1);
    SlugPlot->Draw("AP");
    c2->cd(2);
    TH1D* SlugHist = HistFromMatrix(slugouts_table, Form("PREX Slugs %s", Branch.Data()));
    SlugHist->Fit("gaus");
    gStyle->SetOptFit(1);
    SlugHist->Draw();
    */
    return;
}

vector<vector<Double_t>> OpenRun(TString directory, Int_t runnum, Int_t slugnum, TString ihwps, TString wiens, TString arm, vector<vector<Double_t>> textmatrix, Int_t WindowSize) {

    vector<vector<Double_t>> vals;
    //open the rootfile, crash if not available.
    TFile* runfile = TFile::Open(Form("%s/prexPrompt_pass2_%d.000.root", directory.Data(), runnum), "READ");
    if (runfile == NULL) {
        vals = {{ -1, -1 }};
        return vals;
    }

    //find the trees we need.
    TTree* mul = (TTree*)gROOT->FindObject("mul");
    TTree* mulc_lrb_alldet = (TTree*)gROOT->FindObject("mulc_lrb_alldet");
    TTree* mulc_lrb_burst = (TTree*)gROOT->FindObject("mulc_lrb_burst");
    mul->AddFriend("mulc_lrb_alldet");
    mul->AddFriend("mulc_lrb_burst");
    
    TTree* burst_mulc_lrb_alldet = (TTree*)gROOT->FindObject("burst_mulc_lrb_alldet");
    TTree* burst_mulc_lrb_burst = (TTree*)gROOT->FindObject("burst_mulc_lrb_burst");
    TTree* lrb_alldet = (TTree*)gROOT->FindObject("lrb_alldet");

    TBranch* b1;
    TBranch* b2;
    if(arm.Contains("0")){
        b1 = burst_mulc_lrb_alldet->GetBranch("cor_asym_us_avg");
        b2 = burst_mulc_lrb_burst->GetBranch("cor_asym_us_avg");
    }else if(arm.Contains("1")){
        b1 = burst_mulc_lrb_alldet->GetBranch("cor_usr");
        b2 = burst_mulc_lrb_burst->GetBranch("cor_usr");
    }else if(arm.Contains("2")){
        b1 = burst_mulc_lrb_alldet->GetBranch("cor_usl");
        b2 = burst_mulc_lrb_burst->GetBranch("cor_usl");
    }

    TBranch* b3 = lrb_alldet->GetBranch("good_count");
    //TBranch* b4 = mulc_lrb_burst->GetBranch("cor_usl");
    TBranch* b4 = mul->GetBranch("asym_bcm_an_us");
    TBranch* b5 = mulc_lrb_burst->GetBranch("cor_usr");
    TBranch* b6 = mul->GetBranch("ErrorFlag");

    if(b1 == NULL || b2 == NULL || b3 == NULL){
        //cout << "Cant find branches in tree for run " << runnum << endl;
        vals = {{ -1, -1 }};
        return vals;
    }
    
    TLeaf* values = b1->GetLeaf(ValueLeaf);
    TLeaf* errors = b1->GetLeaf(ErrorLeaf);
    TLeaf* count = b1->GetLeaf("num_samples");

    TLeaf* panvalues = b2->GetLeaf(ValueLeaf);
    TLeaf* panerrors = b2->GetLeaf(ErrorLeaf);
    
    TLeaf* goodcount = b3->GetLeaf("good_count");

    TLeaf* values_usl = b4->GetLeaf(ValueLeaf);
    //TLeaf* errors_usl = b4->GetLeaf(ErrorLeaf);
    //TLeaf* count_usl = b4->GetLeaf("num_samples");

    TLeaf* values_usr = b5->GetLeaf(ValueLeaf);
    //TLeaf* errors_usr = b5->GetLeaf(ErrorLeaf);
    //TLeaf* count_usr = b5->GetLeaf("num_samples");

    TLeaf* errorflag = b6->GetLeaf("ErrorFlag");

    //TLeaf* usl = b4->GetLeaf("hw_sum");
    TH1F* cor;
    mul->Draw("mulc_lrb_burst.cor_usl>>cor()", "ErrorFlag==0", "goff");
    cor = (TH1F*)gROOT->FindObject("cor");
    
    Int_t events = cor->GetEntries();
    //Int_t entries = b1->GetEntries();

    //selects the sign.
    Int_t ihwp = 1;
    Int_t wien = 1;
    if (ihwps.Contains("IN")) {
        ihwp = -1;
    }
    if (wiens.Contains("RIGHT")) {
        wien = -1;
    }
    
    Int_t entries = ceil(events/WindowSize);

    Bool_t combinestragglers = false;
    Bool_t lastmini = false;
    
    Int_t lastcount = 0;
    Int_t thiscount = 0;

    if(entries==0){
        entries = 1;
    }

    for(int j = 0; j < entries; j++){
        TH1F* main = new TH1F();
        Int_t eventcount = 0;

        for(int k = lastcount; k < b4->GetEntries(); k++){
            b4->GetEntry(k);
            b5->GetEntry(k);
            b6->GetEntry(k);
            if(errorflag->GetValue(0) == 0){
                main->Fill(values_usl->GetValue(0));
                eventcount++;
                if(eventcount == WindowSize){
                    thiscount = k;
                    if(j==entries-2){
                        if(combinestragglers==true){
                            thiscount = b4->GetEntries()-1;
                        }else{
                            lastmini = true;
                        }
                    }
                    break;
                }
            }
            if(k == b4->GetEntries()-2){
                thiscount = k;
            }
        }

        vals.push_back({runnum+j/static_cast<double>(entries), main->GetMeanError(), 0, 0});
        main->Delete();

        if(lastmini){
            j++;
        }

        //cout << lastcount << " " << thiscount << endl;
        lastcount = thiscount;
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
