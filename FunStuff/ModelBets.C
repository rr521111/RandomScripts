vector<TString> types_line;
string ins;
vector<vector<Double_t>> games;

vector<vector<Double_t>> Model1(vector<vector<Double_t>> data, Double_t cutoff, Double_t betsize);

TGraphErrors* PlotFromMatrix(vector<vector<Double_t>> data, TString title);

void ModelBets() {
    Int_t i = 0;
    Double_t maxindex = -1;
    Double_t minindex = 99999999999999;
    Double_t maxmult = -1;
    Double_t minmult = 99999999999999;
    Double_t index = 0;
    Double_t mult = 0;
    Double_t random = 0;
    Double_t numcrash = 0;
    srand((unsigned) time(0));
    ifstream infile;
    infile.open("./crash_game_roobet.csv");
    while (!infile.eof()) {
        getline(infile, ins, '\n');
        stringstream line(ins);
        types_line.clear();
        //splits the line at commas
        while (line.good()) {
            string substr;
            getline(line, substr, ',');
            types_line.push_back(substr);
        }

        if (types_line.size() <= 1) {
            break;
        }

        index = stod(types_line[0].Data());
        mult = stod(types_line[1].Data());

        if(i%100000==0){
            cout << "Read " << i << " events." << endl;
        }

        random = (rand() % 100) + 1;
        if(random <= 3){
            games.push_back({index, 0, 0, 0});
            numcrash++;
        }else{
            games.push_back({index, mult, 0, 0 });
        }

        if(mult>maxmult){
            maxmult=mult;
        }
    
        if(mult<minmult){
            minmult=mult;
        }

        if(index>maxindex){
            maxindex=index;
        }
        
        if(index<minindex){
            minindex=index;
        }
        i++;
    }
    cout << "Crash rate " << numcrash/i << endl;

    //Just all crashes in bins.
    TH1D* gamehist = new TH1D("First Million Crashes", "First Million Crashes", 1000, minmult-1, maxmult+1);
    for(Int_t i = 0; i < games.size(); i++){
        gamehist->Fill(games[i][1]);
    }

    TCanvas *c1 = new TCanvas();
    gamehist->Draw();

    //Crashes over time.
    TH2D* gamehist2d = new TH2D("First Million Crashes Over Time", "First Million Crashes Over Time", 1000, minindex-1, maxindex+1, 1000, minmult-1, maxmult+1);
    for(Int_t i = 0; i < games.size(); i++){
        gamehist2d->Fill(games[i][0], games[i][1]);
    }

    TCanvas *c2 = new TCanvas();
    gamehist2d->Draw();
    gamehist2d->Fit("pol0", "Q");
    gStyle->SetOptFit(1);

    //Taylor's game.
    Double_t maxearn = -99999999999;
    Double_t minearn = 99999999999;
    vector<vector<Double_t>> Model1Earnings = Model1(games, 1.36, 1);
    for(Int_t i = 0; i < Model1Earnings.size(); i++){
        if(Model1Earnings[i][1]>maxearn){
                maxearn=Model1Earnings[i][1];
        }
        if(Model1Earnings[i][1]<minearn){
            minearn=Model1Earnings[i][1];
        }
    }

    TH2D* Model1hist2d = new TH2D("Model1 Earnings", "Model1 Earnings", 1000, minindex-1, maxindex+1, 1000, minearn-1, maxearn+1);
    for(Int_t i = 0; i < Model1Earnings.size(); i++){
        Model1hist2d->Fill(Model1Earnings[i][0], Model1Earnings[i][1]);
    }

    TCanvas *c3 = new TCanvas();
    TF1 *fit1 = new TF1("fit1", "pol1");
    Model1hist2d->Fit("fit1", "Q");
    Model1hist2d->Draw();
    gStyle->SetOptFit(1);

    cout << fit1->GetParameter(1) << " " << fit1->GetParError(1) << endl;

    //Optimizing model1.
    //vector<Double_t> cutoffs = {1.05, 1.10, 1.15, 1.20, 1.25, 1.30, 1.35, 1.40, 1.45, 1.50, 1.55, 1.60, 1.65, 1.70, 1.75, 1.80, 1.85, 1.90, 1.95, 2.00};
    vector<Double_t> cutoffs;
    Double_t steps = 64;
    Double_t mincut = 1;
    Double_t maxcut = 4;
    for(Int_t i = 0; i < steps; i++){
        cutoffs.push_back(mincut+((i+1)*(maxcut-mincut)/steps));
    }

    vector<Double_t> bets = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    vector<vector<Double_t>> EarningsTable;
    TCanvas *c4 = new TCanvas();
    c4->Divide(8,8);
    for(Int_t j = 0; j < cutoffs.size(); j++){
        c4->cd(j+1);
        Double_t maxearn = -99999999999;
        Double_t minearn = 99999999999;
        vector<vector<Double_t>> Model1Earnings = Model1(games, cutoffs[j], 1);
        for(Int_t i = 0; i < Model1Earnings.size(); i++){
            if(Model1Earnings[i][1]>maxearn){
                maxearn=Model1Earnings[i][1];
            }
            if(Model1Earnings[i][1]<minearn){
                minearn=Model1Earnings[i][1];
            }
        }

        Model1hist2d = new TH2D("Model1 Earnings", "Model1 Earnings", 1000, minindex-1, maxindex+1, 1000, minearn-1, maxearn+1);
        for(Int_t i = 0; i < Model1Earnings.size(); i++){
            Model1hist2d->Fill(Model1Earnings[i][0], Model1Earnings[i][1]);
        }

        Model1hist2d->Fit("fit1", "Q");
        Model1hist2d->Draw();

        EarningsTable.push_back({cutoffs[j], fit1->GetParameter(1), 0, fit1->GetParError(1)});
    }
    
    TGraphErrors *Optimization = PlotFromMatrix(EarningsTable, "Mean Earnings");
    TCanvas *c5 = new TCanvas();
    Optimization->Draw("AP");
    //Optimization->Fit("pol1", "Q");
    gStyle->SetOptFit(1);
}

vector<vector<Double_t>> Model1(vector<vector<Double_t>> data, Double_t cutoff, Double_t betsize){
    vector<vector<Double_t>> earnings;
    Int_t wins = 0;
    Int_t losses = 0;

    earnings.push_back({0, data[0][1]*betsize, 0, 0});
    for(Int_t i = 1; i < data.size(); i++){
        if(data[i][1]>cutoff){
            earnings.push_back({static_cast<Double_t>(i), earnings[i-1][1]+data[i][1]*betsize, 0, 0});
            wins++;
        }else{
            earnings.push_back({static_cast<Double_t>(i), earnings[i-1][1]-betsize, 0, 0});
            losses++;
        }
    }
    cout << "Win/loss for cut: " << cutoff << " with bet size: " << betsize << " is " << wins << "/" << losses << endl;

    return earnings;
};

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
