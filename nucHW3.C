void nucHW3(){

    Double_t A = 0;
    Double_t Z = 0;
    Double_t BE = 0;
    vector<TString> types_line;
    vector<vector<TString>> type_table;
    //TH3D *h3 = new TH3D("A vs Z vs BE 3d", "A vs Z vs BE 3d", 270, 0, 271, 110, 0, 111, 1000, 1, -1);
    TGraph2DErrors *t3 = new TGraph2DErrors(2438);
    TH2D *h2 = new TH2D("A vs Z vs BE", "A vs Z vs BE", 270, 0, 271, 110, 0, 111);
    TH2D *h2fit = new TH2D("Fit Result A vs Z vs BE", "Fit Result A vs Z vs BE", 270, 0, 271, 110, 0, 111);
    TH2D *h2res = new TH2D("Fit Residual A vs Z vs BE", "Fit Residual A vs Z vs BE", 270, 0, 271, 110, 0, 111);
    TH2D *h2rescent = new TH2D("Fit Residual  A vs Z vs BE Percent Diff", "Fit Residual A vs Z vs BE Percent Diff", 270, 0, 271, 110, 0, 111);

    ifstream infile;
    infile.open("./AME2012.txt");

    string ins;

    bool isemptyline = false;
    Int_t point = 0;

    while (!infile.eof()) {
        //splits the file at line breaks
        getline(infile, ins, '\n');
        stringstream line(ins);
        types_line.clear();
        //splits the line at commas
        while (line.good()) {
            string substr;
            getline(line, substr, ',');
            if(substr == ""){
                isemptyline = true;
                break;
            }
            types_line.push_back(substr);
        }

        if(isemptyline == true){
            break;
        }

        A = stod(types_line[2].Data());
        Z = stod(types_line[1].Data());
        BE = stod(types_line[6].Data())*A;

        //h3->Fill(A, Z, BE);
        t3->SetPoint(point, A, Z, BE);
        t3->SetPointError(point, 0, 0, stod(types_line[7].Data())*A);
        h2->Fill(A, Z, BE);
        point++;
        //cout << types_line[0] << " " << types_line[1] << " " << types_line[6] << "+-" << types_line[7] << " " << endl;
            
    }

    infile.close();

    TCanvas *c1 = new TCanvas();
    c1->Divide(2,2);
    c1->cd(1);
    TF2* f2 = new TF2("f2", "[0]*x + [1]*(x^(2/3)) + [2]*y*(y-1)/(x^(1/3)) + [3]*((y-x/2)^2)/x + [4]*(sqrt(x)^(-1))");
    t3->Fit("f2");
    gStyle->SetOptStat(0);
    //gPad->SetLogz(1);
    h2->Draw("colz");

    infile.open("./AME2012.txt");

    isemptyline = false;

    Double_t BEmat[270][271];
    for(Int_t A = 0; A < 270; A++){
        for(Int_t Z = 0; Z < 271; Z++){
            BEmat[A][Z] = (10^100);
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
            if(substr == ""){
                isemptyline = true;
                break;
            }
            types_line.push_back(substr);

        }

        if(isemptyline == true){
            break;
        }

        A = stod(types_line[2].Data());
        Z = stod(types_line[1].Data());
        BE = stod(types_line[6].Data())*A;

        BEmat[stoi(types_line[2].Data())][stoi(types_line[1].Data())] = stod(types_line[6].Data())*A;

        h2fit->Fill(A, Z, f2->Eval(A, Z));
        h2res->Fill(A, Z, BE - f2->Eval(A, Z));
        h2rescent->Fill(A, Z, abs((BE - f2->Eval(A, Z))/(BE + f2->Eval(A, Z))));   

    }

    infile.close();

    c1->cd(2);
    //gPad->SetLogz(1);
    h2fit->Draw("colz");
    c1->cd(3);
    //gPad->SetLogz(1);
    h2res->Draw("colz");
    c1->cd(4);
    //gPad->SetLogz(1);
    h2rescent->Draw("colz");

    TH2D *h2bp = new TH2D("A vs Z vs B+ Q Values", "A vs Z vs B+ Q Values", 270, 0, 271, 110, 0, 111);
    TH2D *h2bm = new TH2D("A vs Z vs B- Q Values", "A vs Z vs B- Q Values", 270, 0, 271, 110, 0, 111);
    TH2D *h2bpstab = new TH2D("A vs Z vs B+ stability", "A vs Z vs B+ stability", 270, 0, 271, 110, 0, 111);
    TH2D *h2bmstab = new TH2D("A vs Z vs B- stability", "A vs Z vs B- stability", 270, 0, 271, 110, 0, 111);
    TH2D *h2stab = new TH2D("A vs Z vs stability", "A vs Z vs stability", 270, 0, 271, 110, 0, 111);

    for(Int_t A = 0; A < 270; A++){
        for(Int_t Z = 0; Z < 271; Z++){
            if(BEmat[A][Z] != (10^100) && BEmat[A][Z-1] != (10^100)){
                h2bp->Fill(A, Z, BEmat[A][Z] - BEmat[A][Z-1]);
            }
            if(BEmat[A][Z] != (10^100) && BEmat[A][Z+1] != (10^100)){
                h2bm->Fill(A, Z, BEmat[A][Z] - BEmat[A][Z+1]);
            }

            if(BEmat[A][Z] != (10^100)){
                h2bpstab->Fill(A, Z, 1);
                h2bmstab->Fill(A, Z, 1);
                h2stab->Fill(A, Z, 1);
            }

            if(BEmat[A][Z] != (10^100) && BEmat[A][Z-1] != (10^100) && BEmat[A][Z] - BEmat[A][Z-1] > 0){
                h2bpstab->Fill(A, Z, 1);
            }

            if(BEmat[A][Z] != (10^100) && BEmat[A][Z+1] != (10^100) && BEmat[A][Z] - BEmat[A][Z+1] > 0){
                h2bmstab->Fill(A, Z, 1);
            }

            if(BEmat[A][Z] != (10^100) && BEmat[A][Z-1] != (10^100) && BEmat[A][Z+1] != (10^100) && (BEmat[A][Z] - BEmat[A][Z+1] > 0) && (BEmat[A][Z] - BEmat[A][Z+1] > 0)){
                h2stab->Fill(A, Z, 1);
            }
        }
    }

    TCanvas *c2 = new TCanvas();
    c2->Divide(2,2);
    c2->cd(1);
    h2bp->Draw("colz");
    c2->cd(2);
    h2bm->Draw("colz");
    c2->cd(3);
    h2bpstab->Draw("colz");
    c2->cd(4);
    h2bmstab->Draw("colz");

    TCanvas *c3 = new TCanvas();
    h2stab->Draw("colz");

    TF1* f1m = new TF1("f1m", "((192-x)*(x))^2", 0, 10000);
    TF1* f1p = new TF1("f1p", "((-2432-x)*(x))^2", 0, 10000);

    TCanvas *c4 = new TCanvas();

    f1m->Draw();
    f1p->Draw("same");

    return;
}
