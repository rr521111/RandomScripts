void ReadingTextOutputs.C(){

    vector<Double_t> runs1;
    vector<Int_t> counts1;
    vector<Double_t> asymmetries1;
    vector<Double_t> errors1;

    vector<Double_t> runs2;
    vector<Int_t> counts2;
    vector<Double_t> asymmetries2;
    vector<Double_t> errors2;
    
    vector<TString> types_line = { " ", " ", " ", " "};

    ifstream infile;
    infile.open("./outputold.txt");
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
        
        runs1.push_back(stod(types_line[0]));
        counts1.push_back(stoi(types_line[1]));
        asymmetries1.push_back(stod(types_line[2]));
        errors1.push_back(stod(types_line[3]));
        
    }

    infile.close();

    ifstream infile2;
    infile2.open("./outputnew.txt");
    while (!infile2.eof()) {
        getline(infile2, ins, '\n');
        stringstream line(ins);
        types_line.clear();
        //splits the line at commas
        while (line.good()) {
            string substr;
            getline(line, substr, ',');
            types_line.push_back(substr);
        }
        
        runs2.push_back(stod(types_line[0]));
        counts2.push_back(stoi(types_line[1]));
        asymmetries2.push_back(stod(types_line[2]));
        errors2.push_back(stod(types_line[3]));
        
    }

    infile2.close();

    
}
