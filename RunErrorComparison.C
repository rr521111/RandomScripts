TString filepath1 = "/chafs2/work1/apar/japanOutput/";
TString filepath2 = "/volatile/halla/parity/crex-respin1/japanOutput/";

void RunErrorComparison(Int_t Run){
    TFile* file1 = TFile::Open(Form("%s/prexPrompt_pass2_%d.000.root", filepath1.Data(), Run), "READ");
    TFile* file2 = TFile::Open(Form("%s/prexPrompt_pass2_%d.000.root", filepath2.Data(), Run), "READ");

    file1->cd();

    TTree* mul = (TTree*)gROOT->FindObject("mul");

    mul->AddFriend("new = mul",file2);

    TCanvas *c1 = new TCanvas();
    mul->Draw("new.asym_usl:Entry$>>usl_main2D()","new.ErrorFlag==0","goff");
    TH2F *usl_main2D = (TH2F*)gROOT->FindObject("usl_main2D");
    mul->Draw("new.asym_usl:Entry$>>usl_new2D()","ErrorFlag!=0 && new.ErrorFlag==0","goff");
    TH2F *usl_new2D = (TH2F*)gROOT->FindObject("usl_new2D");
    usl_new2D->SetMarkerStyle(7);
    usl_new2D->SetMarkerColor(2);
    usl_main2D->Draw();
    usl_new2D->Draw("same");
    c1->SaveAs("./ComparisonOutputs/Plots/temp/usl2D.pdf");

    TCanvas *c2 = new TCanvas();
    c2->Divide(2,1);
    mul->Draw("new.asym_usl>>usl_main()","new.ErrorFlag==0","goff");
    TH2F *usl_main = (TH2F*)gROOT->FindObject("usl_main");
    mul->Draw("new.asym_usl>>usl_new()","ErrorFlag!=0 && new.ErrorFlag==0","goff");
    TH2F *usl_new = (TH2F*)gROOT->FindObject("usl_new");
    c2->cd(1);
    usl_main->Draw();
    c2->cd(2);
    usl_new->Draw();
    c2->SaveAs("./ComparisonOutputs/Plots/temp/usl.pdf");


    TCanvas *c3 = new TCanvas();
    mul->Draw("new.asym_bcm_an_us:Entry$>>bcm_main2D()","new.ErrorFlag==0","goff");
    TH2F *bcm_main2D = (TH2F*)gROOT->FindObject("bcm_main2D");
    mul->Draw("new.asym_bcm_an_us:Entry$>>bcm_new2D()","ErrorFlag!=0 && new.ErrorFlag==0","goff");
    TH2F *bcm_new2D = (TH2F*)gROOT->FindObject("bcm_new2D");
    bcm_new2D->SetMarkerStyle(7);
    bcm_new2D->SetMarkerColor(2);
    bcm_main2D->Draw();
    bcm_new2D->Draw("same");
    c3->SaveAs("./ComparisonOutputs/Plots/temp/bcm2D.pdf");

    TCanvas *c4 = new TCanvas();
    c4->Divide(2,1);
    mul->Draw("new.asym_bcm_an_us>>bcm_main()","new.ErrorFlag==0","goff");
    TH2F *bcm_main = (TH2F*)gROOT->FindObject("bcm_main");
    mul->Draw("new.asym_bcm_an_us>>bcm_new()","ErrorFlag!=0 && new.ErrorFlag==0","goff");
    TH2F *bcm_new = (TH2F*)gROOT->FindObject("bcm_new");
    c4->cd(1);
    bcm_main->Draw();
    c4->cd(2);
    bcm_new->Draw();
    c4->SaveAs("./ComparisonOutputs/Plots/temp/bcm.pdf");


    TCanvas *c5 = new TCanvas();
    mul->Draw("new.yield_bpm4aX:Entry$>>bpm4aX_main2D()","new.ErrorFlag==0","goff");
    TH2F *bpm4aX_main2D = (TH2F*)gROOT->FindObject("bpm4aX_main2D");
    mul->Draw("new.yield_bpm4aX:Entry$>>bpm4aX_new2D()","ErrorFlag!=0 && new.ErrorFlag==0","goff");
    TH2F *bpm4aX_new2D = (TH2F*)gROOT->FindObject("bpm4aX_new2D");
    bpm4aX_new2D->SetMarkerStyle(7);
    bpm4aX_new2D->SetMarkerColor(2);
    bpm4aX_main2D->Draw();
    bpm4aX_new2D->Draw("same");
    c5->SaveAs("./ComparisonOutputs/Plots/temp/bpm4aX2D.pdf");

    TCanvas *c6 = new TCanvas();
    c6->Divide(2,1);
    mul->Draw("new.yield_bpm4aX>>bpm4aX_main()","new.ErrorFlag==0","goff");
    TH2F *bpm4aX_main = (TH2F*)gROOT->FindObject("bpm4aX_main");
    mul->Draw("new.yield_bpm4aX>>bpm4aX_new()","ErrorFlag!=0 && new.ErrorFlag==0","goff");
    TH2F *bpm4aX_new = (TH2F*)gROOT->FindObject("bpm4aX_new");
    c6->cd(1);
    bpm4aX_main->Draw();
    c6->cd(2);
    bpm4aX_new->Draw();
    c6->SaveAs("./ComparisonOutputs/Plots/temp/bpm4aX.pdf");


    TCanvas *c7 = new TCanvas();
    mul->Draw("new.yield_bpm4aY:Entry$>>bpm4aY_main2D()","new.ErrorFlag==0","goff");
    TH2F *bpm4aY_main2D = (TH2F*)gROOT->FindObject("bpm4aY_main2D");
    mul->Draw("new.yield_bpm4aY:Entry$>>bpm4aY_new2D()","ErrorFlag!=0 && new.ErrorFlag==0","goff");
    TH2F *bpm4aY_new2D = (TH2F*)gROOT->FindObject("bpm4aY_new2D");
    bpm4aY_new2D->SetMarkerStyle(7);
    bpm4aY_new2D->SetMarkerColor(2);
    bpm4aY_main2D->Draw();
    bpm4aY_new2D->Draw("same");
    c7->SaveAs("./ComparisonOutputs/Plots/temp/bpm4aY2D.pdf");

    TCanvas *c8 = new TCanvas();
    c8->Divide(2,1);
    mul->Draw("new.yield_bpm4aY>>bpm4aY_main()","new.ErrorFlag==0","goff");
    TH2F *bpm4aY_main = (TH2F*)gROOT->FindObject("bpm4aY_main");
    mul->Draw("new.yield_bpm4aY>>bpm4aY_new()","ErrorFlag!=0 && new.ErrorFlag==0","goff");
    TH2F *bpm4aY_new = (TH2F*)gROOT->FindObject("bpm4aY_new");
    c8->cd(1);
    bpm4aY_main->Draw();
    c8->cd(2);
    bpm4aY_new->Draw();
    c8->SaveAs("./ComparisonOutputs/Plots/temp/bpm4aY.pdf");


    TCanvas *c9 = new TCanvas();
    mul->Draw("new.yield_bpm12X:Entry$>>bpm12X_main2D()","new.ErrorFlag==0","goff");
    TH2F *bpm12X_main2D = (TH2F*)gROOT->FindObject("bpm12X_main2D");
    mul->Draw("new.yield_bpm12X:Entry$>>bpm12X_new2D()","ErrorFlag!=0 && new.ErrorFlag==0","goff");
    TH2F *bpm12X_new2D = (TH2F*)gROOT->FindObject("bpm12X_new2D");
    bpm12X_new2D->SetMarkerStyle(7);
    bpm12X_new2D->SetMarkerColor(2);
    bpm12X_main2D->Draw();
    bpm12X_new2D->Draw("same");
    c9->SaveAs("./ComparisonOutputs/Plots/temp/bpm12X2D.pdf");

    TCanvas *c10 = new TCanvas();
    c10->Divide(2,1);
    mul->Draw("new.yield_bpm12X>>bpm12X_main()","new.ErrorFlag==0","goff");
    TH2F *bpm12X_main = (TH2F*)gROOT->FindObject("bpm12X_main");
    mul->Draw("new.yield_bpm12X>>bpm12X_new()","ErrorFlag!=0 && new.ErrorFlag==0","goff");
    TH2F *bpm12X_new = (TH2F*)gROOT->FindObject("bpm12X_new");
    c10->cd(1);
    bpm12X_main->Draw();
    c10->cd(2);
    bpm12X_new->Draw();
    c10->SaveAs("./ComparisonOutputs/Plots/temp/bpm12X.pdf");

    gSystem->Exec(Form("pdfunite ./ComparisonOutputs/Plots/temp/* ./ComparisonOutputs/Plots/summary%d.pdf", Run));
    
    return;
}
