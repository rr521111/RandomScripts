TString filepath1 = "/chafs2/work1/apar/japanOutput/";
TString filepath2 = "/volatile/halla/parity/crex-respin1/japanOutput/";

void RunErrorComparison(Int_t Run){
    TFile* file1 = TFile::Open(Form("%s/prexPrompt_pass2_%d.000.root", filepath1.Data(), Run), "READ");
    TFile* file2 = TFile::Open(Form("%s/prexPrompt_pass2_%d.000.root", filepath2.Data(), Run), "READ");

    file1->cd();

    TTree* mul = (TTree*)gROOT->FindObject("mul");
    TTree* mulc_lrb_burst = (TTree*)gROOT->FindObject("mulc_lrb_burst");

    mul->AddFriend("new = mul",file2);
    mul->AddFriend("old_mulc_lrb_burst = mulc_lrb_burst",file1);
    mul->AddFriend("new_mulc_lrb_burst = mulc_lrb_burst",file2);

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

    TCanvas *c1a = new TCanvas();
    mul->Draw("new.yield_usl:Entry$>>yusl_main2D()","new.ErrorFlag==0","goff");
    TH2F *yusl_main2D = (TH2F*)gROOT->FindObject("yusl_main2D");
    mul->Draw("new.yield_usl:Entry$>>yusl_new2D()","ErrorFlag!=0 && new.ErrorFlag==0","goff");
    TH2F *yusl_new2D = (TH2F*)gROOT->FindObject("yusl_new2D");
    yusl_new2D->SetMarkerStyle(7);
    yusl_new2D->SetMarkerColor(2);
    yusl_main2D->Draw();
    yusl_new2D->Draw("same");
    c1a->SaveAs("./ComparisonOutputs/Plots/temp/yusl2D.pdf");

    TCanvas *c2a = new TCanvas();
    c2a->Divide(2,1);
    mul->Draw("new.yield_usl>>yusl_main()","new.ErrorFlag==0","goff");
    TH2F *yusl_main = (TH2F*)gROOT->FindObject("yusl_main");
    mul->Draw("new.yield_usl>>yusl_new()","ErrorFlag!=0 && new.ErrorFlag==0","goff");
    TH2F *yusl_new = (TH2F*)gROOT->FindObject("yusl_new");
    c2a->cd(1);
    yusl_main->Draw();
    c2a->cd(2);
    yusl_new->Draw();
    c2a->SaveAs("./ComparisonOutputs/Plots/temp/yusl.pdf");

    TCanvas *c1b = new TCanvas();
    mul->Draw("new.asym_usr:Entry$>>usr_main2D()","new.ErrorFlag==0","goff");
    TH2F *usr_main2D = (TH2F*)gROOT->FindObject("usr_main2D");
    mul->Draw("new.asym_usr:Entry$>>usr_new2D()","ErrorFlag!=0 && new.ErrorFlag==0","goff");
    TH2F *usr_new2D = (TH2F*)gROOT->FindObject("usr_new2D");
    usr_new2D->SetMarkerStyle(7);
    usr_new2D->SetMarkerColor(2);
    usr_main2D->Draw();
    usr_new2D->Draw("same");
    c1b->SaveAs("./ComparisonOutputs/Plots/temp/usr2D.pdf");

    TCanvas *c2b = new TCanvas();
    c2b->Divide(2,1);
    mul->Draw("new.asym_usr>>usr_main()","new.ErrorFlag==0","goff");
    TH2F *usr_main = (TH2F*)gROOT->FindObject("usr_main");
    mul->Draw("new.asym_usr>>usr_new()","ErrorFlag!=0 && new.ErrorFlag==0","goff");
    TH2F *usr_new = (TH2F*)gROOT->FindObject("usr_new");
    c2b->cd(1);
    usr_main->Draw();
    c2b->cd(2);
    usr_new->Draw();
    c2b->SaveAs("./ComparisonOutputs/Plots/temp/usr.pdf");

    TCanvas *c1ba = new TCanvas();
    mul->Draw("new.yield_usr:Entry$>>yusr_main2D()","new.ErrorFlag==0","goff");
    TH2F *yusr_main2D = (TH2F*)gROOT->FindObject("yusr_main2D");
    mul->Draw("new.yield_usr:Entry$>>yusr_new2D()","ErrorFlag!=0 && new.ErrorFlag==0","goff");
    TH2F *yusr_new2D = (TH2F*)gROOT->FindObject("yusr_new2D");
    yusr_new2D->SetMarkerStyle(7);
    yusr_new2D->SetMarkerColor(2);
    yusr_main2D->Draw();
    yusr_new2D->Draw("same");
    c1ba->SaveAs("./ComparisonOutputs/Plots/temp/yusr2D.pdf");

    TCanvas *c2ba = new TCanvas();
    c2ba->Divide(2,1);
    mul->Draw("new.yield_usr>>yusr_main()","new.ErrorFlag==0","goff");
    TH2F *yusr_main = (TH2F*)gROOT->FindObject("yusr_main");
    mul->Draw("new.yield_usr>>yusr_new()","ErrorFlag!=0 && new.ErrorFlag==0","goff");
    TH2F *yusr_new = (TH2F*)gROOT->FindObject("yusr_new");
    c2ba->cd(1);
    yusr_main->Draw();
    c2ba->cd(2);
    yusr_new->Draw();
    c2ba->SaveAs("./ComparisonOutputs/Plots/temp/yusr.pdf");

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

    TCanvas *c3a = new TCanvas();
    mul->Draw("new.yield_bcm_an_us:Entry$>>ybcm_main2D()","new.ErrorFlag==0","goff");
    TH2F *ybcm_main2D = (TH2F*)gROOT->FindObject("ybcm_main2D");
    mul->Draw("new.yield_bcm_an_us:Entry$>>ybcm_new2D()","ErrorFlag!=0 && new.ErrorFlag==0","goff");
    TH2F *ybcm_new2D = (TH2F*)gROOT->FindObject("ybcm_new2D");
    ybcm_new2D->SetMarkerStyle(7);
    ybcm_new2D->SetMarkerColor(2);
    ybcm_main2D->Draw();
    ybcm_new2D->Draw("same");
    c3a->SaveAs("./ComparisonOutputs/Plots/temp/ybcm2D.pdf");

    TCanvas *c4a = new TCanvas();
    c4a->Divide(2,1);
    mul->Draw("new.yield_bcm_an_us>>ybcm_main()","new.ErrorFlag==0","goff");
    TH2F *ybcm_main = (TH2F*)gROOT->FindObject("ybcm_main");
    mul->Draw("new.yield_bcm_an_us>>ybcm_new()","ErrorFlag!=0 && new.ErrorFlag==0","goff");
    TH2F *ybcm_new = (TH2F*)gROOT->FindObject("ybcm_new");
    c4a->cd(1);
    ybcm_main->Draw();
    c4a->cd(2);
    ybcm_new->Draw();
    c4a->SaveAs("./ComparisonOutputs/Plots/temp/ybcm.pdf");

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

    TCanvas *c11 = new TCanvas();
    mul->Draw("new_mulc_lrb_burst.cor_asym_us_avg:Entry$>>cor_asym_us_avg_main2D()","new.ErrorFlag==0","goff");
    TH2F *cor_asym_us_avg_main2D = (TH2F*)gROOT->FindObject("cor_asym_us_avg_main2D");
    mul->Draw("new_mulc_lrb_burst.cor_asym_us_avg:Entry$>>cor_asym_us_avg_new2D()","ErrorFlag!=0 && new.ErrorFlag==0","goff");
    TH2F *cor_asym_us_avg_new2D = (TH2F*)gROOT->FindObject("cor_asym_us_avg_new2D");
    cor_asym_us_avg_new2D->SetMarkerStyle(7);
    cor_asym_us_avg_new2D->SetMarkerColor(2);
    cor_asym_us_avg_main2D->Draw();
    cor_asym_us_avg_new2D->Draw("same");
    c11->SaveAs("./ComparisonOutputs/Plots/temp/cor_asym_us_avg2D.pdf");

    TCanvas *c12 = new TCanvas();
    c12->Divide(2,1);
    mul->Draw("new_mulc_lrb_burst.cor_asym_us_avg>>cor_asym_us_avg_main()","new.ErrorFlag==0","goff");
    TH2F *cor_asym_us_avg_main = (TH2F*)gROOT->FindObject("cor_asym_us_avg_main");
    mul->Draw("new_mulc_lrb_burst.cor_asym_us_avg>>cor_asym_us_avg_new()","ErrorFlag!=0 && new.ErrorFlag==0","goff");
    TH2F *cor_asym_us_avg_new = (TH2F*)gROOT->FindObject("cor_asym_us_avg_new");
    c12->cd(1);
    cor_asym_us_avg_main->Draw();
    c12->cd(2);
    cor_asym_us_avg_new->Draw();
    c12->SaveAs("./ComparisonOutputs/Plots/temp/cor_asym_us_avg.pdf");

    TCanvas *c13 = new TCanvas();
    mul->Draw("new_mulc_lrb_burst.cor_asym_us_dd:Entry$>>cor_asym_us_dd_main2D()","new.ErrorFlag==0","goff");
    TH2F *cor_asym_us_dd_main2D = (TH2F*)gROOT->FindObject("cor_asym_us_dd_main2D");
    mul->Draw("new_mulc_lrb_burst.cor_asym_us_dd:Entry$>>cor_asym_us_dd_new2D()","ErrorFlag!=0 && new.ErrorFlag==0","goff");
    TH2F *cor_asym_us_dd_new2D = (TH2F*)gROOT->FindObject("cor_asym_us_dd_new2D");
    cor_asym_us_dd_new2D->SetMarkerStyle(7);
    cor_asym_us_dd_new2D->SetMarkerColor(2);
    cor_asym_us_dd_main2D->Draw();
    cor_asym_us_dd_new2D->Draw("same");
    c13->SaveAs("./ComparisonOutputs/Plots/temp/cor_asym_us_dd2D.pdf");

    TCanvas *c14 = new TCanvas();
    c14->Divide(2,1);
    mul->Draw("new_mulc_lrb_burst.cor_asym_us_dd>>cor_asym_us_dd_main()","new.ErrorFlag==0","goff");
    TH2F *cor_asym_us_dd_main = (TH2F*)gROOT->FindObject("cor_asym_us_dd_main");
    mul->Draw("new_mulc_lrb_burst.cor_asym_us_dd>>cor_asym_us_dd_new()","ErrorFlag!=0 && new.ErrorFlag==0","goff");
    TH2F *cor_asym_us_dd_new = (TH2F*)gROOT->FindObject("cor_asym_us_dd_new");
    c14->cd(1);
    cor_asym_us_dd_main->Draw();
    c14->cd(2);
    cor_asym_us_dd_new->Draw();
    c14->SaveAs("./ComparisonOutputs/Plots/temp/cor_asym_us_dd.pdf");

    gSystem->Exec(Form("pdfunite ./ComparisonOutputs/Plots/temp/* ./ComparisonOutputs/Plots/summary%d.pdf", Run));
    
    return;
}
