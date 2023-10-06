#include <iostream>
#include <fstream>
#include <vector>
#include <TFile.h>
#include <TNtuple.h>
#include <TGraph.h>
#include <TF1.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLeaf.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>
#include <numeric>
#include <stdio.h>
#include "TString.h"
#include "TH1D.h"
#include "TMultiGraph.h"
#include "TStyle.h"
#include <chrono>
#include <ctime>
using namespace std;

// print file names vector
void print(std::vector<string> const &input)
{
    for (int i = 0; i < input.size(); i++)
        std::cout << input.at(i) << ' ' << endl;
}
void print_double(std::vector<double> const &input)
{
    for (int i = 0; i < input.size(); i++)
        std::cout << i << " " << input.at(i) << ' ' << endl;
}

vector<double> MovingAverage(const vector<double> &signal, int windowSize)
{
    int n = signal.size();
    vector<double> smoothedSignal(n);
    for (int i = 0; i < n; i++)
    {
        int start = max(0, i - windowSize / 2);
        int end = min(n - 1, i + windowSize / 2);
        int size = end - start + 1;
        vector<double> window(signal.begin() + start, signal.begin() + end + 1);
        smoothedSignal[i] = accumulate(window.begin(), window.end(), 0.0) / size;
    }
    return smoothedSignal;
}

void Calibration()
{
    // current date/time based on current system
    time_t now = time(0);
    // convert now to string form
    char *dt = ctime(&now);
    TString s = dt;
    TString pdf_wf = "test";

    cout << "The local date and time is: " << s << endl;
    // convert now to tm struct for UTC
    tm *gmtm = gmtime(&now);
    dt = asctime(gmtm);
    cout << "The UTC date and time is:" << dt << endl;

    ofstream myfile;
    myfile.open(s + ".txt");

    vector<vector<int>> map_down{{97, 105, 98, 106, 99, 107, 100}, // for horiz mouns
                                 {108, 101, 109, 102, 110, 103, 111},
                                 {104, 112, 65, 73, 66, 74, 67},
                                 {75, 68, 76, 69, 77, 70, 78},
                                 {71, 79, 72, 80, 81, 89, 82},
                                 {90, 83, 91, 84, 92, 85, 93},
                                 {86, 94, 87, 95, 88, 96, 113},
                                 {121, 114, 122, 115, 123, 116, 124},
                                 {117, 125, 118, 126, 119, 127, 120}};

    vector<vector<int>> map_up{{33, 41, 34, 42, 35, 43, 36},
                               {44, 37, 45, 38, 46, 39, 47},
                               {40, 48, 1, 9, 2, 10, 3},
                               {11, 4, 12, 5, 13, 6, 14},
                               {7, 15, 8, 16, 17, 25, 18},
                               {26, 19, 27, 20, 28, 21, 29},
                               {22, 30, 23, 31, 24, 32, 49},
                               {57, 50, 58, 51, 59, 52, 60},
                               {53, 61, 54, 62, 55, 63, 56}};

    // map for data.root
    /*vector<vector<int>> map{{11, 4, 12, 5, 13, 6, 14},           // 4
                            {7, 15, 8, 16, 17, 25, 18},          // 5
                            {26, 19, 27, 20, 28, 21, 29},        // 6
                            {22, 30, 23, 31, 24, 32, 49},        // 7
                            {57, 50, 58, 51, 59, 52, 60},        // 8
                            {53, 61, 54, 62, 55, 63, 56},        // 9
                            {97, 105, 98, 106, 99, 107, 100},    // 1s
                            {108, 101, 109, 102, 110, 103, 111}, // 2s
                            {104, 112, 65, 73, 66, 74, 67}};     // 3s*/

    // new mapping for 27_02_23 cosmic
    vector<vector<int>> map{{11, 4, 12, 5, 13, 6, 14},            // 4
                            {7, 15, 8, 16, 17, 25, 18},           // 5
                            {33, 41, 34, 42, 35, 43, 36},         // 6
                            {44, 37, 45, 38, 46, 39, 47},         // 7
                            {57, 50, 58, 51, 59, 52, 60},         // 8
                            {53, 61, 54, 62, 55, 63, 56},         // 9
                            {97, 105, 98, 106, 99, 107, 100},     // 1s
                            {108, 101, 109, 102, 110, 103, 111},  // 2s
                            {104, 112, 65, 73, 66, 74, 67},       // 3s
                            {75, 68, 76, 69, 77, 70, 78},         // 4s
                            {71, 79, 72, 80, 81, 89, 82},         // 5s
                            {90, 83, 91, 84, 92, 85, 93},         // 6s
                            {86, 94, 87, 95, 88, 96, 113},        // 7s
                            {121, 114, 122, 115, 123, 116, 124},  // 8s
                            {117, 125, 118, 126, 119, 127, 120}}; // 9s

    int cells_arr_for_loop[128] = {0};
    int iiii = 0;
    for (int i = 0; i < 15; i++)
        for (int j = 0; j < 7; j++)
        {
            cells_arr_for_loop[iiii] = map[i][j];
            iiii++;
            // cout << map[i][j] << " ";
        }
    // cout << endl;
    std::sort(std::begin(cells_arr_for_loop), std::end(cells_arr_for_loop));
    // for (int i = 0; i < 63; i++)
    // cout << cells_arr_for_loop[i] << " ";

    vector<vector<int>> mapxyz{{1, 2},
                               {2, 2},
                               {3, 2},
                               {1, 3},
                               {2, 3},
                               {3, 3},
                               {1, 4},
                               {2, 4},
                               {3, 4},
                               {1, 5},
                               {2, 5},
                               {3, 5},
                               {1, 6},
                               {2, 6},
                               {3, 6}};

    vector<int>
        temp,
        temp2;
    vector<double> signal, zero_level_vector;
    int mods5[5] = {0};
    // TFile *_file0 = TFile::Open("/mnt/c/Work/root/builddir/macros/Calibration/30_08_cosmo_large_2.root");
    TFile *_file0 = TFile::Open("/mnt/c/Work/root/builddir/macros/Calibration/cosmic_27_02_23_large.root");
    TTree *adc64_data = (TTree *)_file0->Get("adc64_data");
    Double_t x[2048] = {0.};
    Double_t y[2048] = {0.};
    vector<string> branch_names;
    vector<string> three_cells;
    string xxx = "";
    TH1F *zl_distr = new TH1F("zl_distr", "zl_distr", 100, 30000, 32000);
    TH1F *sss = new TH1F("sss", "sss", 100, 0, 5000);
    TH2F *zl_for_fit = new TH2F("zl_for_fit", "zl_for_fit", 50, 0, 50, 10000, 30000, 32000);
    TH1F *Section[128];
    for (int i = 0; i < 128; i++)
        Section[i] = new TH1F(TString::Format("amp%d", i), TString::Format("amp%d", i), 100, 0., 10000.);
    bool flag = false;
    const int MAX_N_SAMPLES = 2048;
    Short_t wf[MAX_N_SAMPLES];
    int how_many_cells = 0;

    // for (int cell_iter = 0; cell_iter < 63; cell_iter++) //usual
    for (int cell_iter = 0; cell_iter < 128; cell_iter++) // horiz
        if (cells_arr_for_loop[cell_iter] != 0)
            branch_names.push_back("channel_" + std::to_string(cells_arr_for_loop[cell_iter]));
    // print(branch_names);       // check
    TGraph *g1 = new TGraph(); // using the blank constructor
    TGraph *g2 = new TGraph(); // using the blank constructor

    TCanvas *cfs = new TCanvas();
    cfs->Print("/mnt/c/Work/root/builddir/macros/Calibration/pics/" + pdf_wf + ".pdf(", "pdf");
    // cfs->Print("cf_3cells_strict.pdf(", "pdf");
    cout << adc64_data->GetEntries() << endl;
    int n_events_to_loop = 1500000;                                               // adc64_data->GetEntries();       // adc64_data->GetEntries()
    for (int entry = n_events_to_loop * 0; entry < n_events_to_loop * 1; entry++) // adc64_data->GetEntries()
    {
        zl_for_fit->Reset();
        memset(mods5, 0, sizeof mods5);

        string names_fired = "";
        int i_array = 0;
        how_many_cells = 0;
        if (entry % 100000 == 0)
            cout << "Event: " << entry << "/" << n_events_to_loop << endl;
        adc64_data->GetEntry(entry);
        for (auto i_name = branch_names.begin(); i_name != branch_names.end(); i_name++)
        {
            TString name = *i_name;
            TBranch *channel = (TBranch *)adc64_data->GetBranch(name);
            if (name == "channel_22" || name == "channel_29")
                continue;
            // cout << entry << " " << name << endl;
            double ss = channel->GetLeaf("wf_size")->GetValue();
            if (ss != 0)
            {
                string s2(name(8, 10));
                temp2.push_back(stoi(s2));
                how_many_cells++;
            }
        }

        int count_mods = 0;
        if (how_many_cells == 5)
            for (int i = 0; i < 15; i++)
                for (int j = 0; j < 7; j++)
                    for (int it = 0; it < temp2.size(); it++)
                        if (map[i][j] == temp2.at(it))
                        {
                            // cout << entry << " mod " << i << " " << map[i][j] << " coincide " << temp2.at(it) << endl;
                            mods5[count_mods] = i;
                            count_mods++;
                        }

        temp2.clear();
        bool check = true;
        for (int i = 0; i < 5; i++)
            for (int j = 0; j < 5; j++)
            {
                if (mods5[i] == mods5[j] && i != j)
                {
                    // cout << entry << " mods check " << mods5[i] << " " << mods5[j] << " i " << i << " j " << j << endl;
                    check = false;
                    break;
                }
                if (mapxyz[mods5[i]][0] != mapxyz[mods5[j]][0] && i != j)
                {
                    // cout << entry << " map check | mod: " << mods5[i] << " " << mapxyz[mods5[i]][0] << " " << mapxyz[mods5[j]][0] << endl;
                    check = false;
                    break;
                }
            }

        if (how_many_cells == 5 && check)
        {
            // for (int m = 0; m < 5; m++)
            // cout << entry << " mod " << mods5[m] << " x " << mapxyz[mods5[m]][0] << endl;
            for (auto i_name = branch_names.begin(); i_name != branch_names.end(); i_name++)
            {
                int amp = 400000, gate_beg = 80, gate_end = 120;
                TString name = *i_name;
                if (name == "channel_22" || name == "channel_29")
                    continue;
                TBranch *channel = (TBranch *)adc64_data->GetBranch(name);
                double ss = channel->GetLeaf("wf_size")->GetValue();
                // if (name == "channel_52")
                if (ss != 0)
                {
                    double zl_bf = 0;
                    /*for (Int_t i = 0; i < 200; i++)
                    {
                        g2->SetPoint(g2->GetN(), i, channel->GetLeaf("wf")->GetValue(i));
                    }*/
                    for (Int_t i = 0; i < 80; i++)
                    {
                        // g2->SetMarkerColor(kRed);
                        // g2->SetPoint(g2->GetN(), i, channel->GetLeaf("wf")->GetValue(i));
                        double xx = channel->GetLeaf("wf")->GetValue(i);
                        zl_bf += xx / 80;
                        // zl_for_fit->Fill(i, fabs(xx));
                        // zero_level_vector.push_back(channel->GetLeaf("wf")->GetValue(i));
                    }
                    // cout << " bf " << zl << endl;
                    /*double zl_af = 0;
                    for (Int_t i = 100; i <= 199; i++)
                    {
                        // g2->SetMarkerColor(kBlue);
                        // g2->SetPoint(g2->GetN(), i, channel->GetLeaf("wf")->GetValue(i));
                        double xx = channel->GetLeaf("wf")->GetValue(i);
                        zl_af += xx / 100;
                        // zl_for_fit->Fill(i, fabs(xx));
                    }*/
                    double zl = (zl_bf); // + zl_af) / 2;

                    // cout << " af " << zl << endl;

                    /*TCanvas *cf = new TCanvas();
                    zl_for_fit->Draw();
                    zl_for_fit->Fit("pol0", "Q");
                    TF1 *g = (TF1 *)zl_for_fit->GetListOfFunctions()->FindObject("pol0");
                    double zlf = g->GetParameter(0);
                    cout << "f " << zlf << " nf " << fabs(zl) << " " << zlf - fabs(zl) << endl;
                    // zl_distr->Fill(zlf);
                    cf->Print("cf_3cells_strict.pdf", "pdf");
                    zl_for_fit->Reset();*/
                    // double zl = 0;
                    signal.clear();
                    for (Int_t i = gate_beg; i <= gate_end; i++)
                    {
                        // signal.push_back(channel->GetLeaf("wf")->GetValue(i));
                        //  g1->SetPoint(g1->GetN(), i, channel->GetLeaf("wf")->GetValue(i));
                        double xx = channel->GetLeaf("wf")->GetValue(i);
                        // cout << "xx " << xx << " ss " << ss <<endl;
                        if (xx < amp)
                            amp = xx;
                    }
                    // TCanvas *cpeak = new TCanvas();
                    /*g1->SetMarkerStyle(8);
                    g1->Draw("");
                    g1->SetTitle(TString::Format("Peak %i | ", amp) + TString::Format(" ZL %f | ", zl) + TString::Format(" Ampl %f | ", amp - zl));
                    cpeak->Print("/mnt/c/Work/root/builddir/macros/Calibration/pics/" + pdf_wf + ".pdf", "pdf");

                    g1->Set(0);
    */

                    // g2->GetYaxis()->SetRangeUser(zl - 500, amp + 100);
                    // g2->SetTitle(TString::Format("ZL %f", zl));
                    //  gPad->Print(name + TString::Format("_entry%i", entry) + ".png");
                    // cZL->Print("/mnt/c/Work/root/builddir/macros/Calibration/pics/" + pdf_wf + ".pdf", "pdf");

                    int windowSize = 5;
                    // vector<double> smoothedSignal = MovingAverage(zero_level_vector, windowSize);
                    //  cout << entry << " " << name << " amp " << amp << " MAX smoothed ";
                    /*cout << " smooth " << endl;

                    print_double(smoothedSignal);
                    cout << " orig " << endl;
                    print_double(zero_level_vector);*/
                    // auto max = max_element(std::begin(smoothedSignal), std::end(smoothedSignal));
                    // double average = accumulate(smoothedSignal.begin(), smoothedSignal.end(), 0.0) / smoothedSignal.size();
                    //  cout << "The average is" << average << endl;
                    //   cout << *max << " " << amp - *max << endl;
                    // zero_level_vector.clear();
                    //  cout << entry << " " << name << " amp " << amp << endl;
                    //  if (fabs(amp) > 1000)
                    //  Section[i_array]->Fill(amp - (zl_af + zl_bf) / 2);
                    // zl_distr->Fill(fabs(average));

                    Section[i_array]->Fill(zl - amp);
                    // cout << entry << " " << zl - amp << " " << amp << " " << zl << endl;
                    // if (zl - amp < 1000)
                    /*{
                        TCanvas *cZL = new TCanvas();
                        g2->SetMarkerStyle(8);
                        g2->Draw("");
                        g2->SetTitle(name + TString::Format(" Peak %i | ", amp) + TString::Format(" ZL %f | ", zl) + TString::Format(" Ampl %f | ", zl - amp));
                        cZL->Print("/mnt/c/Work/root/builddir/macros/Calibration/pics/" + pdf_wf + ".pdf", "pdf");
                    }
                    g2->Set(0);*/

                    // cout << amp - (zl_af + zl_bf) / 2 << " " <<  amp << " " << endl;
                }
                i_array++;
                flag = false;
            }
        }
    }
    TCanvas *czl_check = new TCanvas();
    zl_distr->Draw("hist");
    czl_check->Print("/mnt/c/Work/root/builddir/macros/Calibration/pics/" + pdf_wf + ".pdf", "pdf");

    cfs->Print("/mnt/c/Work/root/builddir/macros/Calibration/pics/" + pdf_wf + ".pdf)", "pdf");

    cout << "wf saved!" << endl;
    ////////Draw///////////
    TCanvas *c1 = new TCanvas();
    c1->Print(s + ".pdf(", "pdf");
    c1->DivideSquare(128);
    int i_array = 0;
    TFile histoFileFull("calib_5cells.root", "RECREATE");

    for (auto i_name = branch_names.begin(); i_name != branch_names.end(); i_name++)
    {
        TString name = *i_name;
        cout << name << endl;
        c1->cd(i_array + 1);
        double lvl = 0;
        int peak, peak_i, left, right = 0;
        for (int i = 11; i < 100; i++)
            if (Section[i_array]->GetBinContent(i) > lvl)
            {
                lvl = Section[i_array]->GetBinContent(i);
                peak = i * 100;
                peak_i = i;
            }
        for (int i = 1; i < peak_i; i++)
            if (Section[i_array]->GetBinContent(i) > lvl / 3)
            {
                left = i * 100;
                break;
            }

        for (int i = peak_i; i < 45; i++)
            if (Section[i_array]->GetBinContent(i) < lvl / 3)
            {
                right = i * 100;
                break;
            }
        TF1 *fit1 = new TF1("fit1", "gaus", left, right);
        cout << name << " L " << left << " R " << right << " MAX " << peak << endl;
        Section[i_array]->Draw("hist");
        Section[i_array]->SetTitle(name);
        Section[i_array]->Fit(fit1, "R");
        fit1->Draw("same");

        TCanvas *c2 = new TCanvas();
        Section[i_array]->Draw("hist");
        Section[i_array]->SetTitle(name);
        Section[i_array]->GetXaxis()->SetTitle("Ampl [channels]");
        Section[i_array]->GetYaxis()->SetTitle("counts");
        Section[i_array]->Fit(fit1, "QR");
        fit1->Draw("same");
        if (fit1->GetParError(1) / fit1->GetParameter(1) < 0.05 && fit1->GetParameter(1) > 0)
            myfile << name << " " << fit1->GetParameter(1) << " " << fit1->GetParameter(2) << endl;
        else
            myfile << name << " " << 0 << " " << 0 << endl;
        sss->Fill(fit1->GetParameter(2));
        c2->Print(s + ".pdf", "pdf");
        Section[i_array]->Write();

        i_array++;
    }
    c1->Print(s + ".pdf)", "pdf");
    TCanvas *cf = new TCanvas();
    sss->Draw("hist");
    myfile.close();
}