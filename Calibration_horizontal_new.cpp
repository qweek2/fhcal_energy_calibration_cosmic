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
#include <TH3F.h>
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
#include <stdlib.h>

// print file names vector
void print(std::vector<string> const &input)
{
    for (int i = 0; i < input.size(); i++)
        std::cout << i << " | " << input.at(i) << endl;
}

void print_int(std::vector<int> const &input)
{
    for (int i = 0; i < input.size(); i++)
        std::cout << i << " | " << input.at(i) << endl;
}

int **vectorToArray(vector<vector<int>> &vec)
{
    // Get the dimensions of the 2D vector
    int numRows = vec.size();
    int numCols = vec[0].size();

    // Create a 2D array of integers with the same dimensions as the 2D vector
    int **arr = new int *[numRows];
    for (int i = 0; i < numRows; i++)
    {
        arr[i] = new int[numCols];
    }

    // Copy the elements of the 2D vector to the 2D array
    for (int i = 0; i < numRows; i++)
    {
        for (int j = 0; j < numCols; j++)
        {
            arr[i][j] = vec[i][j];
        }
    }

    return arr;
}

void calculateStdDev(const array<double, 10> &data, double &stdev, double &cov)
{
    double sum = 0.0, mean = 0.0, variance = 0.0;

    // Calculate mean
    for (int i = 0; i < data.size(); ++i)
    {
        sum += data[i];
    }
    mean = sum / data.size();

    // Calculate variance
    for (int i = 0; i < data.size(); ++i)
    {
        variance += pow(data[i] - mean, 2);
    }
    variance /= data.size();

    // Calculate standard deviation
    stdev = sqrt(variance);

    // Calculate coefficient of variation
    cov = stdev / mean;
}

bool find_subseq(int ev, int C[], int lengthC, int D[], int lengthD, vector<string> &subsequence)
{
    // Define input arrays C and D

    // Get the lengths of input arrays C and D
    // int lengthA = sizeof(C) / sizeof(C[0]);
    // int lengthB = sizeof(D) / sizeof(D[0]);

    // Initialize variables
    int maxLength = 0;
    int endIndex = -1;

    // Create a vector to store the current subsequence of D that occurs in C in a row
    vector<int> currentSubsequence;

    // Find the longest subsequence of numbers from D (not necessarily using all elements of D) that occurs in C in a row
    for (int i = 0; i < lengthC; i++)
    {
        if (C[i] == D[currentSubsequence.size()])
        {
            currentSubsequence.push_back(C[i]);
            if (currentSubsequence.size() > maxLength)
            {
                maxLength = currentSubsequence.size();
                endIndex = i;
            }
        }
        else if (!currentSubsequence.empty() && C[i] == currentSubsequence[0])
        {
            currentSubsequence.erase(currentSubsequence.begin());
            currentSubsequence.push_back(C[i]);
        }
        else
        {
            currentSubsequence.clear();
        }
    }

    // Output the results
    if (endIndex != -1)
    {
        if (maxLength > 6)
        {
            // cout << ev << " | Longest subsequence of numbers found starting at index " << endIndex - maxLength + 1 << " with length " << maxLength << endl;
            // cout << "The subsequence is: ";
            for (int i = endIndex - maxLength + 1; i <= endIndex; i++)
            {
                // cout << C[i] << " ";
                subsequence.push_back("channel_" + std::to_string(C[i]));
            }
            // cout << endl;
        }
    }
    else
    {
        // cout << "Longest subsequence of numbers from D (not necessarily using all elements of D) that occurs in C in a row not found" << endl;
    }
    if (maxLength > 4)
        return true;
    else
        return false;
}

void Calibration_horizontal_new()
{
    ofstream myfile;
    myfile.open("horiz.txt");
    // map for 30_08 files
    /*vector<vector<int>> map{{32, 40, 33, 41, 34, 42, 35},  // 1
                            {43, 36, 44, 37, 45, 38, 46},  // 2
                            {39, 47, 0, 8, 1, 9, 2},       // 3
                            {10, 3, 11, 4, 12, 5, 13},     // 4
                            {6, 14, 7, 15, 16, 24, 17},    // 5
                            {25, 18, 26, 19, 27, 20, 28},  // 6
                            {21, 29, 22, 30, 23, 31, 48},  // 7
                            {56, 49, 57, 50, 58, 51, 59},  // 8
                            {52, 60, 53, 61, 54, 62, 55}}; // 9
                            */
    // map for data.root (new data collected by Strizhak in a stack config of modules)
    /*vector<vector<int>> map{{11, 4, 12, 5, 13, 6, 14},           // 1
                            {7, 15, 8, 16, 17, 25, 18},          // 2
                            {26, 19, 27, 20, 28, 21, 29},        // 3
                            {22, 30, 23, 31, 24, 32, 49},        // 4
                            {57, 50, 58, 51, 59, 52, 60},        // 5
                            {53, 61, 54, 62, 55, 63, 56},        // 6
                            {97, 105, 98, 106, 99, 107, 100},    // 7
                            {108, 101, 109, 102, 110, 103, 111}, // 8
                            {104, 112, 65, 73, 66, 74, 67}};     // 9

    vector<vector<int>> mapxyz{{0, 0},  // 1
                               {0, 1},  // 2
                               {0, 2},  // 3
                               {1, 0},  // 4
                               {1, 1},  // 5
                               {1, 2},  // 6
                               {2, 0},  // 7
                               {2, 1},  // 8
                               {2, 2}}; // 9*/

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

    int numRows = map.size();
    int numCols = map[0].size();
    int cells_arr_for_loop[numCols * numRows];
    memset(cells_arr_for_loop, -1, sizeof cells_arr_for_loop);
    int iiii = 0;

    for (int i = 0; i < numRows; i++)
        for (int j = 0; j < numCols; j++)
        {
            cells_arr_for_loop[iiii] = map[i][j];
            iiii++;
            // cout << map[i][j] << " ";
        }
    // std::sort(std::begin(cells_arr_for_loop), std::end(cells_arr_for_loop));

    int row_fired[numCols];
    int row_from_map[numCols];

    vector<int> temp, temp2;

    // TFile *_file0 = TFile::Open("/mnt/c/Work/root/builddir/macros/Calibration/30_08_cosmo_large_2.root");
    TFile *_file0 = TFile::Open("/mnt/c/Work/root/builddir/macros/Calibration/cosmic_27_02_23_large.root");
    TTree *adc64_data = (TTree *)_file0->Get("adc64_data");
    Double_t x[2048] = {0.};
    Double_t y[2048] = {0.};
    double zl = 0;
    vector<string> branch_names;
    for (int cell_iter = 0; cell_iter < numCols * numRows; cell_iter++) // horiz
        if (cells_arr_for_loop[cell_iter] != -1)
            branch_names.push_back("channel_" + std::to_string(cells_arr_for_loop[cell_iter]));
    print(branch_names); // check
    vector<string> seven_cells;
    // vector<string> names_fired;
    string xxx = "";
    TGraph *g1 = new TGraph(); // using the blank constructor
    TCanvas *cfs = new TCanvas();
    TString pdf_wf = "test";

    cfs->Print("/mnt/c/Work/root/builddir/macros/Calibration/pics/" + pdf_wf + ".pdf(", "pdf");

    TH1F *cells_fired = new TH1F("cells_fired", "cells_fired", 40, 0, 40);
    TH1F *Section[128];
    for (int i = 0; i < 128; i++)
        Section[i] = new TH1F(TString::Format("amp%d", i), TString::Format("amp%d", i), 100, 0., 10000.);
    const int MAX_N_SAMPLES = 2048;
    Short_t wf[MAX_N_SAMPLES];
    int how_many_cells = 0;
    for (int j = 0; j < 15; j++)
    {
        xxx = "";
        for (int i = 0; i < 7; i++)
        {
            // xxx += "channel_" + std::to_string(map[i][j]);
            temp.push_back(map[j][i]);
        }
        // sort(temp.begin(), temp.end());
        for (int ii = 0; ii < 7; ii++)
            xxx += std::to_string(temp.at(ii)) + "_";
        seven_cells.push_back(xxx);
        temp.clear();
    }
    // print(seven_cells);

    // for (int cell_iter = 0; cell_iter < 63; cell_iter++)
    // branch_names.push_back("channel_" + std::to_string(cell_iter));
    // print(branch_names); // check
    int n_ev_to_fill = 0;
    int n_events_to_loop = 1500000; // 23763050

    // int map_arr[numRows][numCols];
    // memset(map_arr, 0, numRows * numCols * sizeof(int));

    int **map_arr; // pointer to hold address
    map_arr = vectorToArray(map);

    for (int i = 0; i < numRows; i++)
    {
        for (int j = 0; j < numCols; j++)
            cout << map_arr[i][j] << " ";
        cout << endl;
    }
    int names_fired_arr[numRows][numCols];
    cout << " EVENTS : " << adc64_data->GetEntries() << endl;
    for (int entry = 0; entry < n_events_to_loop; entry++) // adc64_data->GetEntries() // 13586 31477 87428 6336
    {
        bool flag = false;
        memset(names_fired_arr, 0, sizeof names_fired_arr);
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
            // cout << name << endl;
            double ss = channel->GetLeaf("wf_size")->GetValue();
            if (ss != 0)
            {
                // cout << entry << " " << name << " wf size " << ss << endl;
                string s2(name(8, 10));
                // names_fired += name;
                temp2.push_back(stoi(s2));
                how_many_cells++;
                // cout << "Channel number: " << s2 << " INT " << ch_for_traction << endl;
            }
        }

        if (how_many_cells > 6)
        {
            // cout << "Entry: " << entry << endl;
            sort(temp2.begin(), temp2.end());
            // print_int(temp2);
            // cout << " SIZE " << temp2.size() << endl;
            for (int i = 0; i < temp2.size(); i++)
            {
                // std::cout << i << " " << input.at(i) << endl;
                // for (int ii = 0; ii < how_many_cells; ii++)
                names_fired += std::to_string(temp2.at(i)) + "_";

                for (int ii = 0; ii < numRows; ii++)
                    for (int j = 0; j < numCols; j++)
                        if (map_arr[ii][j] == temp2.at(i))
                            names_fired_arr[ii][j] = temp2.at(i);
            }
        }
        temp2.clear();
        vector<string> subsequence;
        subsequence.clear();
        // cout << entry << " --! " << names_fired << " " << how_many_cells << endl;
        if (how_many_cells > 6)
        {
            /*for (int i = 0; i < numRows; i++)
            {
                //cout << i << " | ";
                for (int j = 0; j < numCols; j++)
                    //cout << names_fired_arr[i][j] << " ";
                cout << endl;
            }*/
            // cout << names_fired << endl;
            for (int rowIndex = 0; rowIndex < 15; rowIndex++)
            {
                int lengthA = numCols;
                int lengthB = numCols;
                memset(row_fired, -1, sizeof row_fired);
                memset(row_from_map, -1, sizeof row_from_map);
                for (int i = 0; i < numCols; i++)
                {
                    row_fired[i] = names_fired_arr[rowIndex][i];
                    row_from_map[i] = map_arr[rowIndex][i];
                }

                flag = find_subseq(entry, row_from_map, lengthA, row_fired, lengthB, subsequence);
                if (flag)
                    break;
            }
        }

        if (flag)
        {
            /*cout << "ENTRY " << entry << endl;
            for (int i = 0; i < subsequence.size(); i++)
                cout << subsequence[i] << " ";
            cout << endl;*/

            for (auto i_name = subsequence.begin(); i_name != subsequence.end(); ++i_name)
            {
                int amp = 40000, gate_beg = 80, gate_end = 120;
                zl = 0;
                double zl_min = 0, zl_max = 100000;
                int peak = 0;

                TString name = *i_name;
                TBranch *channel = (TBranch *)adc64_data->GetBranch(name);
                double ss = channel->GetLeaf("wf_size")->GetValue();
                if (ss != 0)
                {

                    for (int i = 0; i < 80; i++) // zl measurement
                    {
                        zl += channel->GetLeaf("wf")->GetValue(i) / 80;
                        if (channel->GetLeaf("wf")->GetValue(i) > zl_min)
                            zl_min = channel->GetLeaf("wf")->GetValue(i);
                        if (channel->GetLeaf("wf")->GetValue(i) < zl_max)
                            zl_max = channel->GetLeaf("wf")->GetValue(i);
                    }
                    for (Int_t i = gate_beg; i <= gate_end; i++)
                    {
                        double xx = channel->GetLeaf("wf")->GetValue(i);
                        if (xx < amp)
                        {
                            amp = xx;
                            peak = i;
                        }
                    }
                    string s3(name(8, 10));
                    int sect_to_write = stoi(s3);
                    int iii = 0;
                    for (int kk = peak - 3; kk <= peak + 7; kk++)
                    {
                        if (kk < peak)
                            if (channel->GetLeaf("wf")->GetValue(kk) - channel->GetLeaf("wf")->GetValue(kk - 1) < 0)
                                iii++;
                        if (kk > peak)
                            if (channel->GetLeaf("wf")->GetValue(kk) - channel->GetLeaf("wf")->GetValue(kk - 1) > 0)
                                iii++;
                        // data[iii] = channel->GetLeaf("wf")->GetValue(kk);
                        // sum10 += channel->GetLeaf("wf")->GetValue(kk) - zl;
                    }

                    if (amp > 0 && amp > 2 * abs(zl_min - zl_max) & iii > 8)
                        Section[sect_to_write]->Fill(-amp + zl);
                    // cout << entry << " " << name << " amp " << amp << " zl " << zl << " zl_min " << zl_min << " zl_max " << zl_max << " a-zl " << -amp + zl << " sect " << sect_to_write << " 2 * abs(zl_min - zl_max) " << 2 * abs(zl_min - zl_max) << endl;
                    n_ev_to_fill++;
                    ///////////////////////////////////////Draw wf //////////////////////////////////
                    if (sect_to_write == 93 & iii > 8) // || sect_to_write == 44 || sect_to_write == 44)
                    {
                        TCanvas *cpeak = new TCanvas();
                        for (Int_t i = 0; i <= gate_end; i++)
                        {
                            g1->SetPoint(g1->GetN(), i, channel->GetLeaf("wf")->GetValue(i));
                        }
                        g1->SetMarkerStyle(8);
                        // g1->GetYaxis()->SetRangeUser(zl - 500, amp + 100);

                        g1->Draw("");
                        int size_d = 0;
                        size_d = peak + 5 - (peak - 5);
                        array<double, 10> data;
                        // double data[size_d];
                        int iii = 0;
                        double sum10 = 0;
                        for (int kk = peak - 3; kk <= peak + 7; kk++)
                        {
                            if (kk < peak)
                                if (channel->GetLeaf("wf")->GetValue(kk) - channel->GetLeaf("wf")->GetValue(kk - 1) < 0)
                                    iii++;
                            if (kk > peak)
                                if (channel->GetLeaf("wf")->GetValue(kk) - channel->GetLeaf("wf")->GetValue(kk - 1) > 0)
                                    iii++;
                            // data[iii] = channel->GetLeaf("wf")->GetValue(kk);
                            // sum10 += channel->GetLeaf("wf")->GetValue(kk) - zl;
                        }
                        double stdev = 0.0, cov = 0.0;
                        // calculateStdDev(data, stdev, cov);
                        g1->SetTitle(TString::Format("ev %i | ", entry) + name + TString::Format(" iii %i | ", iii) + TString::Format(" ampl/sum10 %f | ", amp / sum10) + TString::Format(" stddev %f | ", stdev) + TString::Format(" cov %f | ", cov) + TString::Format(" Peak %i | ", amp) + TString::Format(" ZL %f | ", zl) + TString::Format(" Ampl %f | ", -amp + zl));
                        cpeak->Print("/mnt/c/Work/root/builddir/macros/Calibration/pics/" + pdf_wf + ".pdf", "pdf");
                        g1->Set(0);
                    }
                }
                flag = false;
            }

            /*auto c16 = new TCanvas("c16", "c16", 600, 400);
            // gStyle->SetOptStat(kFALSE);
            auto h3box = new TH3F("h3box", "Option BOX", 8, 0, 8, 3, 0, 3, 3, 0, 3);

            for (int i = 0; i < 9; i++)
            {
                for (int j = 0; j < 7; j++)
                    for (int k = 0; k < 9; k++)
                        if (map[i][j] == channels_for_traction[k] && channels_for_traction[k] != 0) // && line_to_draw == i)
                        {
                            {
                                h3box->Fill(j, mapxyz[i][0], mapxyz[i][1], 1);
                            }
                            // cout << i << " " << j << " " << k << " MAP " << map[i][j] << "  channels_for_traction[k] " << channels_for_traction[k] << " xyz " << j << " xyz " << mapxyz[i][0] << " xyz " << mapxyz[i][1] << endl;
                            //   line[i] += 1;
                        }
                        else
                        {
                            h3box->Fill(j, mapxyz[i][0], mapxyz[i][1], 10);
                        }
                h3box->Draw("BOX2 Z");
            }*/
        }
    }

    ////////Draw///////////
    cout << n_ev_to_fill << endl;
    TCanvas *c1 = new TCanvas();
    c1->DivideSquare(144);
    int i_array = 0;
    TFile histoFileFull("calib_hor.root", "RECREATE");
    cout << " SIZE " << branch_names.size() << endl;
    for (auto i_name = branch_names.begin(); i_name != branch_names.end(); i_name++)
    {
        TString name = *i_name;
        cout << name << endl;
        string s3(name(8, 10));
        int sect_to_write = stoi(s3);

        cout << sect_to_write << " " << Section[sect_to_write]->Integral() << endl;
        c1->cd(sect_to_write);
        double lvl = 0;
        int peak = 0;
        for (int i = 10; i < 80; i++)
            if (Section[sect_to_write]->GetBinContent(i) > lvl)
            {
                lvl = Section[sect_to_write]->GetBinContent(i);
                peak = i * 100;
            }
        TF1 *fit1 = new TF1("fit1", "gaus", peak - 600, peak + 750);
        Section[sect_to_write]->Draw();
        Section[sect_to_write]->Write();
        Section[sect_to_write]->SetTitle(name);
        Section[sect_to_write]->Fit(fit1, "R");
        fit1->Draw("same");
        myfile << sect_to_write << " " << fit1->GetParameter(1) << " " << fit1->GetParameter(2) << endl;
        c1->Update();
        // TCanvas *c2 = new TCanvas();
        /*Section[i_array]->Draw("hist");
        Section[i_array]->SetTitle(name);
        Section[i_array]->GetXaxis()->SetTitle("Ampl [channels]");
        Section[i_array]->GetYaxis()->SetTitle("counts");
        Section[i_array]->Fit(fit1, "QR");
        fit1->Draw("same");*/
        // myfile << i_array << " " << fit1->GetParameter(1) << " " << fit1->GetParameter(2) << endl;
        // c2->Print(s + ".pdf", "pdf");
    }
    c1->Print("c2.pdf", "pdf");
    // c1->Write();
    TCanvas *c2 = new TCanvas();
    cells_fired->Scale(1. / n_events_to_loop);
    cells_fired->Draw("hist");
    cells_fired->GetXaxis()->SetTitle("cells fired");
    cells_fired->GetYaxis()->SetTitle("counts scaled to n_{events}");
    myfile.close();
    cfs->Print("/mnt/c/Work/root/builddir/macros/Calibration/pics/" + pdf_wf + ".pdf)", "pdf");
}