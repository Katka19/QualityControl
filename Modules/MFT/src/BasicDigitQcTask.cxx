// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

///
/// \file   BasicDigitQcTask.cxx
/// \author Tomas Herman
/// \author Guillermo Contreras
/// \author Katarina Krizkova Gajdosova
///

// ROOT
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
// O2
#include <DataFormatsITSMFT/Digit.h>
#include <Framework/InputRecord.h>
// Quality Control
#include "QualityControl/QcInfoLogger.h"
#include "MFT/BasicDigitQcTask.h"
// C++
#include <fstream>

using namespace o2;
using namespace o2::itsmft;

namespace o2::quality_control_modules::mft
{

const int gnZones = 4;
const int gnChips = 9;
int gHalf1Disk0Layer0[gnZones][gnChips] = {
  {485, 484, 483, 482, 481, 480, -1, 469, 468},//zone0
  {494, 493, 492, 491, 490, 489, 488, 487, 486},//zone1
  {503, 502, 501, 500, 499, 498, 497, 496, 495},//zone2
  {-1, 473, 472, -1, 471, 470, 506, 505, 504}};//zone3};
int gHalf0Disk0Layer0[gnZones][gnChips] = {
  {0, 1, -1, 12, 13, 14, 15, 16, 17},//zone0
  {18, 19, 20, 21, 22, 23, 24, 25, 26},//zone1
  {27, 28, 29, 30, 31, 32, 33, 34, 35},//zone2
  {36, 37, 38, 2, 3, -1, 4, 5, -1}};//zone3};
int gHalf1Disk0Layer1[gnZones][gnChips] = {
  {-1, 479, 478, 533, 532, 531, 530, 529, 528},//zone0
  {527, 526, 525, 524, 523, 522, 521, 520, 519},//zone1
  {518, 517, 516, 515, 514, 513, 512, 511, 510},//zone2
  {509, 50, 507, -1, 477, 476, -1, 475, 474}};//zone3
int gHalf0Disk0Layer1[gnZones][gnChips] = {
  {60, 61, 62, 63, 64, 65, 10, 11, -1},//zone0
  {51, 52, 53, 54, 55, 56, 57, 58, 59},//zone1
  {42, 43, 44, 45, 46, 47, 48, 49, 50},//zone2
  {6, 7, -1, 8, 9, -1, 39, 40, 41}};//zone3
int gHistBin[gnChips][2] = {{0,2}, {0,1}, {0,0}, {1,2}, {1,1}, {1,0}, {2,2}, {2,1}, {2,0}};

const double gPixelHitMapsMaxBinX = 1023.5;
const double gPixelHitMapsMaxBinY = 511.5;
const double gPixelHitMapsMinBin = -0.5;
const int gPixelHitMapsBinWidth = 8;

BasicDigitQcTask::~BasicDigitQcTask()
{
  /*
    not needed for unique pointers
  */
}

void BasicDigitQcTask::initialize(o2::framework::InitContext& /*ctx*/)
{
  ILOG(Info, Support) << "initialize BasicDigitQcTask" << ENDM; // QcInfoLogger is used. FairMQ logs will go to there as well.

  // this is how to get access to custom parameters defined in the config file at qc.tasks.<task_name>.taskParameters
  if (auto param = mCustomParameters.find("FLP"); param != mCustomParameters.end()) {
    ILOG(Info, Support) << "Custom parameter - myOwnKey: " << param->second << ENDM;
    FLP = stoi(param->second);
    minChipID = (FLP - 1) * nchip / 4;
    maxChipID = FLP * nchip / 4;
  }

  //  -------------------
  mMFT_chip_index_H = std::make_unique<TH1F>("ChipHitMaps/mMFT_chip_index_H", "mMFT_chip_index_H", 936, -0.5, 935.5);
  getObjectsManager()->startPublishing(mMFT_chip_index_H.get());
  getObjectsManager()->addMetadata(mMFT_chip_index_H->GetName(), "custom", "34");

 
  //  MFT BOTTOM HALF -> HALF_0
  for(int iZone = 0; iZone < gnZones; iZone++)
  {
    auto hitmaphalf0disk0layer0 = std::make_unique<TH2F>(Form("Half_0/Disk_0/Layer_0/mMFTHitMap_Zone%d", iZone), Form("MFT ChipHitMap: Half 0, Disk 0, Layer 0, Zone %d", iZone), 3, 0, 3, 4, 0, 4);
    hitmaphalf0disk0layer0->GetXaxis()->SetNdivisions(300);
    hitmaphalf0disk0layer0->GetYaxis()->SetNdivisions(300);
    hitmaphalf0disk0layer0->SetStats(0);
    mMFTHitMap_Half0Disk0Layer0.push_back(std::move(hitmaphalf0disk0layer0));
    getObjectsManager()->startPublishing(mMFTHitMap_Half0Disk0Layer0[iZone].get());

    auto hitmaphalf0disk0layer1 = std::make_unique<TH2F>(Form("Half_0/Disk_0/Layer_1/mMFTHitMap_Zone%d", iZone), Form("MFTHitMap: Half 0, Disk 0, Layer 1, Zone %d", iZone), 3, 0, 3, 4, 0, 4);
    hitmaphalf0disk0layer1->GetXaxis()->SetNdivisions(300);
    hitmaphalf0disk0layer1->GetYaxis()->SetNdivisions(300);
    hitmaphalf0disk0layer1->SetStats(0);
    mMFTHitMap_Half0Disk0Layer1.push_back(std::move(hitmaphalf0disk0layer1));
    getObjectsManager()->startPublishing(mMFTHitMap_Half0Disk0Layer1[iZone].get());
  }//	end of loop over zones: defining histograms of integrated hits per chip


  //  MFT TOP HALF -> HALF_1
  for(int iZone = 0; iZone < gnZones; iZone++)
  {
    auto hitmaphalf1disk0layer0 = std::make_unique<TH2F>(Form("Half_1/Disk_0/Layer_0/mMFTHitMap_Zone%d", iZone), Form("MFTHitMap: Half 1, Disk 0, Layer 0, Zone %d", iZone), 3, 0, 3, 4, 0, 4);
    hitmaphalf1disk0layer0->GetXaxis()->SetNdivisions(300);
    hitmaphalf1disk0layer0->GetYaxis()->SetNdivisions(300);
    hitmaphalf1disk0layer0->SetStats(0);
    mMFTHitMap_Half1Disk0Layer0.push_back(std::move(hitmaphalf1disk0layer0));
    getObjectsManager()->startPublishing(mMFTHitMap_Half1Disk0Layer0[iZone].get());

    auto hitmaphalf1disk0layer1 = std::make_unique<TH2F>(Form("Half_1/Disk_0/Layer_1/mMFTHitMap_Zone%d", iZone), Form("MFTHitMap: Half 1, Disk 0, Layer 1, Zone %d", iZone), 3, 0, 3, 4, 0, 4);
    hitmaphalf1disk0layer1->GetXaxis()->SetNdivisions(300);
    hitmaphalf1disk0layer1->GetYaxis()->SetNdivisions(300);
    hitmaphalf1disk0layer1->SetStats(0);
    mMFTHitMap_Half1Disk0Layer1.push_back(std::move(hitmaphalf1disk0layer1));
    getObjectsManager()->startPublishing(mMFTHitMap_Half1Disk0Layer1[iZone].get());
  }//	end of loop over zones: defining histograms of integrated hits per chip

 
  //	canvases to plot the 4 zones on it
  mMFT_canvas_HitMap_Half0Disk0Layer0 = std::make_unique<TCanvas>("Half_0/Disk_0/Layer_0/mMFTHitMap", "mMFTHitMap_Half0_Disk0_Layer0", 700, 700);
  mMFT_canvas_HitMap_Half0Disk0Layer0->Divide(4,1,0,0);
  getObjectsManager()->startPublishing(mMFT_canvas_HitMap_Half0Disk0Layer0.get());

  mMFT_canvas_HitMap_Half1Disk0Layer0 = std::make_unique<TCanvas>("Half_1/Disk_0/Layer_0/mMFTHitMap", "mMFTHitMap_Half1_Disk0_Layer0", 700, 700);
  mMFT_canvas_HitMap_Half1Disk0Layer0->Divide(4,1,0,0);
  getObjectsManager()->startPublishing(mMFT_canvas_HitMap_Half1Disk0Layer0.get());

  mMFT_canvas_HitMap_Half0Disk0Layer1 = std::make_unique<TCanvas>("Half_0/Disk_0/Layer_1/mMFTHitMap", "mMFTHitMap_Half0_Disk0_Layer1", 700, 700);
  mMFT_canvas_HitMap_Half0Disk0Layer1->Divide(4,1,0,0);
  getObjectsManager()->startPublishing(mMFT_canvas_HitMap_Half0Disk0Layer1.get());

  mMFT_canvas_HitMap_Half1Disk0Layer1 = std::make_unique<TCanvas>("Half_1/Disk_0/Layer_1/mMFTHitMap", "mMFTHitMap_Half1_Disk0_Layer1", 700, 700);
  mMFT_canvas_HitMap_Half1Disk0Layer1->Divide(4,1,0,0);
  getObjectsManager()->startPublishing(mMFT_canvas_HitMap_Half1Disk0Layer1.get());


  //	------------------------
  //	draw histos on canvas
  for(int iZone = 0; iZone < gnZones; iZone++)
  {
    mMFT_canvas_HitMap_Half1Disk0Layer0->cd(iZone+1);
    mMFTHitMap_Half1Disk0Layer0[gnZones-1-iZone]->Draw("text colz");

    mMFT_canvas_HitMap_Half0Disk0Layer0->cd(iZone+1);
    mMFTHitMap_Half0Disk0Layer0[iZone]->Draw("text colz");

    //	note: rear layer has opposite zone numbering
    mMFT_canvas_HitMap_Half1Disk0Layer1->cd(iZone+1);
    mMFTHitMap_Half1Disk0Layer1[iZone]->Draw("text colz");

    mMFT_canvas_HitMap_Half0Disk0Layer1->cd(iZone+1);
    mMFTHitMap_Half0Disk0Layer1[gnZones-1-iZone]->Draw("text colz");
  }

  mMFT_chip_std_dev_H = std::make_unique<TH1F>("ChipHitMaps/mMFT_chip_std_dev_H", "mMFT_chip_std_dev_H", 936, -0.5, 935.5);
  getObjectsManager()->startPublishing(mMFT_chip_std_dev_H.get());
  getObjectsManager()->addMetadata(mMFT_chip_std_dev_H->GetName(), "custom", "34");

  //==============================================
  //  chip hit maps
  readTable();
  for (int iHitMap = 0; iHitMap < nhitmaps; iHitMap++) {
    //  generate folder and histogram name using the mapping table
    TString FolderName = "";
    TString HistogramName = "";
    getChipName(FolderName, HistogramName, iHitMap);

    auto chiphitmap = std::make_unique<TH2F>(
      FolderName, HistogramName,
      binsChipHitMaps[iHitMap][0], binsChipHitMaps[iHitMap][1], binsChipHitMaps[iHitMap][2],
      binsChipHitMaps[iHitMap][3], binsChipHitMaps[iHitMap][4], binsChipHitMaps[iHitMap][5]);
    chiphitmap->SetStats(0);
    mMFTChipHitMap.push_back(std::move(chiphitmap));
    getObjectsManager()->startPublishing(mMFTChipHitMap[iHitMap].get());
    getObjectsManager()->addMetadata(mMFTChipHitMap[iHitMap]->GetName(), "custom", "34");
  }

  //==============================================
  //  pixel hit maps
  for (int iChipID = 0; iChipID < nchip; iChipID++) {
    //  generate folder and histogram name using the mapping table
    TString FolderName = "";
    TString HistogramName = "";
    getPixelName(FolderName, HistogramName, iChipID);

    //  create pixel hit map
    auto pxlhitmap = std::make_unique<TH2F>(
      FolderName, HistogramName,
      gPixelHitMapsMaxBinX / gPixelHitMapsBinWidth, gPixelHitMapsMinBin, gPixelHitMapsMaxBinX,
      gPixelHitMapsMaxBinY / gPixelHitMapsBinWidth, gPixelHitMapsMinBin, gPixelHitMapsMaxBinY);
    pxlhitmap->SetStats(0);
    mMFTPixelHitMap.push_back(std::move(pxlhitmap));

    if ((iChipID >= minChipID) && (iChipID < maxChipID)) {
      getObjectsManager()->startPublishing(mMFTPixelHitMap[iChipID].get());
      getObjectsManager()->addMetadata(mMFTPixelHitMap[iChipID]->GetName(), "custom", "34");
    }
  }
}

void BasicDigitQcTask::startOfActivity(Activity& /*activity*/)
{
  ILOG(Info, Support) << "startOfActivity" << ENDM;

  mMFT_chip_index_H->Reset();
  mMFT_chip_std_dev_H->Reset();

  for (int iHitMap = 0; iHitMap < nhitmaps; iHitMap++) {
    mMFTChipHitMap[iHitMap]->Reset();
  }

  for (int iChipID = 0; iChipID < nchip; iChipID++) {
    mMFTPixelHitMap[iChipID]->Reset();
  }
}

void BasicDigitQcTask::startOfCycle()
{
  ILOG(Info, Support) << "startOfCycle" << ENDM;
}

void BasicDigitQcTask::monitorData(o2::framework::ProcessingContext& ctx)
{
  // get the digits
  const auto digits = ctx.inputs().get<gsl::span<o2::itsmft::Digit>>("randomdigit");
  if (digits.size() < 1)
    return;

  // fill the histograms
  for (auto& one_digit : digits) {
    int chipIndex = one_digit.getChipIndex();

    //  fill pixel hit maps
    mMFTPixelHitMap[chipIndex]->Fill(one_digit.getColumn(), one_digit.getRow());
    // fill number of entries and standard dev for all chips
    mMFT_chip_index_H->SetBinContent(chipIndex, mMFTPixelHitMap[chipIndex]->GetEntries());
    mMFT_chip_std_dev_H->SetBinContent(chipIndex, mMFTPixelHitMap[chipIndex]->GetStdDev(1));
  }

  //  fill the chip hit maps
  for (int iChipID = 0; iChipID < nchip; iChipID++) {
    int nEntries = mMFTPixelHitMap[iChipID]->GetEntries();
    mMFTChipHitMap[layer[iChipID] + half[iChipID] * 10]->SetBinContent(binx[iChipID], biny[iChipID], nEntries);
  }

  //for(auto& histogram : mMFTD0FrontLowerZ0HitMap)
  //	histogram->Draw("colz");


}

void BasicDigitQcTask::endOfCycle()
{
  ILOG(Info, Support) << "endOfCycle" << ENDM;
}

void BasicDigitQcTask::endOfActivity(Activity& /*activity*/)
{
  ILOG(Info, Support) << "endOfActivity" << ENDM;
}

void BasicDigitQcTask::reset()
{
  // clean all the monitor objects here
  ILOG(Info, Support) << "Resetting the histogram" << ENDM;

  mMFT_chip_index_H->Reset();
  mMFT_chip_std_dev_H->Reset();

  for (int iHitMap = 0; iHitMap < nhitmaps; iHitMap++) {
    mMFTChipHitMap[iHitMap]->Reset();
  }

  for (int iChipID = 0; iChipID < nchip; iChipID++) {
    mMFTPixelHitMap[iChipID]->Reset();
  }
}

void BasicDigitQcTask::getChipName(TString& FolderName, TString& HistogramName, int iHitMap)
{
  FolderName = Form("ChipHitMaps/Half_%d/Disk_%d/Face_%d/mMFTChipHitMap",
                    int(iHitMap / 10), int((iHitMap % 10) / 2), (iHitMap % 10) % 2);

  HistogramName = Form("h%d-d%d-f%d;x (cm);y (cm)",
                       int(iHitMap / 10), int((iHitMap % 10) / 2), (iHitMap % 10) % 2);
}

void BasicDigitQcTask::getPixelName(TString& FolderName, TString& HistogramName, int iChipID)
{
  FolderName = Form("PixelHitMaps/Half_%d/Disk_%d/Face_%d/mMFTPixelHitMap-z%d-l%d-s%d-tr%d",
                    half[iChipID], disk[iChipID], face[iChipID], zone[iChipID], ladder[iChipID], sensor[iChipID], transID[iChipID]);

  HistogramName = Form("h%d-d%d-f%d-z%d-l%d-s%d-tr%d",
                       half[iChipID], disk[iChipID], face[iChipID], zone[iChipID], ladder[iChipID], sensor[iChipID], transID[iChipID]);
}

void BasicDigitQcTask::readTable()
{
  //const int nchip = 936;

  //  reset arrays
  for (int i = 0; i < nchip; i++) {
    half[i] = 0;
    disk[i] = 0;
    face[i] = 0;
    zone[i] = 0;
    ladder[i] = 0;
    sensor[i] = 0;
    transID[i] = 0;
    layer[i] = 0;
    x[i] = 0;
    y[i] = 0;
    z[i] = 0;
    binx[i] = 0;
    biny[i] = 0;
  }

  // read file
  std::ifstream read_table;
  read_table.open("./table_file_binidx.txt");
  for (int i = 0; i < nchip; ++i) {
    read_table >> half[i];
    read_table >> disk[i];
    read_table >> face[i];
    read_table >> zone[i];
    read_table >> ladder[i];
    read_table >> sensor[i];
    read_table >> transID[i];
    read_table >> layer[i];
    read_table >> x[i];
    read_table >> y[i];
    read_table >> z[i];
    read_table >> binx[i];
    read_table >> biny[i];
  }
  read_table.close();
}

} // namespace o2::quality_control_modules::mft
