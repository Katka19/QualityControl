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

// temporary header file with mapping table
#include "MFT/BasicDigitQcTaskConversionTable.h"

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

namespace o2::quality_control_modules::mft
{

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
    //minChipID = (FLP - 1) * nchip / 4;
    //maxChipID = FLP * nchip / 4;
  }

  //readTable();

  //  -------------------
  mMFT_chip_index_H = std::make_unique<TH1F>("ChipHitMaps/mMFT_chip_index_H", "mMFT_chip_index_H", 936, -0.5, 935.5);
  getObjectsManager()->startPublishing(mMFT_chip_index_H.get());
  getObjectsManager()->addMetadata(mMFT_chip_index_H->GetName(), "custom", "34");

  mMFT_chip_std_dev_H = std::make_unique<TH1F>("ChipHitMaps/mMFT_chip_std_dev_H", "mMFT_chip_std_dev_H", 936, -0.5, 935.5);
  getObjectsManager()->startPublishing(mMFT_chip_std_dev_H.get());
  getObjectsManager()->addMetadata(mMFT_chip_std_dev_H->GetName(), "custom", "34");

  //==============================================
  //  chip hit maps
  for (int iHitMap = 0; iHitMap < nhitmaps; iHitMap++) {
    //  generate folder and histogram name using the mapping table
    TString FolderName = "";
    TString HistogramName = "";
    getChipName(FolderName, HistogramName, iHitMap);

    auto chiphitmap = std::make_unique<TH2F>(
      FolderName, HistogramName,
      binsChipHitMaps[iHitMap][0], binsChipHitMaps[iHitMap][1], binsChipHitMaps[iHitMap][2],
      binsChipHitMaps[iHitMap][3], binsChipHitMaps[iHitMap][4], binsChipHitMaps[iHitMap][5]);
    //chiphitmap->SetStats(0);
    mMFTChipHitMap.push_back(std::move(chiphitmap));
    getObjectsManager()->startPublishing(mMFTChipHitMap[iHitMap].get());
    getObjectsManager()->addMetadata(mMFTChipHitMap[iHitMap]->GetName(), "custom", "34");
  }

  //==============================================
  //  pixel hit maps
  for (int iVectorID = 0; iVectorID < (nMaps[FLP] + nMaps[4-FLP]); iVectorID++) {
    //  generate folder and histogram name using the mapping table
    TString FolderName = "";
    TString HistogramName = "";
    getPixelName(FolderName, HistogramName, iVectorID);

    //  create pixel hit map
    auto pxlhitmap = std::make_unique<TH2F>(
      FolderName, HistogramName,
      gPixelHitMapsMaxBinX / gPixelHitMapsBinWidth, gPixelHitMapsMinBin, gPixelHitMapsMaxBinX,
      gPixelHitMapsMaxBinY / gPixelHitMapsBinWidth, gPixelHitMapsMinBin, gPixelHitMapsMaxBinY);
    pxlhitmap->SetStats(0);
    mMFTPixelHitMap.push_back(std::move(pxlhitmap));
    getObjectsManager()->startPublishing(mMFTPixelHitMap[iVectorID].get());
    getObjectsManager()->addMetadata(mMFTPixelHitMap[iVectorID]->GetName(), "custom", "34");
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

  for (int iVectorID = 0; iVectorID < (nMaps[FLP] + nMaps[4-FLP]); iVectorID++) {
    mMFTPixelHitMap[iVectorID]->Reset();
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

    //  simulate QC task on specific FLP (i.e. only digits from disk X will arrive); this can be removed once the code is on FLP
    //if(disk[chipIndex] != FLP) continue;
    if((disk[chipIndex] == 1 && half[chipIndex] == 0) || (disk[chipIndex] == 3 && half[chipIndex] == 1))// corresponds to FLP 1 (FLP=(0,1,2,3,4))
    {

      int vectorIndex = getVectorIndex(chipIndex);
    //  fill pixel hit maps
      mMFTPixelHitMap[vectorIndex]->Fill(one_digit.getColumn(), one_digit.getRow());
    // fill number of entries and standard dev for all chips
      mMFT_chip_index_H->SetBinContent(chipIndex, mMFTPixelHitMap[vectorIndex]->GetEntries());
      mMFT_chip_std_dev_H->SetBinContent(chipIndex, mMFTPixelHitMap[vectorIndex]->GetStdDev(1));
    }
  }

  //  fill the chip hit maps
  for (int iVectorID = 0; iVectorID < (nMaps[FLP] + nMaps[4-FLP]); iVectorID++) {
    int nEntries = mMFTPixelHitMap[iVectorID]->GetEntries();
    int chipID = getChipIndex(iVectorID);
    mMFTChipHitMap[layer[chipID] + half[chipID] * 10]->SetBinContent(binx[chipID], biny[chipID], nEntries);
  }
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

  for (int iVectorID = 0; iVectorID < (nMaps[FLP] + nMaps[4-FLP]); iVectorID++) {
    mMFTPixelHitMap[iVectorID]->Reset();
  }
}

void BasicDigitQcTask::getChipName(TString& FolderName, TString& HistogramName, int iHitMap)
{
  FolderName = Form("ChipHitMaps/Half_%d/Disk_%d/Face_%d/mMFTChipHitMap",
                    int(iHitMap / 10), int((iHitMap % 10) / 2), (iHitMap % 10) % 2);

  HistogramName = Form("h%d-d%d-f%d;x (cm);y (cm)",
                       int(iHitMap / 10), int((iHitMap % 10) / 2), (iHitMap % 10) % 2);
}

void BasicDigitQcTask::getPixelName(TString& FolderName, TString& HistogramName, int iVectorID)
{
  int iChipID = getChipIndex(iVectorID);

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

int BasicDigitQcTask::getVectorIndex(int chipID)
{
/*
  int vectorID = chipID + half[chipID]*(-nchip/2 + nMaps[disk[chipID]]/2);
  
  for(int idisk=0; idisk<disk[chipID]; idisk++)
    vectorID = vectorID - nMaps[idisk]/2;
  
  return vectorID;
*/

  int vectorID = chipID + half[chipID]*(-nchip/2 + nMaps[4-disk[chipID]]);
  
  for(int idisk=0; idisk < disk[chipID]; idisk++)
    vectorID = vectorID - nMaps[idisk];
  
  return vectorID;

}

int BasicDigitQcTask::getChipIndex(int vectorID)
{
/*
  int half = int(vectorID/(nMaps[FLP]/2));
  int chipID = vectorID + half*(-nMaps[FLP]/2 + nchip/2);

  for(int idisk=0; idisk<FLP; idisk++)
    chipID = chipID + nMaps[idisk]/2;

  return chipID;
*/

  int half = 0;
  if(int(vectorID/nMaps[FLP]) < 1)
    half = 0;
  else half = 1;

  int chipID = vectorID + half*(-nMaps[FLP] + nchip/2);

  int maxDisk = 0;
  if(half == 0) maxDisk = FLP;
  else maxDisk = 4-FLP;

  for(int idisk=0; idisk < maxDisk; idisk++)
    chipID = chipID + nMaps[idisk];

  return chipID;

}

} // namespace o2::quality_control_modules::mft
