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
#include <ITSMFTReconstruction/ChipMappingMFT.h> 
// Quality Control
#include "QualityControl/QcInfoLogger.h"
#include "MFT/BasicDigitQcTask.h"
#include "MFT/BasicDigitQcTaskConversionTable.h"
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
    ILOG(Info, Support) << "Custom parameter - FLP: " << param->second << ENDM;
    FLP = stoi(param->second);
  }
  if (auto param = mCustomParameters.find("TaskLevel"); param != mCustomParameters.end()) {
    ILOG(Info, Support) << "Custom parameter - TaskLevel: " << param->second << ENDM;
    TaskLevel = stoi(param->second);
  }
  
  getChipMapData();

  //  reset arrays of vector and chip IDs
  resetArrays(vectorIDOfChips, vectorHitMapIDOfChips, vectorIDOfChipMaps);

  // Defining histograms
  //==============================================
  mMFT_chip_index_H = std::make_unique<TH1F>("ChipHitMaps/mMFT_chip_index_H", "mMFT_chip_index_H;Chip ID;#Entries", 936, -0.5, 935.5);
  if (TaskLevel == 3 || TaskLevel == 4)
    getObjectsManager()->startPublishing(mMFT_chip_index_H.get());

  mMFT_chip_std_dev_H = std::make_unique<TH1F>("ChipHitMaps/mMFT_chip_std_dev_H", "mMFT_chip_std_dev_H;Chip ID;Chip std dev", 936, -0.5, 935.5);
  if (TaskLevel == 3 || TaskLevel == 4)
    getObjectsManager()->startPublishing(mMFT_chip_std_dev_H.get());

  // --Chip hit maps
  //==============================================
  for (int iVectorHitMapID = 0; iVectorHitMapID < 4; iVectorHitMapID++) {
    // create only hit maps corresponding to the FLP
    int iHitMapID = getHitMapIndex(iVectorHitMapID);

    //  generate folder and histogram name using the mapping table
    TString FolderName = "";
    TString HistogramName = "";
    getChipName(FolderName, HistogramName, iHitMapID);

    auto chiphitmap = std::make_unique<TH2F>(
      FolderName, HistogramName,
      binsChipHitMaps[iHitMapID][0], binsChipHitMaps[iHitMapID][1], binsChipHitMaps[iHitMapID][2],
      binsChipHitMaps[iHitMapID][3], binsChipHitMaps[iHitMapID][4], binsChipHitMaps[iHitMapID][5]);
    // chiphitmap->SetStats(0);
    mMFTChipHitMap.push_back(std::move(chiphitmap));
    if (TaskLevel == 0 || TaskLevel == 1 || TaskLevel == 4)
      getObjectsManager()->startPublishing(mMFTChipHitMap[iVectorHitMapID].get());
  }

  // --Pixel hit maps
  //==============================================
  for (int iVectorID = 0; iVectorID < (nMaps[FLP] + nMaps[4 - FLP]); iVectorID++) {
    // create only hit maps corresponding to the FLP
    int iChipID = getChipIndex(iVectorID);

    //  generate folder and histogram name using the mapping table
    TString FolderName = "";
    TString HistogramName = "";
    getPixelName(FolderName, HistogramName, iChipID);

    auto pxlhitmap = std::make_unique<TH2F>(
      FolderName, HistogramName,
      gPixelHitMapsMaxBinX / gPixelHitMapsBinWidth, gPixelHitMapsMinBin - gPixelHitMapsShift, gPixelHitMapsMaxBinX - gPixelHitMapsShift,
      gPixelHitMapsMaxBinY / gPixelHitMapsBinWidth, gPixelHitMapsMinBin - gPixelHitMapsShift, gPixelHitMapsMaxBinY - gPixelHitMapsShift);
    pxlhitmap->SetStats(0);
    mMFTPixelHitMap.push_back(std::move(pxlhitmap));
    if (TaskLevel == 2 || TaskLevel == 4)
      getObjectsManager()->startPublishing(mMFTPixelHitMap[iVectorID].get());
  }
}

void BasicDigitQcTask::startOfActivity(Activity& /*activity*/)
{
  ILOG(Info, Support) << "startOfActivity" << ENDM;

  mMFT_chip_index_H->Reset();
  mMFT_chip_std_dev_H->Reset();

  for (int iVectorHitMapID = 0; iVectorHitMapID < 4; iVectorHitMapID++) {
    mMFTChipHitMap[iVectorHitMapID]->Reset();
  }

  for (int iVectorID = 0; iVectorID < (nMaps[FLP] + nMaps[4 - FLP]); iVectorID++) {
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

  // fill the pixel hit maps and overview histograms
  for (auto& one_digit : digits) {

    int chipIndex = one_digit.getChipIndex();

    int vectorIndex = getVectorIndex(chipIndex);
    if(vectorIndex < 0) // if the chip is not from wanted FLP, the array will give -1
      continue;

    // fill pixel hit maps
    mMFTPixelHitMap[vectorIndex]->Fill(one_digit.getColumn(), one_digit.getRow());

    // fill overview histograms
    mMFT_chip_index_H->SetBinContent(chipIndex + 1, mMFTPixelHitMap[vectorIndex]->GetEntries());
    mMFT_chip_std_dev_H->SetBinContent(chipIndex + 1, mMFTPixelHitMap[vectorIndex]->GetStdDev(1));

    // fill integrated chip hit maps 
    int VectorHitMapID = getVectorHitMapIndex(chipIndex);
    if(VectorHitMapID < 0)
      continue;
    mMFTChipHitMap[VectorHitMapID]->Fill(x[chipIndex], y[chipIndex]);

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

  for (int iVectorHitMapID = 0; iVectorHitMapID < 4; iVectorHitMapID++) {
    mMFTChipHitMap[iVectorHitMapID]->Reset();
  }

  for (int iVectorID = 0; iVectorID < (nMaps[FLP] + nMaps[4 - FLP]); iVectorID++) {
    mMFTPixelHitMap[iVectorID]->Reset();
  }
}

void BasicDigitQcTask::getChipName(TString& FolderName, TString& HistogramName, int iHitMapID)
{
  FolderName = Form("ChipHitMaps/Half_%d/Disk_%d/Face_%d/mMFTChipHitMap",
                    int(iHitMapID / 10), int((iHitMapID % 10) / 2), (iHitMapID % 10) % 2);

  HistogramName = Form("h%d-d%d-f%d;x (cm);y (cm)",
                       int(iHitMapID / 10), int((iHitMapID % 10) / 2), (iHitMapID % 10) % 2);
}

void BasicDigitQcTask::getPixelName(TString& FolderName, TString& HistogramName, int iChipID)
{
  FolderName = Form("PixelHitMaps/Half_%d/Disk_%d/Face_%d/mMFTPixelHitMap-z%d-l%d-s%d-tr%d",
                    half[iChipID], disk[iChipID], face[iChipID], zone[iChipID], ladder[iChipID], sensor[iChipID], transID[iChipID]);

  HistogramName = Form("h%d-d%d-f%d-z%d-l%d-s%d-tr%d",
                       half[iChipID], disk[iChipID], face[iChipID], zone[iChipID], ladder[iChipID], sensor[iChipID], transID[iChipID]);
}

void BasicDigitQcTask::getChipMapData()
{

  const o2::itsmft::ChipMappingMFT map;
  auto chipMapData = map.getChipMappingData();

  //  reset arrays
  for (int i = 0; i < nChip; i++) {
    half[i] = 0;
    disk[i] = 0;
    face[i] = 0;
    zone[i] = 0;
    sensor[i] = 0;
    transID[i] = 0;
    layer[i] = 0;
  }

  for(int i = 0; i < nChip; i++)  {
    half[i] = chipMapData[i].half; 
    disk[i] = chipMapData[i].disk; 
    layer[i] = chipMapData[i].layer; 
    face[i] = layer[i]%2; 
    zone[i] = chipMapData[i].zone; 
    //ladder[i] = chipMapData[i].ladder; // not in Bogdan's mapping, I will leave it in our table header file
    sensor[i] = chipMapData[i].localChipSWID; 
    transID[i] = chipMapData[i].cable; 
  }
  
}

int BasicDigitQcTask::getVectorHitMapIndex(int chipID)
{

  int HitMapID = vectorHitMapIDOfChips[chipID];

  int VectorHitMapID = vectorIDOfChipMaps[HitMapID];

  return VectorHitMapID;
}

int BasicDigitQcTask::getHitMapIndex(int VectorHitMapID)
{
  int VectorHitMapHalf = int(VectorHitMapID / 2);

  int HitMapID;
  if (VectorHitMapHalf == 0)
    HitMapID = VectorHitMapID + FLP * 2;
  else
    HitMapID = (VectorHitMapID % 2) + (4 - FLP) * 2 + nHitMaps / 2;

  //  fill the array of vector ID for corresponding hit map
  //  (opposite matching)
  vectorIDOfChipMaps[HitMapID] = VectorHitMapID;

  return HitMapID;
}

int BasicDigitQcTask::getVectorIndex(int chipID)
{

  int vectorID = vectorIDOfChips[chipID];
  
  return vectorID;
}

int BasicDigitQcTask::getChipIndex(int vectorID)
{
  int VectorHalf = 0;
  if (int(vectorID / nMaps[FLP]) < 1)
    VectorHalf = 0;
  else
    VectorHalf = 1;

  int chipID = vectorID + VectorHalf * (-nMaps[FLP] + nChip / 2);

  int maxDisk = 0;
  if (VectorHalf == 0)
    maxDisk = FLP;
  else
    maxDisk = 4 - FLP;

  for (int idisk = 0; idisk < maxDisk; idisk++)
    chipID = chipID + nMaps[idisk];

  //  fill the array of vector index for corresponding chipID 
  //  (opposite matching)
  vectorIDOfChips[chipID] = vectorID; 
  //  fill the array of hit map ID for corresponding chipID
  vectorHitMapIDOfChips[chipID] = layer[chipID] + half[chipID] * nHitMaps / 2;

  return chipID;
}

void BasicDigitQcTask::resetArrays(int array[], int array2[], int array3[])
{

  for(int iChip = 0; iChip < nChip; iChip++) {
    array[iChip] = -1;
    array2[iChip] = -1;
  }

  for(int iMap = 0; iMap < nHitMaps; iMap++)
    array3[iMap] = -1;

}

} // namespace o2::quality_control_modules::mft
