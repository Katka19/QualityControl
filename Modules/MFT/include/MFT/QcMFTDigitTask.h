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
/// \file   QcMFTDigitTask.h
/// \author Tomas Herman
/// \author Guillermo Contreras
/// \author Katarina Krizkova Gajdosova
///

#ifndef QC_MODULE_MFT_DIGIT_TASK_H
#define QC_MODULE_MFT_DIGIT_TASK_H

// ROOT
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
// Quality Control
#include "QualityControl/TaskInterface.h"

using namespace o2::quality_control::core;

namespace o2::quality_control_modules::mft
{

/// \brief MFT Digit QC task
///
/// \author Katarina Krizkova Gajdosova
/// \author Tomas Herman
/// \author Guillermo Contreras
class QcMFTDigitTask final : public TaskInterface
{
 public:
  /// \brief Constructor
  QcMFTDigitTask() = default;
  /// Destructor
  ~QcMFTDigitTask() override;

  // Definition of the methods for the template method pattern
  void initialize(o2::framework::InitContext& ctx) override;
  void startOfActivity(Activity& activity) override;
  void startOfCycle() override;
  void monitorData(o2::framework::ProcessingContext& ctx) override;
  void endOfCycle() override;
  void endOfActivity(Activity& activity) override;
  void reset() override;

 private:
  //  variables
  const double maxBinXPixelHitMap = 1024;
  const double maxBinYPixelHitMap = 512;
  const double minBinPixelHitMap = 0;
  const double shiftPixelHitMap = 0.5;
  const int binWidthPixelHitMap = 1;

  const int numberOfHitMaps = 20;
  const int numberOfChips = 936;

  int mCurrentFLP;
  int mTaskLevel;
  int mNumberOfMapsPerFLP[5] = { 66, 66, 82, 118, 136 };

  //  bin numbers for chip hit maps
  double binsChipHitMaps[20][6] = {
    // half0
    { 12, -10, 10, 4, -12, 0 }, // disk0, face 0
    { 12, -10, 10, 4, -12, 0 }, // disk0, face 1

    { 12, -10, 10, 4, -12, 0 },
    { 12, -10, 10, 4, -12, 0 },

    { 13, -11, 10, 4, -12, 0 },
    { 13, -10, 11, 4, -12, 0 },

    { 16, -13, 14, 5, -15, 0 },
    { 16, -14, 13, 5, -15, 0 },

    { 17, -14, 14, 5, -15, 0 },
    { 17, -14, 14, 5, -15, 0 },

    // half1
    { 12, -10, 10, 4, 0, 12 },
    { 12, -10, 10, 4, 0, 12 },

    { 12, -10, 10, 4, 0, 12 },
    { 12, -10, 10, 4, 0, 12 },

    { 13, -10, 11, 4, 0, 12 },
    { 13, -11, 10, 4, 0, 12 },

    { 16, -14, 13, 5, 0, 15 },
    { 16, -13, 14, 5, 0, 15 },

    { 17, -14, 14, 5, 0, 15 },
    { 17, -14, 14, 5, 0, 15 },
  };

  std::unique_ptr<TH1F> mChipOccupancy = nullptr;
  std::unique_ptr<TH1F> mChipOccupancyStdDev = nullptr;

  std::vector<std::unique_ptr<TH2F>> mChipOccupancyMap;
  std::vector<std::unique_ptr<TH2F>> mPixelOccupancyMap;

  //  functions
  // void readTable();
  int getVectorChipOccupancyMapIndex(int ChipOccupancyMapIndex);
  int getChipOccupancyMapIndex(int VectorChipOccupancyMapIndex);
  int getVectorIndex(int chipIndex);
  int getChipIndex(int vectorIndex);
  void getNameOfChipOccupancyMap(TString& FolderName, TString& HistogramName, int iChipOccupancyMapIndex);
  void getNameOfPixelOccupancyMap(TString& FolderName, TString& HistogramName, int iChipIndex);
};

} // namespace o2::quality_control_modules::mft

#endif // QC_MODULE_MFT_DIGIT_TASK_H
