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
/// \file    runMFTClustersRootFileReader.cxx
/// \author  Guillermo Contreras
/// \author  Tomas Herman
/// \author  Katarina Krizkova Gajdosova
/// \author  Diana Maria Krupova
///
/// \brief This is an executable that reads clusters from a root file from disk and sends the data to QC via DPL.
///
/// This is an executable that reads clusters from a root file from disk and sends the data to QC via the Data Processing Layer.
/// The idea is based in a similar reader by Andrea Ferrero and the DigitReaderSpec definition in MFT
//  It can be used as a data source for QC development. For example, one can do:
/// \code{.sh}
/// o2-qc-mft-clusters-root-file-reader --mft-cluster-infile=some_data_file | o2-qc --config json://${QUALITYCONTROL_ROOT}/etc/your_config.json
/// \endcode
///

// C++
#include <vector>
// ROOT
#include <TFile.h>
#include <TTree.h>
// O2
#include <Framework/CallbackService.h>
#include <Framework/ControlService.h>
#include <Framework/runDataProcessing.h>
#include <Framework/Task.h>
#include <DataFormatsITSMFT/Cluster.h>
#include <DataFormatsITSMFT/CompCluster.h>
#include <DataFormatsITSMFT/ROFRecord.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::itsmft;

class MFTClustersRootFileReader : public o2::framework::Task
{
 public:
  //_________________________________________________________________________________________________
  void init(framework::InitContext& ic)
  {
    LOG(INFO) << " In MFTClustersRootFileReader::init ... entering ";

    // open the input file
    auto filename = ic.options().get<std::string>("mft-cluster-infile");
    mFile = std::make_unique<TFile>(filename.c_str(), "OLD");
    if (!mFile->IsOpen()) {
      LOG(ERROR) << "MFTClustersRootFileReader::init. Cannot open the file: " << filename.c_str();
      ic.services().get<ControlService>().endOfStream();
      ic.services().get<ControlService>().readyToQuit(QuitRequest::Me);
      return;
    }
  }

  //_________________________________________________________________________________________________

  void run(framework::ProcessingContext& pc)
  {
    // get vector of ROF
    std::unique_ptr<TTree> tree((TTree*)mFile->Get("o2sim"));
    std::vector<o2::itsmft::ROFRecord> rofs, *profs = &rofs;
    tree->SetBranchAddress("MFTClustersROF", &profs);
    tree->GetEntry(0);

    // Check if there is a new ROF
    auto numberOfROF = rofs.size();
    if (mCurrentROF >= numberOfROF) {
      // if (mCurrentROF >= 50) {
      LOG(INFO) << " MFTClustersRootFileReader::run. End of file reached";
      pc.services().get<ControlService>().endOfStream();
      pc.services().get<ControlService>().readyToQuit(QuitRequest::Me);
      return;
    }
    // prepare the rof output
    std::vector<o2::itsmft::ROFRecord>* oneROFvec = new std::vector<o2::itsmft::ROFRecord>();
    std::copy(rofs.begin() + mCurrentROF, rofs.begin() + mCurrentROF + 1, std::back_inserter(*oneROFvec));
    // get the current ROF
    auto& rof = rofs[mCurrentROF];
    mCurrentROF++;

    // get the clusters in current ROF
    // --> get cluster branch
    std::vector<o2::itsmft::CompClusterExt> clusters, *pclusters = &clusters;
    tree->SetBranchAddress("MFTClusterComp", &pclusters);
    tree->GetEntry(0);
    // --> find the ranges
    int index = rof.getFirstEntry();        // first cluster position
    int numberOfClustersInROF = rof.getNEntries(); 
    int lastIndex = index + numberOfClustersInROF;
    // --> fill in the corresponding clusters
    std::vector<o2::itsmft::CompClusterExt>* clustersInROF = new std::vector<o2::itsmft::CompClusterExt>();
    std::copy(clusters.begin() + index, clusters.begin() + lastIndex, std::back_inserter(*clustersInROF));

    // fill in the message
    // LOG(INFO) << " MFTClustersRootFileReader::run. In this ROF there are  " << clustersInROF.size() << " clusters";
    pc.outputs().snapshot(Output{ "MFT", "CLUSTERS", 0, Lifetime::Timeframe }, *clustersInROF);
    pc.outputs().snapshot(Output{ "MFT", "MFTClusterROF", 0, Lifetime::Timeframe }, *oneROFvec);
  }

 private:
  std::unique_ptr<TFile> mFile = nullptr;
  unsigned long mCurrentROF = 0;

}; // end class definition

WorkflowSpec defineDataProcessing(const ConfigContext&)
{
  WorkflowSpec specs; // to return the work flow

  // define the outputs
  std::vector<OutputSpec> outputs;
  outputs.emplace_back("MFT", "CLUSTERS", 0, Lifetime::Timeframe);
  outputs.emplace_back("MFT", "MFTClusterROF", 0, Lifetime::Timeframe);

  // The producer to generate some data in the workflow
  DataProcessorSpec producer{
    "clusters-root-file-reader-mft",
    Inputs{},
    outputs,
    AlgorithmSpec{ adaptFromTask<MFTClustersRootFileReader>() },
    Options{ { "mft-cluster-infile", VariantType::String, "mftclusters.root", { "Name of the input file" } } }
  };
  specs.push_back(producer);

  return specs;
}
