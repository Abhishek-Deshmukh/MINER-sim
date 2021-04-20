#include "G4MTRunManager.hh"
#include "G4RunManager.hh"

#include "G4UImanager.hh"

//#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
//#endif

#include "globals.hh"
#include "G4UIExecutive.hh"
#include "G4VPhysicsConstructor.hh"
#include "GeometryConstruction.hh"
#include "Shielding.hh"
#include "ActionInitialization.hh"
#include "RootIO.hh"
#include "HistManager.hh"
// stuff for importance sampling
#include "G4GeometryManager.hh"
#include "ImportanceGeometryConstruction.hh"
#include "MonitorGeometryConstruction.hh"
#include "G4ImportanceBiasing.hh"
#include "G4ParallelWorldPhysics.hh"
#include "G4GeometrySampler.hh"
#include "G4IStore.hh"
#include "G4VWeightWindowStore.hh"
#include "G4WeightWindowAlgorithm.hh"
#include "G4ScoringManager.hh"

int main(int argc,char** argv) {

    G4Random::setTheEngine(new CLHEP::RanecuEngine);
    G4long seeds[2];
    time_t systime = time(NULL);
    seeds[0] = (long) systime;
    seeds[1] = (long) (systime*G4UniformRand());
    G4Random::setTheSeeds(seeds);
    auto * runManager = new G4RunManager;
    G4ScoringManager* scoringManager = G4ScoringManager::GetScoringManager();

    // set mandatory initialization classes
    //MINERMaterials::Instance();
    auto* detector = new GeometryConstruction;
    runManager->SetUserInitialization(detector);

    //detector->RegisterParallelWorld(pdet);

    //MonitorGeometryConstruction *mondet = new MonitorGeometryConstruction("ParallelMonitoringWorld");
    //detector->RegisterParallelWorld(mondet);

    //G4GeometrySampler pgs(pdet->GetWorldVolume(),"gamma");
    //pgs.SetParallel(true);

   // G4GeometrySampler pgs2(pdet->GetWorldVolume(),"neutron");
    //pgs2.SetParallel(true);

    G4VModularPhysicsList *physicsList = new Shielding;
    //physicsList->RegisterPhysics(new G4ParallelWorldPhysics("ParallelBiasingWorld"));
    //physicsList->RegisterPhysics(new G4ImportanceBiasing(&pgs,"ParallelBiasingWorld"));
    //physicsList->RegisterPhysics(new G4ImportanceBiasing(&pgs2,"ParallelBiasingWorld"));
    //physicsList->RegisterPhysics(new G4ParallelWorldPhysics("ParallelMonitoringWorld"));


    /*
     // add thermal neutron model
     G4HadronElasticProcess* theNeutronElasticProcess = new G4HadronElasticProcess;
     // Cross Section Data set
     G4NeutronHPElasticData* theHPElasticData = new G4NeutronHPElasticData;
     theNeutronElasticProcess->AddDataSet(theHPElasticData);
     G4NeutronHPThermalScatteringData* theHPThermalScatteringData = new G4NeutronHPThermalScatteringData;
     theNeutronElasticProcess->AddDataSet(theHPThermalScatteringData);
     // Models
     G4NeutronHPElastic* theNeutronElasticModel = new G4NeutronHPElastic;
     theNeutronElasticModel->SetMinEnergy(4.0*eV);
     theNeutronElasticProcess->RegisterMe(theNeutronElasticModel);
     G4NeutronHPThermalScattering* theNeutronThermalElasticModel = new G4NeutronHPThermalScattering;
     theNeutronThermalElasticModel->SetMaxEnergy(4.0*eV);
     theNeutronElasticProcess->RegisterMe(theNeutronThermalElasticModel);
     // Apply Processes to Process Manager of Neutron
     G4ProcessManager* pmanager = G4Neutron::Neutron()->GetProcessManager();
     pmanager->AddDiscreteProcess(theNeutronElasticProcess);
     */

    runManager->SetUserInitialization(physicsList);
    runManager->SetUserInitialization(new ActionInitialization);
    RootIO::GetInstance();

    runManager->Initialize();

    //pdet->CreateImportanceStore();

    //pgs2.PrepareImportanceSampling(G4IStore::GetInstance(pdet->GetName()),0);
    //pgs2.Configure();



    // visualization manager
    G4VisManager* visManager = new G4VisExecutive;
    visManager->Initialize();


    // get the pointer to the User Interface manager
    G4UImanager* UImanager = G4UImanager::GetUIpointer();

        // Define (G)UI terminal for interactive mode
        if ( argc==1 )   // Automatically run default macro for writing...
        {
            auto* ui = new G4UIExecutive(argc, argv);
            UImanager->ApplyCommand("/control/execute vis.mac");
            ui->SessionStart();
            delete ui;
        }
        else             // Interactive, provides macro in input
        {
            auto* ui = new G4UIExecutive(argc, argv);
            G4String command = "/control/execute ";
            G4String fileName = argv[1];
            UImanager->ApplyCommand(command+fileName);
            ui->SessionStart();
            delete ui;
        }


    // job termination
    //
    G4GeometryManager::GetInstance()->OpenGeometry();
   // pgs.ClearSampling();
   // pgs2.ClearSampling();
    delete visManager;
    delete runManager;

    return 0;
}
