#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDProducer.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
 
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/InputSourceMacros.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Sources/interface/ProducerSourceBase.h"

#include "FWCore/Concurrency/interface/FunctorTask.h"

#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimDataFormats/GeneratorProducts/interface/LesHouches.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEXMLStringProduct.h"

#include "GeneratorInterface/LHEInterface/interface/LHERunInfo.h"
#include "GeneratorInterface/LHEInterface/interface/LHEEvent.h"
#include "GeneratorInterface/LHEInterface/interface/LHEReader.h"
 
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
 
#include "FWCore/MessageLogger/interface/MessageLogger.h"

/*
typical HNL event
//header: NUP (number of particles in event), IDPRUP (process ID), XWGTUP (xsec weight), 
          SCALUP (scale in GeV), AQEDUP (alphaEW), AQCDUP (alphaS) 
//particle: IDUP (pdgid), ISTUP (status), MOTHUP[0] (1st mother index), MOTHUP[1] (2nd mother index), 
            ICOLUP[0] (color flow id of particle), ICOLUP[1] (color flow id of antiparticle), 
            PUP[0..4] (px,py,pz,E,mass (not on-shell mass)), 
            VTIMUP (invariant proper lifetime in mm), SPINUP (helicity)
            
 8      2 +6.3200181e-04 8.99545700e+01 7.81616400e-03 1.18245400e-01
       -1 -1    0    0    0  501 -0.0000000000e+00 +0.0000000000e+00 +2.5426129881e+02 2.5426129886e+02 5.0400000000e-03 0.0000e+00 1.0000e+00
        1 -1    0    0  501    0 +0.0000000000e+00 -0.0000000000e+00 -7.9562093420e+00 7.9562109384e+00 5.0400000000e-03 0.0000e+00 -1.0000e+00
       23  2    1    2    0    0 +0.0000000000e+00 -1.7763568394e-15 +2.4630508947e+02 2.6221750980e+02 8.9954573798e+01 0.0000e+00 0.0000e+00
  9900012  2    3    3    0    0 +7.1325341723e+00 -1.0096813457e+01 +3.0110936723e-01 1.5903123472e+01 9.9999991814e+00 3.0871e-01 0.0000e+00
      -16  1    3    3    0    0 -7.1325341723e+00 +1.0096813457e+01 +2.4600398010e+02 2.4631438633e+02 0.0000000000e+00 0.0000e+00 1.0000e+00
      -11  1    4    4    0    0 +7.3382156799e-01 -4.7149386384e+00 +8.3504240076e-01 4.8442168128e+00 5.1100000000e-04 0.0000e+00 1.0000e+00
        3  1    4    4  502    0 +1.3117896098e+00 -5.2509093443e+00 -1.7235399759e+00 5.6809886438e+00 1.0100000000e-01 0.0000e+00 -1.0000e+00
       -4  1    4    4    0  502 +5.0869229945e+00 -1.3096547415e-01 +1.1896069424e+00 5.3779180158e+00 1.2700000000e+00 0.0000e+00 1.0000e+00
*/
class HNLGun: 
    public edm::ProducerSourceBase
{
    private:
        const double cme = 13000;
    
        
        
    public:
        explicit HNLGun(const edm::ParameterSet& iConfig,  const edm::InputSourceDescription& desc):
            edm::ProducerSourceBase(iConfig, desc, false)
        {
            //setNewRun();
            //setNewLumi();
 
            produces<LHERunInfoProduct, edm::Transition::BeginRun>();
            produces<LHEEventProduct>();
        }
        
        ~HNLGun()
        {
        }
        
        bool setRunAndEventInfo(edm::EventID&, edm::TimeValue_t&, edm::EventAuxiliary::ExperimentType&) override
        {
            return true;
        }
        
        void beginRun(edm::Run &run) override
        {
            lhef::HEPRUP heprup;
            heprup.resize(1);
            heprup.IDBMUP = {2212,2212}; //beam particle pdgids in +z/-z
            heprup.EBMUP = {cme/2.,cme/2.}; //energy of beam particles
            heprup.PDFGUP = {1,1}; //pdf author group?
            heprup.PDFSUP = {1001,1001}; //pdf id
            heprup.IDWTUP = 4; //treatment of event weights (4 = weighted in/ weighted out)
            
            heprup.XSECUP[0] = 1.; // xsec of process in pb
            heprup.XERRUP[0] = 1e-3; //error of xsec in pb
            heprup.XMAXUP[0] = 1.; //maximum weight of process
            heprup.LPRUP[0] = 1; //process ID
        
            std::unique_ptr<LHERunInfoProduct> runInfo(new LHERunInfoProduct(heprup));
            run.put(std::move(runInfo));
        }
        
        
        void produce(edm::Event& event) override
        {
            
            
            double x1 = 0.1;
            double x2 = 0.1;
            double scale = std::sqrt(x1*x2)*cme;
            
            
            lhef::HEPEUP hepeup;
            hepeup.resize(4); //number of particles
            hepeup.AQCDUP = 0.12; //alphaS
            hepeup.AQEDUP = 1./137.; //alphaEW
            hepeup.IDPRUP = 1; //process ID
            hepeup.SCALUP = scale; //scale in GeV of process
            hepeup.XWGTUP = 1.; //xsec weight
            
            hepeup.IDUP[0] = 1; // particle ID
            hepeup.ISTUP[0] = -1; //status code (-1 incomming, +1 outgoing, +2 intermediate resonance)
            hepeup.MOTHUP[0] = {0,0}; // mother indices (0 = no mother)
            hepeup.ICOLUP[0] = {1,0}; //color tag to indicate color flow ({particle, antiparticle})
            hepeup.PUP[0][0]= 0.; //px
            hepeup.PUP[0][1]= 0.; //py
            hepeup.PUP[0][2]= x1*cme/2.; //pz
            hepeup.PUP[0][3]= x1*cme/2.; //E
            hepeup.PUP[0][4]= 0.; //mass
            hepeup.VTIMUP[0] = 0.; //ctau in mm
            hepeup.SPINUP[0] = 0.; //helicity
            
            hepeup.IDUP[1] = -1; // particle ID
            hepeup.ISTUP[1] = -1; //status code (-1 incomming, +1 outgoing, +2 intermediate resonance)
            hepeup.MOTHUP[1] = {0,0}; // mother indices (0 = no mother)
            hepeup.ICOLUP[1] = {0,1}; //color tag to indicate color flow ({particle, antiparticle})
            hepeup.PUP[1][0]= 0.; //px
            hepeup.PUP[1][1]= 0.; //py
            hepeup.PUP[1][2]= -x2*cme/2.; //pz
            hepeup.PUP[1][3]= x2*cme/2.; //E
            hepeup.PUP[1][4]= 0.; //mass
            hepeup.VTIMUP[1] = 0.; //ctau in mm
            hepeup.SPINUP[1] = 0.; //helicity

            
            //per particle info
            hepeup.IDUP[2] = 1; // particle ID
            hepeup.ISTUP[2] = 1; //status code (-1 incomming, +1 outgoing, +2 intermediate resonance)
            hepeup.MOTHUP[2] = {1,2}; // mother indices (0 = no mother)
            hepeup.ICOLUP[2] = {2,0}; //color tag to indicate color flow ({particle, antiparticle})
            hepeup.PUP[2][0]= -scale/2.; //px
            hepeup.PUP[2][1]= 0.; //py
            hepeup.PUP[2][2]= hepeup.PUP[0][2]+hepeup.PUP[1][2]; //pz
            hepeup.PUP[2][3]= scale/2.; //E
            hepeup.PUP[2][4]= 0.; //mass
            hepeup.VTIMUP[2] = 1.; //ctau in mm
            hepeup.SPINUP[2] = 0.; //helicity
            
            hepeup.IDUP[3] = -1; // particle ID
            hepeup.ISTUP[3] = 1; //status code (-1 incomming, +1 outgoing, +2 intermediate resonance)
            hepeup.MOTHUP[3] = {1,2}; // mother indices (0 = no mother)
            hepeup.ICOLUP[3] = {0,2}; //color tag to indicate color flow ({particle, antiparticle})
            hepeup.PUP[3][0]= scale/2.; //px
            hepeup.PUP[3][1]= 0.; //py
            hepeup.PUP[3][2]= hepeup.PUP[0][2]+hepeup.PUP[1][2]; //pz
            hepeup.PUP[3][3]= scale/2.; //E
            hepeup.PUP[3][4]= 0.; //mass
            hepeup.VTIMUP[3] = 1.; //ctau in mm
            hepeup.SPINUP[3] = 0.; //helicity

            std::unique_ptr<LHEEventProduct> product(new LHEEventProduct(hepeup,1.0));
            event.put(std::move(product));
        }
        
        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions)
        {
            edm::ParameterSetDescription desc;
            descriptions.addDefault(desc);
        }
};

DEFINE_FWK_INPUT_SOURCE(HNLGun);

