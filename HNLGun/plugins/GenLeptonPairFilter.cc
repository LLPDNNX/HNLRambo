// system include files
#include <memory>
#include <vector>
#include <unordered_map>
// user include files

#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/angle.h"

#include <ctime>
#include <unordered_set>

#include "TH3F.h"
#include "TH1D.h"

#include "TRandom3.h"



class GenLeptonPairFilter:
    public edm::one::EDFilter<>

{
    public:
        
    private:
        edm::InputTag _genParticleInputTag;
        edm::EDGetTokenT<edm::View<reco::GenParticle>> _genParticleToken;

		TRandom3 rng;
        bool filter(edm::Event& iEvent, const edm::EventSetup& iSetup) override;
        
        std::vector<int> acceptedPerClass;
        std::vector<int> totalPerClass;
        
        std::vector<std::vector<int>> pdgIdComb;

    public:
        explicit GenLeptonPairFilter(const edm::ParameterSet&);
        ~GenLeptonPairFilter();

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
        
        const reco::GenParticle* findLeptonFromTau(const reco::GenParticle* genParticle)
        {
            int absId = abs(genParticle->pdgId());
            if (absId==11 or absId==13) return genParticle;
            for (const auto daughter: genParticle->daughterRefVector())
            {
                if (daughter->isDirectPromptTauDecayProductFinalState())
                {
                    const reco::GenParticle* particle = findLeptonFromTau(daughter.get());
                    if (particle) return particle;
                }
            }
            return nullptr;
        }

};


//
// constructors and destructor

//
GenLeptonPairFilter::GenLeptonPairFilter(const edm::ParameterSet& iConfig):
    _genParticleInputTag(iConfig.getParameter<edm::InputTag>("srcGenParticles")),
    _genParticleToken(consumes<edm::View<reco::GenParticle>>(_genParticleInputTag)),
    rng(time(NULL)),
    acceptedPerClass(8,0),
    totalPerClass(8,0),
    pdgIdComb(6,std::vector<int>(6,0))
{
}


GenLeptonPairFilter::~GenLeptonPairFilter()
{
    std::cout<<"========== genlepton filter result ================"<<std::endl;
    for (size_t i=0; i <totalPerClass.size();++i)
    {
        std::cout<<i<<": "<<acceptedPerClass[i]<<"/"<<totalPerClass[i]<<" ("<<(100.*acceptedPerClass[i]/totalPerClass[i])<<"%)"<<std::endl;
    }
    std::cout<<"----------- total events per lepton FS ------------"<<std::endl;
    printf("   ");
    for (size_t i=0; i <pdgIdComb.size();++i)
    {
       printf(" %5lu",i+11);
    }
    printf("\n");
    for (size_t i=0; i <pdgIdComb.size();++i)
    {   
        printf("%2lu:",i+11);
        for (size_t j=0; j <pdgIdComb[i].size();++j)
        {
            printf(" %5i",pdgIdComb[i][j]);
        }
        printf("\n");
    }
    
    std::cout<<"==================================================="<<std::endl;
}




// ------------ method called to produce the data  ------------
bool
GenLeptonPairFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    edm::Handle<edm::View<reco::GenParticle>> genParticleCollection;
    iEvent.getByToken(_genParticleToken, genParticleCollection);

    const reco::GenParticle* l1 = nullptr;
    const reco::GenParticle* l2 = nullptr;

    for (const auto& genParticle: *genParticleCollection)
    {
        if (genParticle.isHardProcess())
        {
            int absId = abs(genParticle.pdgId());
            if (absId>=11 and absId<=16)
            {
                int motherId = abs(genParticle.mother(0)->pdgId());
                if (motherId==9900012 or motherId==9990012)
                {
                    if (l2) std::cout<<"WARNING - l2 already found"<<std::endl;
                    l2 = &genParticle;
                }
                else
                {
                    if (l1) std::cout<<"WARNING - l1 already found"<<std::endl;
                    l1 = &genParticle;
                }
            }
        }
    }
    
    
    if (not l1)
    {
        std::cout<<"ERROR - l1 not found"<<std::endl;
        return false;
    }
    if (not l2)
    {
        std::cout<<"ERROR - l2 not found"<<std::endl;
        return false;
    }
    
    int l1Id = abs(l1->pdgId());
    int l2Id = abs(l2->pdgId());
    
    pdgIdComb[l1Id-11][l2Id-11] += 1;
    
    //std::cout<<l1Id<<","<<l2Id<<std::endl;
    
    int classId = 0; //neutrino FS
    if (l1Id==11 and l2Id==11) classId=1; //ee
    if (l1Id==13 and l2Id==13) classId=2; //mumu
    if (l1Id==15 and l2Id==15) classId=3; //tautau
    if ((l1Id==11 and l2Id==13) or (l1Id==13 and l2Id==11)) classId=4; //emu
    if ((l1Id==11 and l2Id==15) or (l1Id==15 and l2Id==11)) classId=5; //etau
    if ((l1Id==15 and l2Id==13) or (l1Id==13 and l2Id==15)) classId=6; //mutau
    
    //std::cout<<l1Id;
    
    if (l1Id==15)
    {
        l1 = findLeptonFromTau(l1);
        //if (l1) std::cout<<"->"<<abs(l1->pdgId());
        if (not l1) classId = 7; //hadronic tau
    }
    //std::cout<<", "<<l2Id;
    if (l2Id==15)
    {
        l2 = findLeptonFromTau(l2);
        //if (l2) std::cout<<"->"<<abs(l2->pdgId());
        if (not l2) classId = 7; //hadronic tau
    }
    //std::cout<<std::endl;

    totalPerClass[classId]+=1;
    
    float maxPt = std::max<float>(l1?l1->pt():0.f,l2?l2->pt():0.f);
    if (maxPt<20.) return false;
    
    int minClass =  std::distance(acceptedPerClass.begin(), std::min_element(acceptedPerClass.begin(),acceptedPerClass.end()));
    
    
    if (acceptedPerClass[classId]<=(acceptedPerClass[minClass]+2))
    {
        acceptedPerClass[classId]+=1;
        return true;
    }
    
    
    
    return false;
}



// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GenLeptonPairFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}



//define this as a plug-in
DEFINE_FWK_MODULE(GenLeptonPairFilter);
