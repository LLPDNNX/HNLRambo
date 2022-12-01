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
        
        std::unordered_map<int,int> acceptedPerClass;
        std::unordered_map<int,int> totalPerClass;
        

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
        
        float filterProb(int l1Id, int l2Id) const
        {
            if ((l1Id>=11 and l1Id<=16) and (l2Id==11 or l2Id==13)) return 0.01;
            if ((l1Id==11 or l1Id==13) and (l2Id>=11 and l2Id<=16)) return 0.01;
            
            
            if ((l1Id>=11 and l1Id<=16) and (l2Id==12 or l2Id==14 or l2Id==16)) return 0.002;
            if ((l1Id==12 or l1Id==14 or l1Id==16) and (l2Id>=11 and l2Id<=16)) return 0.002;
            
            if ((l1Id>=11 and l1Id<=16) and l2Id==100) return 0.002;
            if ((l1Id>=11 and l1Id<=16) and (l2Id==111 or l2Id==113)) return 0.1;
            
            if ((l1Id==100) and (l2Id==11 or l2Id==13)) return 0.01;
            if ((l1Id==100) and (l2Id==12 or l2Id==14 or l2Id==16)) return 0.002;
            
            if ((l1Id==111 or l1Id==113) and (l2Id==11 or l2Id==13)) return 0.1;
            if ((l1Id==111 or l1Id==113) and (l2Id==12 or l2Id==14 or l2Id==16)) return 0.02;
            
            if ((l1Id>=100) and l2Id==100) return 0.05;
            
            return 1.0;
        }

};


//
// constructors and destructor

//
GenLeptonPairFilter::GenLeptonPairFilter(const edm::ParameterSet& iConfig):
    _genParticleInputTag(iConfig.getParameter<edm::InputTag>("srcGenParticles")),
    _genParticleToken(consumes<edm::View<reco::GenParticle>>(_genParticleInputTag)),
    rng(time(NULL))
{
}


GenLeptonPairFilter::~GenLeptonPairFilter()
{
    std::vector<int> ids;
    for (const auto& idPair: totalPerClass)
    {
        ids.push_back(idPair.first);
    }
    std::sort(ids.begin(), ids.end());
    
    int acceptedTotal = 0;
    
    
    std::cout<<"========== genlepton filter result ================"<<std::endl;
    for (int id: ids)
    {
        acceptedTotal+=acceptedPerClass[id];
        printf("%9i: %5i/%5i = %6.2f%%\n",id,acceptedPerClass[id],totalPerClass[id],100.*acceptedPerClass[id]/totalPerClass[id]);
    }
    std::cout<<"=========== fraction / total accepted ============="<<std::endl;
    printf("Total accepted: %i\n",acceptedTotal);
    for (int id: ids)
    {
        printf("%9i: %5i => %6.2f%%\n",id,acceptedPerClass[id],100.*acceptedPerClass[id]/acceptedTotal);
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
    
    

    if (l1Id==15)
    {
        auto l1TauDecay = findLeptonFromTau(l1);
        if (not l1TauDecay) 
        {
            l1Id=100; //hadronic tau
        }
        else 
        {
            l1 = l1TauDecay;
            l1Id=100+abs(l1->pdgId());
        }
    }
    if (l2Id==15)
    {
        auto l2TauDecay = findLeptonFromTau(l2);
        if (not l2TauDecay) 
        {
            l2Id=100; //hadronic tau
        }
        else 
        {
            l2 = l2TauDecay;
            l2Id=100+abs(l2->pdgId());
        }
    }

    totalPerClass[l1Id*1000+l2Id]+=1;
    
    float lMaxPt = 0.f;
    //int lMaxId = 0;
    //int lMax = 0;
    //l1 is e or mu
    if (l1Id==11 or l1Id==13 or l1Id==111 or l1Id==113)
    {
        //lMaxId=l1Id;
        //lMax = 1;
        lMaxPt = l1->pt();
        //l2 is e or mu; take pT max of l1/l2
        if ((l2Id==11 or l2Id==13 or l2Id==111 or l2Id==113) and l1->pt()<l2->pt())
        {
            //lMax = 2;
            //lMaxId=l2Id;
            lMaxPt = l2->pt();
        }
    }
    //no l1; only l2 is e or mu
    else if (l2Id==11 or l2Id==13 or l2Id==111 or l2Id==113)
    {
        //lMax = 2;
        //lMaxId=l2Id;
        lMaxPt = l2->pt();
    }
    
    
    if (lMaxPt<20.) return false;
    if (rng.Uniform()>filterProb(l1Id,l2Id)) return false;
    
    acceptedPerClass[l1Id*1000+l2Id]+=1;
    return true;
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
