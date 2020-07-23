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

#include "LLPReco/DataFormats/interface/DisplacedGenVertex.h"

#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/angle.h"

#include <ctime>

#include "TH2F.h"
#include "TRandom3.h"

struct DisplacedVertex
{
    const reco::GenParticle* motherLLP;
    reco::Candidate::Point vertex;
    double displacement;

    DisplacedVertex(
        const reco::GenParticle* motherLLP,
        reco::Candidate::Point vertex,
        double displacement
    ):
        motherLLP(motherLLP),
        vertex(vertex),
        displacement(displacement)
    {
    }
};

class DisplacedGenJetFilter:
    public edm::one::EDFilter<>
    
{
    public:
        constexpr static double MIN_DISPLACEMENT = 1e-10;
    private:
        static double distance(const reco::Candidate::Point& p1, const reco::Candidate::Point& p2)
        {
            return std::sqrt((p1-p2).mag2());
        }
        
        static bool ignoreDisplacement(const reco::Candidate& genParticle)
        {
            
            int absPdgId = std::abs(genParticle.pdgId());
            if (absPdgId==111) //ignore pi0 to gamma gamma
            {
                return true;
            }
            
            return false;
        }
        
        static reco::Candidate::Point correctedDisplacement(const reco::Candidate& genParticle)
        {
            //return mother vertex if displacement is ignored
            if (genParticle.mother() and ignoreDisplacement(*genParticle.mother()) and  (distance(genParticle.mother()->vertex(),genParticle.vertex())>MIN_DISPLACEMENT))
            {
                return correctedDisplacement(*genParticle.mother()); //call recursively
            }
            return genParticle.vertex();
        }
        
        edm::InputTag _genParticleInputTag;
        edm::EDGetTokenT<edm::View<reco::GenParticle>> _genParticleToken;
        
        edm::InputTag _genJetInputTag;
        edm::EDGetTokenT<edm::View<reco::GenJet>> _genJetToken;
        
    	TRandom3 rng;
        TH2F ptDisplacementHist;
        std::vector<size_t> njets;
        
        bool acceptJet(double pt, double displacement);
        void fillJet(double pt, double displacement);

        bool filter(edm::Event& iEvent, const edm::EventSetup& iSetup) override;

    public:
        explicit DisplacedGenJetFilter(const edm::ParameterSet&);
        ~DisplacedGenJetFilter();

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

};


//
// constructors and destructor

//
DisplacedGenJetFilter::DisplacedGenJetFilter(const edm::ParameterSet& iConfig):
    _genParticleInputTag(iConfig.getParameter<edm::InputTag>("srcGenParticles")),
    _genParticleToken(consumes<edm::View<reco::GenParticle>>(_genParticleInputTag)),
    _genJetInputTag(iConfig.getParameter<edm::InputTag>("srcGenJets")),
    _genJetToken(consumes<edm::View<reco::GenJet>>(_genJetInputTag)), 
    rng(time(NULL)),
    ptDisplacementHist("pt","",20,1,3,10,-3,3),
    njets(10,0)
{
}


DisplacedGenJetFilter::~DisplacedGenJetFilter()
{
    std::cout<<"total jets: "<<ptDisplacementHist.GetEntries()<<std::endl;
    std::cout<<"max jets: "<<ptDisplacementHist.GetMaximum()<<std::endl;

    std::cout<<"           ";
    for (int dbin = 0; dbin < ptDisplacementHist.GetNbinsY(); ++dbin)
    {
        printf("%5.1e ",std::pow(10.,ptDisplacementHist.GetYaxis()->GetBinCenter(dbin+1)));
    }
    std::cout<<std::endl;
    std::cout<<"--------------------------------------------------------"<<std::endl;
    for (int ptbin = 0; ptbin < ptDisplacementHist.GetNbinsX(); ++ptbin)
    {
        printf("pt=%5.1f: ",std::pow(10.,ptDisplacementHist.GetXaxis()->GetBinCenter(ptbin+1)));
        for (int dbin = 0; dbin < ptDisplacementHist.GetNbinsY(); ++dbin)
        {
            printf("%5.0f ",ptDisplacementHist.GetBinContent(ptbin+1,dbin+1));
        }
        std::cout<<std::endl;
    }
    for (size_t i = 0; i < njets.size(); ++i)
    {
        printf("%3i ",int(i));
    }
    std::cout<<std::endl;
    for (size_t i = 0; i < njets.size(); ++i)
    {
        printf("%3i ",int(njets[i]));
    }
    std::cout<<std::endl;
}

bool DisplacedGenJetFilter::acceptJet(double pt, double displacement)
{
    //try to achieve uniform sampling in log(pt) and log(displacement)
    int ptBin = ptDisplacementHist.GetXaxis()->FindBin(std::log10(pt));
    int displacementBin = ptDisplacementHist.GetYaxis()->FindBin(std::log10(displacement));
    double binContent = ptDisplacementHist.GetBinContent(ptBin,displacementBin);
    double maximum = ptDisplacementHist.GetMaximum();
    if (ptDisplacementHist.GetEntries()!=0 and maximum<rng.Gaus(binContent+1,1))
    {
        return false;
    }
    return true;
}

void DisplacedGenJetFilter::fillJet(double pt, double displacement)
{
    ptDisplacementHist.Fill(std::log10(pt),std::log10(displacement));
}


// ------------ method called to produce the data  ------------
bool
DisplacedGenJetFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    edm::Handle<edm::View<reco::GenParticle>> genParticleCollection;
    iEvent.getByToken(_genParticleToken, genParticleCollection);
    
    edm::Handle<edm::View<reco::GenJet>> genJetCollection;
    iEvent.getByToken(_genJetToken, genJetCollection);
    
    std::vector<DisplacedVertex> vertices;
    for (const auto& genParticle: *genParticleCollection)
    {
        if (genParticle.pdgId()==9900012 and genParticle.numberOfDaughters()>0)
        {
            double distance = std::sqrt((genParticle.vertex()-genParticle.daughterRef(0)->vertex()).mag2());
            if (distance>1e-3 and distance<1e3) //10mum .. 10m
            {
                vertices.emplace_back(
                    &genParticle,
                    genParticle.daughterRef(0)->vertex(),
                    distance
                );
            }
            
        }
    }
    
    size_t nacceptedJets = 0;
    std::vector<double> acceptedJetPts;
    std::vector<double> acceptedJetDisplacements;
    
    for (const auto& genJet: *genJetCollection)
    {
        if (genJet.pt()<10. or std::fabs(genJet.eta())>2.4)
        {
            continue;
        }
        std::vector<reco::Candidate::LorentzVector> momentaPerVertex(vertices.size(),reco::Candidate::LorentzVector(0,0,0,0));
        reco::Candidate::LorentzVector restMomenta(0,0,0,0);
        for (const auto& genParticle: genJet.getGenConstituents())
        {
            bool foundVertex = false;
            for (size_t ivertex = 0; ivertex < vertices.size(); ++ivertex)
            {
                double distance = std::sqrt((vertices[ivertex].vertex-correctedDisplacement(*genParticle)).mag2());
                if (distance<MIN_DISPLACEMENT) 
                {
                    momentaPerVertex[ivertex]+=genParticle->p4();
                    foundVertex = true;
                    break;
                }
            }
            if (not foundVertex) restMomenta+=genParticle->p4();
            
        }
        double maxSharedVertexFraction = restMomenta.Vect().Dot(genJet.p4().Vect())/(genJet.p()*genJet.p());
        int vertexIndex = -1;
        for (size_t ivertex = 0; ivertex < vertices.size(); ++ivertex)
        {
            double sharedVertexFraction = momentaPerVertex[ivertex].Vect().Dot(genJet.p4().Vect())/(genJet.p()*genJet.p());
            if (sharedVertexFraction>maxSharedVertexFraction)
            {
                maxSharedVertexFraction = sharedVertexFraction;
                vertexIndex = ivertex;
            }
        }
        
        //most of the genjet momentum comes from an LLP decay vertex
        if (vertexIndex>=0)
        {
            if (acceptJet(genJet.pt(),vertices[vertexIndex].displacement))
            {
                nacceptedJets+=1;
                acceptedJetPts.push_back(genJet.pt());
                acceptedJetDisplacements.push_back(vertices[vertexIndex].displacement);
            }
            else
            {
                return false; //reject whole event
            }
        }
    }
    
    //if (nacceptedJets<1) return false; 
    if (nacceptedJets<(0.1+std::fabs(rng.Gaus(0.,0.5)))) return false;
    
    for (size_t i = 0; i < acceptedJetPts.size(); ++i)
    {
        fillJet(acceptedJetPts[i],acceptedJetDisplacements[i]);
    }
    if (nacceptedJets<njets.size()) njets[nacceptedJets]+=1;
    
    return true;
    
}



// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DisplacedGenJetFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}



//define this as a plug-in
DEFINE_FWK_MODULE(DisplacedGenJetFilter);


