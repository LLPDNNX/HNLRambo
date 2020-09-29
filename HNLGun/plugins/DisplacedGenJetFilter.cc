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

namespace
{
    //numpy.concatenate([numpy.linspace(10,20,3),numpy.logspace(math.log10(25),3,18)])
	constexpr std::array<float,21> ptBins{{
          10.        ,    15.        ,    20.        ,    25.        ,
          31.05838232,    38.58492449,    47.93541347,    59.55185592,
          73.98337237,    91.91215457,   114.18571346,   141.85694176,
         176.23388528,   218.94157546,   271.99884626,   337.9137663 ,
         419.8021978 ,   521.53508632,   647.92144416,   804.93567703,
        1000.
	}};

	//numpy.logspace(-3,3,13)
	constexpr std::array<float,13> displacementBins{{
	    1.00000000e-03,   3.16227766e-03,   1.00000000e-02,
        3.16227766e-02,   1.00000000e-01,   3.16227766e-01,
        1.00000000e+00,   3.16227766e+00,   1.00000000e+01,
        3.16227766e+01,   1.00000000e+02,   3.16227766e+02,
        1.00000000e+03
	}};

	constexpr std::array<float,8> classBins{{
			0.5,
			1.5, //e
			2.5, //mu
			3.5, //tau
			4.5, //eq
			5.5, //muq
			6.5, //tauq
			7.5 //q
	}};
}

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

				static bool isDaughter(const reco::Candidate* genParticle, int pdgId)
				{
						if (abs(genParticle->pdgId())==pdgId) return true;
						if (genParticle->numberOfMothers()>0)
						{
								return isDaughter(genParticle->mother(0),pdgId);
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

				static int getClassBin(int nE, int nMu, int nTau, int nQ)
				{
						int classBin = 0;
						if (nE>0)
								if (nQ>0) classBin = 4;
								else classBin = 1;
						else if (nMu>0)
								if (nQ>0) classBin = 5;
								else classBin = 2;
						else if (nTau>0)
								if (nQ>0) classBin = 6;
								else classBin = 3;
						else if (nQ>0)
								classBin = 7;
						return classBin;
				}

        edm::InputTag _genParticleInputTag;
        edm::EDGetTokenT<edm::View<reco::GenParticle>> _genParticleToken;

        edm::InputTag _genJetInputTag;
        edm::EDGetTokenT<edm::View<reco::GenJet>> _genJetToken;

  			TRandom3 rng;

        TH3F ptDisplacementHist;
        std::vector<size_t> njets;

				const std::unordered_set<int> _daughterIds;

        bool acceptJet(double pt, double displacement, int classBin);
        void fillJet(double pt, double displacement, int classBin);

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
    ptDisplacementHist(
        "pt","",
        ptBins.size()-1,ptBins.data(),
        displacementBins.size()-1,displacementBins.data(),
				classBins.size()-1,classBins.data()
    ),
    njets(10,0),
		_daughterIds({1,2,3,4,5,11,13,15})
{
}


DisplacedGenJetFilter::~DisplacedGenJetFilter()
{
    for (int i = 0; i < 110; ++i) std::cout<<"=";
    std::cout<<std::endl;
    std::cout<<"total jets: "<<ptDisplacementHist.GetEntries()<<std::endl;
    std::cout<<"max jets: "<<ptDisplacementHist.GetMaximum()<<std::endl;

    std::cout<<"           ";
    for (int dbin = 0; dbin < ptDisplacementHist.GetNbinsY(); ++dbin)
    {
        printf("%7.1e ",ptDisplacementHist.GetYaxis()->GetBinCenter(dbin+1));
    }
    std::cout<<std::endl;
    for (int i = 0; i < 110; ++i) std::cout<<"-";
    std::cout<<std::endl;
		TH1D* projectedPtDisplacementHist = ptDisplacementHist.ProjectionZ();

    for (int ptbin = 0; ptbin < ptDisplacementHist.GetNbinsX(); ++ptbin)
    {
				for (int iclass = 1; iclass < 8; ++iclass)
				{
		        printf("pt=%5.1f, c=%1i: ",ptDisplacementHist.GetXaxis()->GetBinCenter(ptbin+1),iclass);
		        for (int dbin = 0; dbin < ptDisplacementHist.GetNbinsY(); ++dbin)
		        {
		            printf("%7.0f ",ptDisplacementHist.GetBinContent(ptbin+1,dbin+1,iclass));
		        }
		        std::cout<<std::endl;
				}
    }
    for (int i = 0; i < 110; ++i) std::cout<<"-";
    std::cout<<std::endl;

		for (int i = 0; i < projectedPtDisplacementHist->GetNbinsX(); ++i)
		{
				printf("c=%1i: %7.0f\n",i+1,projectedPtDisplacementHist->GetBinContent(i+1));
		}
		for (int i = 0; i < 110; ++i) std::cout<<"-";
    std::cout<<std::endl;

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
    for (int i = 0; i < 110; ++i) std::cout<<"=";
    std::cout<<std::endl;
}



bool DisplacedGenJetFilter::acceptJet(double pt, double displacement, int classBin)
{
    //try to achieve uniform sampling in log(pt) and log(displacement)
    int ptBin = ptDisplacementHist.GetXaxis()->FindBin(pt);
    int displacementBin = ptDisplacementHist.GetYaxis()->FindBin(displacement);

		if (classBin==0) return false;

    double binContent = ptDisplacementHist.GetBinContent(ptBin,displacementBin,classBin);
    double maximum = ptDisplacementHist.GetMaximum();
    if (ptDisplacementHist.GetEntries()!=0 and maximum<rng.Gaus(binContent+1,1))
    {
        return false;
    }
    return true;
}

void DisplacedGenJetFilter::fillJet(double pt, double displacement, int classBin)
{
    ptDisplacementHist.Fill(pt,displacement,classBin);
}


// ------------ method called to produce the data  ------------
bool
DisplacedGenJetFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    edm::Handle<edm::View<reco::GenParticle>> genParticleCollection;
    iEvent.getByToken(_genParticleToken, genParticleCollection);

    edm::Handle<edm::View<reco::GenJet>> genJetCollection;
    iEvent.getByToken(_genJetToken, genJetCollection);

    std::vector<DisplacedVertex> llpVertices;
		std::vector<const reco::GenParticle*> llpDaughters;
    for (const auto& genParticle: *genParticleCollection)
    {
        if (genParticle.pdgId()==9900012 and genParticle.numberOfDaughters()>0)
        {
            double distance = std::sqrt((genParticle.vertex()-genParticle.daughterRef(0)->vertex()).mag2());
            if (distance>1e-3 and distance<1e3) //10mum .. 10m
            {
                llpVertices.emplace_back(
                    &genParticle,
                    genParticle.daughterRef(0)->vertex(),
                    distance
                );
            }
        }
				if (_daughterIds.find(abs(genParticle.pdgId()))!=_daughterIds.cend())
				{
						if (genParticle.numberOfDaughters()>0 and genParticle.daughterRef(0)->pdgId()==genParticle.pdgId()) continue;
						if (isDaughter(&genParticle,9900012))
						{
								llpDaughters.push_back(&genParticle);
						}
				}
    }

    size_t nacceptedJets = 0;
    std::vector<double> acceptedJetPts;
    std::vector<double> acceptedJetDisplacements;
		std::vector<int> acceptedJetClasses;

    for (const auto& genJet: *genJetCollection)
    {
				int nE = 0;
				int nMu = 0;
				int nTau = 0;
				int nQ = 0;

				reco::Candidate::LorentzVector leptonP4Sum(0,0,0,0);

				for (const auto llpDaughter: llpDaughters)
				{
						double dR = reco::deltaR(*llpDaughter,genJet);
						if (dR<0.3) //use 0.3 instead of 0.4 to exclude cases where decay products of close-by jet overlaps
						{
								if (abs(llpDaughter->pdgId())==11)
								{
										nE+=1;
										leptonP4Sum+=llpDaughter->p4();
								}
								if (abs(llpDaughter->pdgId())==13)
								{
										nMu+=1;
										leptonP4Sum+=llpDaughter->p4();
								}
								if (abs(llpDaughter->pdgId())==15)
								{
										nTau+=1;
										leptonP4Sum+=llpDaughter->p4();
								}
								if (abs(llpDaughter->pdgId())<6) nQ+=1;
						}
				}

        if (genJet.pt()<10. or (genJet.p4()-leptonP4Sum).pt()<10. or std::fabs(genJet.eta())>2.4)
        {
            continue;
        }
        std::vector<reco::Candidate::LorentzVector> momentaPerVertex(llpVertices.size(),reco::Candidate::LorentzVector(0,0,0,0));
        reco::Candidate::LorentzVector restMomenta(0,0,0,0);
        for (const auto& genParticle: genJet.getGenConstituents())
        {
            bool foundVertex = false;
            for (size_t ivertex = 0; ivertex < llpVertices.size(); ++ivertex)
            {
                double distance = std::sqrt((llpVertices[ivertex].vertex-correctedDisplacement(*genParticle)).mag2());
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
        for (size_t ivertex = 0; ivertex < llpVertices.size(); ++ivertex)
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
						int classBin = getClassBin(nE, nMu, nTau, nQ);
            if (acceptJet(genJet.pt(),llpVertices[vertexIndex].displacement,classBin))
            {
                nacceptedJets+=1;
                acceptedJetPts.push_back(genJet.pt());
                acceptedJetDisplacements.push_back(llpVertices[vertexIndex].displacement);
								acceptedJetClasses.push_back(classBin);
            }
            else
            {
                return false; //reject whole event
            }
        }
    }

    //if (nacceptedJets<1) return false;
    if (nacceptedJets<(0.5+std::fabs(rng.Gaus(0.,1.5)))) return false;

    for (size_t i = 0; i < acceptedJetPts.size(); ++i)
    {
        fillJet(acceptedJetPts[i],acceptedJetDisplacements[i],acceptedJetClasses[i]);
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
