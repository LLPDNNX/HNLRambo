#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDProducer.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"

#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/InputSourceMacros.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Sources/interface/ProducerSourceBase.h"

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

#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/MTwistEngine.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"

#include "HNLRambo/Rambo/interface/Rambo.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "TLorentzVector.h"

#include <ctime>
#include <vector>

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

class ParticleSystem;
typedef std::shared_ptr<ParticleSystem> ParticleSystemPtr;

class ParticleSystem
{
    protected:
        int pdgId_;
        double ctau_;
        int colorTag_;
        TLorentzVector momentum_;
        std::vector<ParticleSystemPtr> subsystems_;

        ParticleSystem(
            int pdgId,
            TLorentzVector momentum = TLorentzVector(0.,0.,0.,0.)
        ):
            pdgId_(pdgId),
            ctau_(0),
            colorTag_(0),
            momentum_(momentum)
        {
        }

    public:
        static ParticleSystemPtr create(int pdgId, double mass=1e-6)
        {
            TLorentzVector momentum;
            momentum.SetXYZM(0,0,0,mass);
            return std::shared_ptr<ParticleSystem>(new ParticleSystem(pdgId,momentum));
        }

        const ParticleSystemPtr& add(const ParticleSystemPtr& ParticleSystem)
        {
            subsystems_.push_back(ParticleSystem);
            return ParticleSystem;
        }

        inline int getPdgId() const
        {
            return pdgId_;
        }

        inline void setCtau(double ctau)
        {
            ctau_ = ctau;
        }

        inline void setColorTag(int colorTag)
        {
            colorTag_ = colorTag;
        }

        inline double mass() const
        {
            return momentum_.Mag();
        }

        inline void setP4(const TLorentzVector& momentum)
        {
            momentum_ = momentum;
        }

        inline TLorentzVector getP4() const
        {
            return momentum_;
        }

        void generateRandomMomenta(const std::shared_ptr<CLHEP::HepRandomEngine>& rngEngine_)
        {
            if (subsystems_.size()==0)
            {
                return;
            }
            Rambo rambo(CLHEP::RandFlat::shootInt(rngEngine_.get(),1e12));

            std::vector<double> masses(subsystems_.size(),0);
            for (size_t i = 0; i < masses.size(); ++i)
            {
                masses[i] = subsystems_[i]->mass();
            }

            double weight = 0;
            std::vector<Rambo4Vector> randomMomenta = rambo.generate(
                this->mass(), //on-shell
                masses,
                weight
            );


            //randomize rotation
            constexpr double PI = std::atan(1.0)*4;
            double rotateX = CLHEP::RandFlat::shoot(rngEngine_.get(),0,2*PI);
            double rotateY = CLHEP::RandFlat::shoot(rngEngine_.get(),0,2*PI);
            double rotateZ = CLHEP::RandFlat::shoot(rngEngine_.get(),0,2*PI);

            TVector3 boostVector = momentum_.BoostVector();

            for (size_t i = 0; i < randomMomenta.size(); ++i)
            {
                TLorentzVector randomMomentum = randomMomenta[i].toLorentzVector();
                randomMomentum.RotateX(rotateX);
                randomMomentum.RotateY(rotateY);
                randomMomentum.RotateZ(rotateZ);

                randomMomentum.Boost(boostVector);

                subsystems_[i]->setP4(randomMomentum);
                subsystems_[i]->generateRandomMomenta(rngEngine_); //call recursively
            }
        }

        int nParticles() const
        {
            int n = 1;
            for (size_t i = 0; i < subsystems_.size(); ++i)
            {
                n += subsystems_[i]->nParticles();
            }
            return n;
        }

        void fillLHE(lhef::HEPEUP& hepeup, size_t& index, const std::array<int,2>& motherIndices=std::array<int,2>{{0,0}})
        {
            hepeup.IDUP[index] = pdgId_; // particle ID
            hepeup.ISTUP[index] = subsystems_.size()>0 ? 2 : 1; //status code (-1 incomming, +1 outgoing, +2 intermediate resonance)

            hepeup.MOTHUP[index] = {motherIndices[0],motherIndices[1]}; // mother indices (0 = no mother)
            hepeup.ICOLUP[index] = {
                pdgId_>0 ? colorTag_ : 0,
                pdgId_<0 ? colorTag_ : 0,
            }; //color tag to indicate color flow ({particle, antiparticle})


            hepeup.PUP[index][0]= momentum_.Px(); //px
            hepeup.PUP[index][1]= momentum_.Py(); //py
            hepeup.PUP[index][2]= momentum_.Pz(); //pz
            hepeup.PUP[index][3]= momentum_.E(); //E
            hepeup.PUP[index][4]= momentum_.Mag(); //mass
            hepeup.VTIMUP[index] = ctau_; //ctau in mm
            hepeup.SPINUP[index] = 0.; //helicity

            //set fermions to be always left-handed; antifermions right-handed
            if (pdgId_>=1 and pdgId_<=16) hepeup.SPINUP[index] = pdgId_>0 ? -1 : 1;


            int motherIndex = index+1; //counts from 1

            index+=1;

            for (size_t i = 0; i < subsystems_.size(); ++i)
            {
                subsystems_[i]->fillLHE(hepeup,index,{{motherIndex,motherIndex}});
            }
        }
};




class HNLGun:
    public edm::ProducerSourceBase
{
    private:
        const double cme = 13000;

        const std::array<double,3> leptonMasses{{0.511e-3,0.105,1.776}};
        const std::array<double,4> quarkMasses{{4.8e-3,2.3e-3,95.e-3,1.28}};

        std::shared_ptr<CLHEP::HepRandomEngine> rngEngine_;
        Rambo rambo_;

    public:
        explicit HNLGun(const edm::ParameterSet& iConfig,  const edm::InputSourceDescription& desc):
            edm::ProducerSourceBase(iConfig, desc, false),
            //source is not allowed to use the RandomNumberService so just create an engine
            rngEngine_(new CLHEP::MTwistEngine(12345))//time(NULL)))
        {


            produces<LHERunInfoProduct,edm::InRun>();
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

        bool decayQQ(std::vector<ParticleSystemPtr>& displacedParticles, int colorTag, ParticleSystemPtr llp)
        {
            int quarkType1 = CLHEP::RandFlat::shootInt(rngEngine_.get(),0,quarkMasses.size());
            int quarkType2 = (quarkType1+2*CLHEP::RandFlat::shootInt(rngEngine_.get(),0,quarkMasses.size()))%quarkMasses.size();

            int quarkId1 = (quarkType1+1);
            int quarkId2 = (quarkType2+1);

            if (2*(quarkMasses[quarkType1]+quarkMasses[quarkType2])>llp->mass())
            {
                return false;
            }

            auto p1 = llp->add(ParticleSystem::create(quarkId1,quarkMasses[quarkType1]));
            p1->setColorTag(100+colorTag);
            displacedParticles.push_back(p1);
            auto p2 = llp->add(ParticleSystem::create(-quarkId2,quarkMasses[quarkType2]));
            p2->setColorTag(100+colorTag);
            displacedParticles.push_back(p2);

            return true;
        }

        bool decayQQNu(std::vector<ParticleSystemPtr>& displacedParticles, int colorTag, ParticleSystemPtr llp)
        {
            int quarkType1 = CLHEP::RandFlat::shootInt(rngEngine_.get(),0,quarkMasses.size());
            int quarkType2 = (quarkType1+2*CLHEP::RandFlat::shootInt(rngEngine_.get(),0,quarkMasses.size()))%quarkMasses.size();

            int quarkId1 = (quarkType1+1);
            int quarkId2 = (quarkType2+1);

            if (2*(quarkMasses[quarkType1]+quarkMasses[quarkType2])>llp->mass())
            {
                return false;
            }

            //auto zboson = llp->add(ParticleSystem::create(23,3*quarkMasses[quarkType1]));

            auto p1 = llp->add(ParticleSystem::create(quarkId1,quarkMasses[quarkType1]));
            p1->setColorTag(100+colorTag);
            displacedParticles.push_back(p1);

            auto p2 = llp->add(ParticleSystem::create(-quarkId2,quarkMasses[quarkType2]));
            p2->setColorTag(100+colorTag);
            displacedParticles.push_back(p2);

            //add neutrino
            llp->add(ParticleSystem::create(12));

            return true;
        }

        bool decayQQL(std::vector<ParticleSystemPtr>& displacedParticles, int colorTag, ParticleSystemPtr llp, int leptonGeneration)
        {
            int leptonFlavor = std::copysign(
                leptonGeneration*2+11,
                CLHEP::RandFlat::shoot(rngEngine_.get(),-1,1)
            );
            double leptonMass = leptonMasses[leptonGeneration];

            if (2*leptonMass>llp->mass()) return false;

            auto lepton = llp->add(ParticleSystem::create(leptonFlavor,leptonMass));
            displacedParticles.push_back(lepton);

            int quarkTypeDown = (2*CLHEP::RandFlat::shootInt(rngEngine_.get(),0,quarkMasses.size()))%quarkMasses.size();
            int quarkTypeUp = (1+2*CLHEP::RandFlat::shootInt(rngEngine_.get(),0,quarkMasses.size()))%quarkMasses.size();

            int quarkIdDown = (quarkTypeDown+1);
            int quarkIdUp = (quarkTypeUp+1);

            //match EM charge
            if (leptonFlavor>0)
            {
                //auto wboson = llp->add(ParticleSystem::create(24,3*quarkMasses[quarkTypeUp]));

                auto p1 = llp->add(ParticleSystem::create(quarkIdUp,quarkMasses[quarkTypeUp]));
                p1->setColorTag(100+colorTag);
                displacedParticles.push_back(p1);

                auto p2 = llp->add(ParticleSystem::create(-quarkIdDown,quarkMasses[quarkTypeDown]));
                p2->setColorTag(100+colorTag);
                displacedParticles.push_back(p2);
            }
            else
            {
                //auto wboson = llp->add(ParticleSystem::create(-24,10*quarkMasses[quarkTypeUp]));

                auto p1 = llp->add(ParticleSystem::create(-quarkIdUp,quarkMasses[quarkTypeUp]));
                p1->setColorTag(100+colorTag);
                displacedParticles.push_back(p1);

                auto p2 = llp->add(ParticleSystem::create(quarkIdDown,quarkMasses[quarkTypeDown]));
                p2->setColorTag(100+colorTag);
                displacedParticles.push_back(p2);
            }

            return true;
        }

        bool decayLL(std::vector<ParticleSystemPtr>& displacedParticles, int colorTag, ParticleSystemPtr llp, int leptonGeneration1, int leptonGeneration2)
        {

            double leptonMass1 = leptonMasses[leptonGeneration1];
            double leptonMass2 = leptonMasses[leptonGeneration2];

            if (2*(leptonMass1+leptonMass2)>llp->mass()) return false;

            int leptonSign = CLHEP::RandFlat::shootInt(rngEngine_.get(),0,2)*2-1;

            auto p1 = llp->add(ParticleSystem::create((leptonGeneration1*2+11)*leptonSign,leptonMass1));
            displacedParticles.push_back(p1);

            auto p2 = llp->add(ParticleSystem::create(-(leptonGeneration2*2+11)*leptonSign,leptonMass2));
            displacedParticles.push_back(p2);

            llp->add(ParticleSystem::create(12));

            return true;
        }

        LHEEventProduct* generate()
        {
            double scale = std::pow(10.,CLHEP::RandGauss::shoot(rngEngine_.get(),2.5,0.25));

            if (scale<60.) return nullptr;
            if (scale>cme) return nullptr;

            double pz = std::copysign(
                std::min(cme/2.1,std::pow(10.,CLHEP::RandGauss::shoot(rngEngine_.get(),2.,0.75))),
                CLHEP::RandFlat::shoot(rngEngine_.get(), -1, 1)
            );

            //take only positive solution
            double x2 = -pz/cme+std::sqrt((pz*pz+scale*scale)/cme/cme);
            double x1 = scale*scale/cme/cme/x2;

            if (x1>0.9 or x2>0.9) return nullptr;
            if (x1<1e-4 or x2<1e-4) return nullptr;


            int nHNL = CLHEP::RandFlat::shootInt(rngEngine_.get(),1,4);
            int nNonLLP = 1+std::abs(CLHEP::RandGauss::shoot(rngEngine_.get(),0,2.5)); //need at least 1 non HNL; otherwise HNL pT = 0
            //int nParticles = nHNL+nNonLLP;


            double massLLP = CLHEP::RandFlat::shoot(rngEngine_.get(), 1., 5)+std::fabs(CLHEP::RandGauss::shoot(rngEngine_.get(),0.,10.)); //in GeV
            double ctauLLP = std::pow(10.,CLHEP::RandFlat::shoot(rngEngine_.get(), -3., 5.)); //in mm

            if (scale<(massLLP*nHNL*1.1)) return nullptr;



            ParticleSystemPtr ps = ParticleSystem::create(1000022,scale);
            ps->setP4(TLorentzVector(0,0,pz,std::sqrt(scale*scale+pz*pz)));
            for (int i = 0; i < (nNonLLP/2); ++i)
            {
                int quarkType1 = CLHEP::RandFlat::shootInt(rngEngine_.get(),0,quarkMasses.size());
                int quarkType2 = (quarkType1+2*CLHEP::RandFlat::shootInt(rngEngine_.get(),0,quarkMasses.size()))%quarkMasses.size();

                int quarkId1 = (quarkType1+1);
                int quarkId2 = (quarkType2+1);

                auto p1 = ps->add(ParticleSystem::create(
                    quarkId1,
                    quarkMasses[quarkType1]
                ));
                p1->setColorTag(i+10);
                auto p2 = ps->add(ParticleSystem::create(
                    -quarkId2,
                    quarkMasses[quarkType2]
                ));
                p2->setColorTag(i+10);
            }

            if (nNonLLP%2==1)
            {
                //add prompt neutrino in case of odd number of particles
                ps->add(ParticleSystem::create(12,1e-6));
            }

            std::vector<ParticleSystemPtr> llps;
            std::vector<ParticleSystemPtr> displacedParticles;
            for (int i = 0; i < nHNL; ++i)
            {
                auto llp = ps->add(ParticleSystem::create(9900012,massLLP));
                llp->setCtau(ctauLLP);
                llps.push_back(llp);


                bool sucess = false;
                do
                {
                    double llpDecayType = CLHEP::RandFlat::shoot(rngEngine_.get(),0,1);

                    if (llpDecayType<0.01) sucess = decayQQ(displacedParticles,i,llp); //1.% QQ
                    else if (llpDecayType<0.01) sucess = decayQQNu(displacedParticles,i,llp); //1.% QQnu
                    else if (llpDecayType<0.05) sucess = decayQQL(displacedParticles,i,llp,0); //3% QE
                    else if (llpDecayType<0.10) sucess = decayQQL(displacedParticles,i,llp,1); //5% QMu
                    else if (llpDecayType<0.13) sucess = decayLL(displacedParticles,i,llp,0,0); //3% EE
                    else if (llpDecayType<0.18) sucess = decayLL(displacedParticles,i,llp,1,1); //5% MUMU
                    else if (llpDecayType<0.55) sucess = decayQQL(displacedParticles,i,llp,2); //37% QTau
                    else if (llpDecayType<0.65) sucess = decayLL(displacedParticles,i,llp,2,0); //10% TAUE
                    else if (llpDecayType<0.75) sucess = decayLL(displacedParticles,i,llp,2,1); //10% TAUMU
                    else if (llpDecayType<1.00) sucess = decayLL(displacedParticles,i,llp,2,2); //25% TAUTAU
                }
                while (not sucess);
            }


            ps->generateRandomMomenta(rngEngine_);

            //require at least one LLP to be within eta acceptance
            //and require one LLP within acceptance to have an energy above 10 GeV
            double minAbsEta = 100.;
            double maxE = 0.;
            for (const auto& llp: llps)
            {
                double absEta = std::fabs(llp->getP4().Eta());
                minAbsEta = std::min(minAbsEta,absEta);
                if (absEta<2.4)
                {
                    maxE = std::max(maxE,std::fabs(llp->getP4().E()));
                }
            }
            //if (minAbsEta>2.4) return nullptr;
            //if (maxE<10.) return nullptr;


            //require at least 50% of all displaced particles to be within acceptance
            int nInAcceptance = 0;
            for (const auto& particle: displacedParticles)
            {
                if (std::fabs(particle->getP4().Eta())<2.4 and particle->getP4().Pt()>10.)
                {
                    nInAcceptance+=1;
                }
            }
            if ((1.*nInAcceptance/displacedParticles.size())<0.5) return nullptr;



            int nTotalParticles = ps->nParticles();

            lhef::HEPEUP hepeup;
            hepeup.resize(nTotalParticles+2); //number of particles
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

            size_t index = 2;
            ps->fillLHE(hepeup, index, {{1,2}});

            LHEEventProduct* product = new LHEEventProduct(hepeup,1.0);
            product->addWeight(gen::WeightsInfo("ctau",ctauLLP));
            product->addWeight(gen::WeightsInfo("llpmass",massLLP));

            return product;
        }

        void produce(edm::Event& event) override
        {
            std::unique_ptr<LHEEventProduct> product(nullptr);
            do
            {
                try
                {
                    product.reset(generate());
                }
                catch (...)
                {
                    product.reset(nullptr);
                }
            }
            while (not product);
            event.put(std::move(product));
        }

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions)
        {
            edm::ParameterSetDescription desc;
            desc.setAllowAnything();
            descriptions.addDefault(desc);
        }
};

DEFINE_FWK_INPUT_SOURCE(HNLGun);
