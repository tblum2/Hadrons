/*
 * A2AVectors.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Fionn O hOgain <fionn.o.hogain@ed.ac.uk>
 * Author: Fionn Ó hÓgáin <fionnoh@gmail.com>
 * Author: fionnoh <fionnoh@gmail.com>
 *
 * Hadrons is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * Hadrons is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Hadrons.  If not, see <http://www.gnu.org/licenses/>.
 *
 * See the full license in the file "LICENSE" in the top level distribution 
 * directory.
 */

/*  END LEGAL */
#ifndef Hadrons_MSolver_A2AVectors_hpp_
#define Hadrons_MSolver_A2AVectors_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Solver.hpp>
#include <Hadrons/EigenPack.hpp>
#include <Hadrons/A2AVectors.hpp>
#include <Hadrons/DilutedNoise.hpp>

#ifdef USE_QLATTICE
#include <qlat/grid.h>
#endif 

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                       Create all-to-all V & W vectors                      *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSolver)

class A2AVectorsPar: Serializable
{
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(A2AVectorsPar,
                                  std::string, noise,
                                  std::string, action,
                                  std::string, eigenPack,
                                  std::string, solver,
                                  std::string, output,
                                  double, mass,
                                  bool,        multiFile);
};

template <typename FImpl, typename Pack>
class TA2AVectors : public Module<A2AVectorsPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SOLVER_TYPE_ALIASES(FImpl,);
    typedef HADRONS_DEFAULT_SCHUR_A2A<FImpl> A2A;
public:
    // constructor
    TA2AVectors(const std::string name);
    // destructor
    virtual ~TA2AVectors(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    std::string  solverName_;
    unsigned int Nl_{0};
};

MODULE_REGISTER_TMP(A2AVectors, 
    ARG(TA2AVectors<FIMPL, BaseFermionEigenPack<FIMPL>>), MSolver);
MODULE_REGISTER_TMP(ZA2AVectors, 
    ARG(TA2AVectors<ZFIMPL, BaseFermionEigenPack<ZFIMPL>>), MSolver);

/******************************************************************************
 *                       TA2AVectors implementation                           *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
TA2AVectors<FImpl, Pack>::TA2AVectors(const std::string name)
: Module<A2AVectorsPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
std::vector<std::string> TA2AVectors<FImpl, Pack>::getInput(void)
{
    std::string              sub_string;
    std::vector<std::string> in;

    if (!par().eigenPack.empty())
    {
        in.push_back(par().eigenPack);
        sub_string = (!par().eigenPack.empty()) ? "_subtract" : "";
    }
    in.push_back(par().solver + sub_string);
    in.push_back(par().noise);

    return in;
}

template <typename FImpl, typename Pack>
std::vector<std::string> TA2AVectors<FImpl, Pack>::getOutput(void)
{
    std::vector<std::string> out = {getName() + "_v", getName() + "_w"};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
void TA2AVectors<FImpl, Pack>::setup(void)
{
    bool        hasLowModes = (!par().eigenPack.empty());
    std::string sub_string  = (hasLowModes) ? "_subtract" : "";
    auto        &noise      = envGet(SpinColorDiagonalNoise<FImpl>, par().noise);
    auto        &action     = envGet(FMat, par().action);
    auto        &solver     = envGet(Solver, par().solver + sub_string);
    int         Ls          = env().getObjectLs(par().action);

    if (hasLowModes)
    {
        auto &epack = envGet(Pack, par().eigenPack);
        Nl_ = epack.evec.size();
    }
    envCreate(std::vector<FermionField>, getName() + "_v", 1, 
              Nl_ + noise.fermSize(), envGetGrid(FermionField));
    envCreate(std::vector<FermionField>, getName() + "_w", 1, 
              Nl_ + noise.fermSize(), envGetGrid(FermionField));
    if (Ls > 1)
    {
        envTmpLat(FermionField, "f5", Ls);
    }
    envTmp(A2A, "a2a", 1, action, solver);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
void TA2AVectors<FImpl, Pack>::execute(void)
{
    std::string sub_string = (Nl_ > 0) ? "_subtract" : "";
    auto        &action    = envGet(FMat, par().action);
    auto        &solver    = envGet(Solver, par().solver + sub_string);
    auto        &noise     = envGet(SpinColorDiagonalNoise<FImpl>, par().noise);
    auto        &v         = envGet(std::vector<FermionField>, getName() + "_v");
    auto        &w         = envGet(std::vector<FermionField>, getName() + "_w");
    int         Ls         = env().getObjectLs(par().action);

    envGetTmp(A2A, a2a);

    if (Nl_ > 0)
    {
        LOG(Message) << "Computing all-to-all vectors "
                     << " using eigenpack '" << par().eigenPack << "' ("
                     << Nl_ << " low modes) and noise '"
                     << par().noise << "' (" << noise.fermSize() 
                     << " noise vectors)" << std::endl;
    }
    else
    {
        LOG(Message) << "Computing all-to-all vectors "
                     << " using noise '" << par().noise << "' (" << noise.fermSize() 
                     << " noise vectors)" << std::endl;
    }
    // Low modes
    for (unsigned int il = 0; il < Nl_; il++)
    {
        auto &epack  = envGet(Pack, par().eigenPack);

        startTimer("V low mode");
        LOG(Message) << "V vector i = " << il << " (low mode)" << std::endl;
        if (Ls == 1)
        {
            a2a.makeLowModeV(v[il], epack.evec[il], epack.eval[il]);
        }
        else
        {
            envGetTmp(FermionField, f5);
            a2a.makeLowModeV5D(v[il], f5, epack.evec[il], epack.eval[il]);
        }
        stopTimer("V low mode");
        startTimer("W low mode");
        LOG(Message) << "W vector i = " << il << " (low mode)" << std::endl;
        if (Ls == 1)
        {
            a2a.makeLowModeW(w[il], epack.evec[il], epack.eval[il]);
        }
        else
        {
            envGetTmp(FermionField, f5);
            a2a.makeLowModeW5D(w[il], f5, epack.evec[il], epack.eval[il]);
        }
        stopTimer("W low mode");
    }

    // High modes
    for (unsigned int ih = 0; ih < noise.fermSize(); ih++)
    {
        startTimer("V high mode");
        LOG(Message) << "V vector i = " << Nl_ + ih
                     << " (" << ((Nl_ > 0) ? "high " : "") 
                     << "stochastic mode)" << std::endl;
        if (Ls == 1)
        {
            a2a.makeHighModeV(v[Nl_ + ih], noise.getFerm(ih));
        }
        else
        {
            envGetTmp(FermionField, f5);
            a2a.makeHighModeV5D(v[Nl_ + ih], f5, noise.getFerm(ih));
        }
        stopTimer("V high mode");
        startTimer("W high mode");
        LOG(Message) << "W vector i = " << Nl_ + ih
                     << " (" << ((Nl_ > 0) ? "high " : "") 
                     << "stochastic mode)" << std::endl;
        if (Ls == 1)
        {
            a2a.makeHighModeW(w[Nl_ + ih], noise.getFerm(ih));
        }
        else
        {
            envGetTmp(FermionField, f5);
            a2a.makeHighModeW5D(w[Nl_ + ih], f5, noise.getFerm(ih));
        }
        stopTimer("W high mode");
    }

    // I/O if necessary
    if (!par().output.empty())
    {
        startTimer("V I/O");
        A2AVectorsIo::write(par().output + "_v", v, par().multiFile, vm().getTrajectory());
        stopTimer("V I/O");
        startTimer("W I/O");
        A2AVectorsIo::write(par().output + "_w", w, par().multiFile, vm().getTrajectory());
        stopTimer("W I/O");
    }
}

// Staggered A2A vectors
//
class StagSparseA2AVectorsPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(StagSparseA2AVectorsPar,
                                    std::string, action,
                                    std::string, gauge,
                                    std::string, eigenPack,
                                    std::string, output,
                                    int, inc,
                                    int, tinc,
                                    double, mass,
                                    bool, multiFile);
};

//
// Sparsened A2A vectors for staggered conserved current
//
template <typename FImpl, typename Pack>
class TStagSparseA2AVectors : public Module<StagSparseA2AVectorsPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef A2AVectorsLowStaggered<FImpl> A2A;
    typedef typename Grid::NaiveStaggeredFermionD::FermionField SparseFermionField;
    
public:
    // constructor
    TStagSparseA2AVectors(const std::string name);
    // destructor
    virtual ~TStagSparseA2AVectors(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    unsigned int Nl_{0};
};

MODULE_REGISTER_TMP(StagSparseA2AVectors,
                    ARG(TStagSparseA2AVectors<STAGIMPL, BaseFermionEigenPack<STAGIMPL>>),MSolver);

/******************************************************************************
 *                       TStagSparseA2AVectors implementation                           *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
TStagSparseA2AVectors<FImpl, Pack>::TStagSparseA2AVectors(const std::string name)
: Module<StagSparseA2AVectorsPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
std::vector<std::string> TStagSparseA2AVectors<FImpl, Pack>::getInput(void)
{
    std::vector<std::string> in;
    
    in.push_back(par().eigenPack);
    in.push_back(par().gauge);
    in.push_back(par().action);
    
    return in;
}

template <typename FImpl, typename Pack>
std::vector<std::string> TStagSparseA2AVectors<FImpl, Pack>::getOutput(void)
{
    std::vector<std::string> out = {getName() +"_v", getName() +"_w0", getName() +"_w1", getName() +"_w2"};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
void TStagSparseA2AVectors<FImpl, Pack>::setup(void)
{
    auto        &action     = envGet(FMat, par().action);
    
    auto &epack = envGet(Pack, par().eigenPack);
    Nl_ = epack.evec.size();
    envTmp(A2A, "a2a", 1, action);
    
    // Sparse Grid
    std::vector<int> blocksize(4);
    blocksize[0] = par().inc;
    blocksize[1] = par().inc;
    blocksize[2] = par().inc;
    blocksize[3] = par().tinc;
    envCreate(std::vector<SparseFermionField>, getName() + "_v", 1,
              2*Nl_, envGetCoarseGrid(SparseFermionField,blocksize));
    envCreate(std::vector<SparseFermionField>, getName() + "_w0", 1,
              2*Nl_, envGetCoarseGrid(SparseFermionField,blocksize));
    envCreate(std::vector<SparseFermionField>, getName() + "_w1", 1,
              2*Nl_, envGetCoarseGrid(SparseFermionField,blocksize));
    envCreate(std::vector<SparseFermionField>, getName() + "_w2", 1,
              2*Nl_, envGetCoarseGrid(SparseFermionField,blocksize));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
void TStagSparseA2AVectors<FImpl, Pack>::execute(void)
{
    auto    &action = envGet(FMat, par().action);
    auto    &U      = envGet(LatticeGaugeField, par().gauge);
    auto    &epack  = envGet(Pack, par().eigenPack);
    double  mass    = par().mass;
    uint64_t     nt      = env().getDim(Tp);
    uint64_t     ns      = env().getDim(Xp);
    uint64_t     glbsize = ns*ns*ns * nt;
    envGetTmp(A2A, a2a);
    
    int orthogdim=env().getNd()-1; // time dir
        
    auto &v = envGet(std::vector<SparseFermionField>, getName() + "_v");
    auto &w0 = envGet(std::vector<SparseFermionField>, getName() + "_w0");
    auto &w1 = envGet(std::vector<SparseFermionField>, getName() + "_w1");
    auto &w2 = envGet(std::vector<SparseFermionField>, getName() + "_w2");
    ScidacWriter binWriter(U.Grid()->IsBoss());
    A2AVectorsIo::Record record;
    const int traj=vm().getTrajectory();
    
    LOG(Message) << "Computing all-to-all vectors using eigenpack " << par().eigenPack << " with " << 2*Nl_ << " low modes " << std::endl;

    // Staggered Phases. Do spatial gamma only
    Lattice<iScalar<vInteger> > x(U.Grid()); LatticeCoordinate(x,0);
    Lattice<iScalar<vInteger> > y(U.Grid()); LatticeCoordinate(y,1);
    Lattice<iScalar<vInteger> > lin_z(U.Grid()); lin_z=x+y;
        
    ComplexField phases(U.Grid());
    int shift;
    FermionField temp(U.Grid());
    FermionField temp2(U.Grid());
    //SparseFermionField temp3(&sparseGrid);
    
    // step size for hypercube loop
    int step=2*par().inc;
    //std::random_device rd;  // a seed source for the random number engine
    //std::mt19937 gen(rd()); // mersenne_twister_engine seeded with rd()
    // sparsen on time slice starting anywhere between 0 and ns-1, inclusive
    // seed with rngSerial (below)
    std::uniform_int_distribution<uint32_t> uid(0, step-1);
    std::vector<uint32_t> xshift(nt);
    std::vector<uint32_t> yshift(nt);
    std::vector<uint32_t> zshift(nt);
    if(par().inc != 1){
        for(int t=0;t<nt;t++){
            xshift[t]=uid(rngSerial()._generators[0]);
            yshift[t]=uid(rngSerial()._generators[0]);
            zshift[t]=uid(rngSerial()._generators[0]);
        }
    }else{//will miss site 0000 otherwise
        for(int t=0;t<nt;t++){
            xshift[t]=0;
            yshift[t]=0;
            zshift[t]=0;
        }
    }
    CartesianCommunicator::BroadcastWorld(0,(void *)&xshift[0],sizeof(uint32_t)*xshift.size());
    CartesianCommunicator::BroadcastWorld(0,(void *)&yshift[0],sizeof(uint32_t)*yshift.size());
    CartesianCommunicator::BroadcastWorld(0,(void *)&zshift[0],sizeof(uint32_t)*zshift.size());
    
    //save for later
    std::vector<complex<double>> evalM(2*Nl_);
    
    // global to local coord shifts
    int locx=U.Grid()->_ldimensions[0];
    int locy=U.Grid()->_ldimensions[1];
    int locz=U.Grid()->_ldimensions[2];
    int loct=U.Grid()->_ldimensions[3];
    int lstartx=U.Grid()->_lstart[0];
    int lstarty=U.Grid()->_lstart[1];
    int lstartz=U.Grid()->_lstart[2];
    int lstartt=U.Grid()->_lstart[3];
    printf("my node %d local tsites= %d\n", U.Grid()->_processor_coor[3],loct);
    printf("my node %d global t start= %d\n", U.Grid()->_processor_coor[3], lstartt);
    LOG(Message) << "xshift" << xshift << std::endl;
    LOG(Message) << "yshift" << yshift << std::endl;
    LOG(Message) << "zshift" << zshift << std::endl;

    for (unsigned int il = 0; il < 2*Nl_; il++)
    {
        // eval of unpreconditioned Dirac op
        std::complex<double> eval(mass,sqrt(epack.eval[il/2]-mass*mass));
        
        startTimer("W low mode");
        LOG(Message) << "W vector i = " << il << " (low modes)" << std::endl;
        // don't divide by lambda. Do it in contraction since it is complex
        a2a.makeLowModeW(temp, epack.evec[il/2], eval, il%2);
        //LOG(Message) << "evec "<< il <<" norm^2 " << norm2(epack.evec[il/2]) << std::endl;
        stopTimer("W low mode");
        
        // debug //////////////////////////////////
//        ColourVector cvec;
//        int lxrb=epack.evec[il/2].Grid()->_ldimensions[0];
//        int lyrb=epack.evec[il/2].Grid()->_ldimensions[1];
//        int lzrb=epack.evec[il/2].Grid()->_ldimensions[2];
//        int ltrb=epack.evec[il/2].Grid()->_ldimensions[3];
//        for(int tl=0;tl<ltrb;tl++){
//            for(int zl=0;zl<lzrb;zl++){
//                for(int yl=0;yl<lyrb;yl++){
//                    for(int xl=0;xl<lxrb;xl++){
//                        Coordinate site(Nd);
//                        site[0]=xl;site[1]=yl;site[2]=zl;site[3]=tl;
//                        peekLocalSite(cvec,epack.evec[il/2],site);
//                        LOG(Message) << "evec " << il << " " << site << " " << cvec << std::endl;
//                    }
//                }
//            }
//        }
//        for(int tl=0;tl<loct;tl++){
//            for(int zl=0;zl<locz;zl++){
//                for(int yl=0;yl<locy;yl++){
//                    for(int xl=0;xl<locx;xl++){
//                        Coordinate site(Nd);
//                        site[0]=xl;site[1]=yl;site[2]=zl;site[3]=tl;
//                        peekLocalSite(cvec,temp,site);
//                        LOG(Message) << "W vec " << il << " " << site << " " << cvec << std::endl;
//                    }
//                }
//            }
//        }
        // end debug /////////////////////////////
        
        il%2 ? eval=conjugate(eval) : eval ;
        evalM[il]=eval;
        
        for (int mu=0;mu<3;mu++){
            
            phases=1.0;
            if(mu==1){
                phases = where( mod(x    ,2)==(Integer)0, phases,-phases);
            } else if(mu==2){
                phases = where( mod(lin_z,2)==(Integer)0, phases,-phases);
            }
            LatticeColourMatrix Umu(U.Grid());
            Umu = PeekIndex<LorentzIndex>(U,mu);
            Umu *= phases;
        
            // v vec is shifted and * link for conserved current
            temp2 = Umu*Cshift(temp, mu, 1);
                        
            thread_for(t,loct,{
            //for(int t=0; t<loct;t++){//debug
                int tglb=t+lstartt;
                // same random shift for t, t+1 in same hypercube
                if(t%2 == 1) continue;
                
                Coordinate site(Nd);
                Coordinate sparseSite(Nd);
                ColourVector vec;
                
                
                // loop over hypercubes on time slice
                // and sparsen, keep all sites in hypercube
                for(int z=0;z<ns;z+=step){
                    int zg=(zshift[tglb]+z)%ns;
                    for(int zl=0;zl<locz;zl++){
                        int zgp=zl+lstartz;
                        if(zgp==zg || zgp==(zg+1)%ns){
                            site[2]=zl;
                            if(par().inc==1){
                                sparseSite[2]=site[2];
                            }else if(zshift[tglb]!=step-1){
                                sparseSite[2]=2*int(site[2]/step) + (site[2]+zshift[tglb])%2;
                            }else{
                                sparseSite[2]=2*int(site[2]/step) + (site[2])%2;
                            }
                            for(int y=0;y<ns;y+=step){
                                int yg=(yshift[tglb]+y)%ns;
                                for(int yl=0;yl<locy;yl++){
                                    int ygp=yl+lstarty;
                                    if(ygp==yg || ygp==(yg+1)%ns){
                                        site[1]=yl;
                                        if(par().inc==1){
                                            sparseSite[1]=site[1];
                                        }else if(yshift[tglb]!=step-1){
                                            sparseSite[1]=2*int(site[1]/step) + (site[1]+yshift[tglb])%2;
                                        }else{
                                            sparseSite[1]=2*int(site[1]/step) + (site[1])%2;
                                        }
                                        for(int x=0;x<ns;x+=step){
                                            int xg=(xshift[tglb]+x)%ns;
                                            for(int xl=0;xl<locx;xl++){
                                                int xgp=xl+lstartx;
                                                if(xgp==xg || xgp==(xg+1)%ns){
                                                    site[0]=xl;
                                                    if(par().inc==1){
                                                        sparseSite[0]=site[0];
                                                    }else if(xshift[tglb]!=step-1){
                                                        sparseSite[0]=2*int(site[0]/step) + (site[0]+xshift[tglb])%2;
                                                    }else{
                                                        sparseSite[0]=2*int(site[0]/step) + (site[0])%2;
                                                    }
                                                    for(int that=0;that<2;that++){
                                                        
                                                        site[3]=t+that;
                                                        sparseSite[3]=site[3];
                                                        
                                                        if(mu==0){// do v once
                                                            peekLocalSite(vec,temp,site);
//LOG(Message) << "Sparse V vec " << il << " " << sparseSite << " " << vec << std::endl;
                                                            pokeLocalSite(vec,v[il],sparseSite);
                                                            peekLocalSite(vec,temp2,site);
                                                            pokeLocalSite(vec,w0[il],sparseSite);
//LOG(Message) << "Sparse W0 vec " << il << " " << sparseSite << " " << vec << std::endl;
                                                        }else if(mu==1){
                                                            peekLocalSite(vec,temp2,site);
                                                            pokeLocalSite(vec,w1[il],sparseSite);
                                                        }else if(mu==2){
                                                            peekLocalSite(vec,temp2,site);
                                                            pokeLocalSite(vec,w2[il],sparseSite);
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            });
        }// end mu
    }// end evecs
    
    std::string dir = dirname(par().output);
    if (!par().output.empty())
    {
        int status = mkdir(dir);
        if (status)
        {
            HADRONS_ERROR(Io, "cannot create directory '" + dir
                          + "' ( " + std::strerror(errno) + ")");
        }
        startTimer("V I/O");
        A2AVectorsIo::write(par().output + "_v", v, par().multiFile, vm().getTrajectory());
        stopTimer("V I/O");
        startTimer("W I/O");
        A2AVectorsIo::write(par().output + "_w0", w0, par().multiFile, vm().getTrajectory());
        A2AVectorsIo::write(par().output + "_w1", w1, par().multiFile, vm().getTrajectory());
        A2AVectorsIo::write(par().output + "_w2", w2, par().multiFile, vm().getTrajectory());
        stopTimer("W I/O");
    }
    
    if ( env().getGrid()->IsBoss() ) {
        
        std::string eval_filename;
        if (!par().output.empty()){
            eval_filename=A2AVectorsIo::evalFilename(par().output,vm().getTrajectory());
        }else{
            eval_filename=A2AVectorsIo::evalFilename("evals",vm().getTrajectory());
        }
        A2AVectorsIo::initEvalFile(eval_filename,
                                   evalM.size());// total size
        A2AVectorsIo::saveEvalBlock(eval_filename,
                                    evalM.data(),
                                    0,// start of chunk
                                    2*Nl_);// size of chunk saved
    }
}


//
// Sparsened A2A vectors for staggered conserved current, streaming evecs from Grid LIME files
//

class StagSparseA2AVectorsGridIoPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(StagSparseA2AVectorsGridIoPar,
                                    std::string, action,
                                    std::string, gauge,
                                    std::string, evecPath,
                                    std::string, output,
                                    int, numEvecs,
                                    int, inc,
                                    int, tinc,
                                    double, mass,
                                    bool, multiFile);
};


template <typename FImpl>
class TStagSparseA2AVectorsGridIo : public Module<StagSparseA2AVectorsGridIoPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef A2AVectorsLowStaggered<FImpl> A2A;
    typedef typename Grid::NaiveStaggeredFermionD::FermionField SparseFermionField;

public:
    // constructor
    TStagSparseA2AVectorsGridIo(const std::string name);
    // destructor
    virtual ~TStagSparseA2AVectorsGridIo(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    unsigned int Nl_{0};
};

MODULE_REGISTER_TMP(StagSparseA2AVectorsGridIo,
                    ARG(TStagSparseA2AVectorsGridIo<STAGIMPL>), MSolver);

/******************************************************************************
 *                       TStagSparseA2AVectorsGridIo implementation           *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TStagSparseA2AVectorsGridIo<FImpl>::TStagSparseA2AVectorsGridIo(const std::string name)
: Module<StagSparseA2AVectorsGridIoPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TStagSparseA2AVectorsGridIo<FImpl>::getInput(void)
{
    std::vector<std::string> in;

    in.push_back(par().gauge);
    in.push_back(par().action);

    return in;
}

template <typename FImpl>
std::vector<std::string> TStagSparseA2AVectorsGridIo<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName() +"_v", getName() +"_w0", getName() +"_w1", getName() +"_w2"};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TStagSparseA2AVectorsGridIo<FImpl>::setup(void)
{
    auto &action = envGet(FMat, par().action);

    Nl_ = par().numEvecs;
    envTmp(A2A, "a2a", 1, action);

    // Sparse Grid
    std::vector<int> blocksize(4);
    blocksize[0] = par().inc;
    blocksize[1] = par().inc;
    blocksize[2] = par().inc;
    blocksize[3] = par().tinc;
    envCreate(std::vector<SparseFermionField>, getName() + "_v", 1,
              2*Nl_, envGetCoarseGrid(SparseFermionField,blocksize));
    envCreate(std::vector<SparseFermionField>, getName() + "_w0", 1,
              2*Nl_, envGetCoarseGrid(SparseFermionField,blocksize));
    envCreate(std::vector<SparseFermionField>, getName() + "_w1", 1,
              2*Nl_, envGetCoarseGrid(SparseFermionField,blocksize));
    envCreate(std::vector<SparseFermionField>, getName() + "_w2", 1,
              2*Nl_, envGetCoarseGrid(SparseFermionField,blocksize));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TStagSparseA2AVectorsGridIo<FImpl>::execute(void)
{
    auto     &action = envGet(FMat, par().action);
    auto     &U      = envGet(LatticeGaugeField, par().gauge);
    double    mass   = par().mass;
    uint64_t  nt     = env().getDim(Tp);
    uint64_t  ns     = env().getDim(Xp);
    envGetTmp(A2A, a2a);

    auto &v  = envGet(std::vector<SparseFermionField>, getName() + "_v");
    auto &w0 = envGet(std::vector<SparseFermionField>, getName() + "_w0");
    auto &w1 = envGet(std::vector<SparseFermionField>, getName() + "_w1");
    auto &w2 = envGet(std::vector<SparseFermionField>, getName() + "_w2");
    const int traj = vm().getTrajectory();

    // scratch space: only one eigenvector in memory at a time
    FermionField tempEvec(env().getRbGrid());
    tempEvec.Checkerboard() = Odd;
    assert(tempEvec.Checkerboard() == Odd);
    LOG(Message) << " Checkerboard grid: " << std::endl;
    tempEvec.Grid()->show_decomposition();
    LOG(Message) << " Checkerboard dimension: "<< tempEvec.Grid()->_checker_dim  << std::endl;
    RealD currentEval = 0.;
    PackRecord packRecord;
    ScidacReader binReader;

    // build the base filename from the stem + trajectory (Grid EigenPack convention)
    std::string t       = "." + std::to_string(traj);
    std::string stem    = par().evecPath + t;          // directory for multiFile
    std::string binFile = par().evecPath + t + ".bin"; // single-file path

   

    LOG(Message) << "Computing sparse A2A vectors streaming " << 2*Nl_
                 << " low modes from " << (par().multiFile ? stem : binFile) << std::endl;

    // Staggered Phases. Do spatial gamma only
    Lattice<iScalar<vInteger> > x(U.Grid()); LatticeCoordinate(x,0);
    Lattice<iScalar<vInteger> > y(U.Grid()); LatticeCoordinate(y,1);
    Lattice<iScalar<vInteger> > lin_z(U.Grid()); lin_z=x+y;

    ComplexField phases(U.Grid());
    FermionField temp(U.Grid());
    FermionField temp2(U.Grid());

    int step = 2*par().inc;
    std::uniform_int_distribution<uint32_t> uid(0, step-1);
    std::vector<uint32_t> xshift(nt);
    std::vector<uint32_t> yshift(nt);
    std::vector<uint32_t> zshift(nt);
    if (par().inc != 1)
    {
        for (int tt=0; tt<(int)nt; tt++)
        {
            xshift[tt]=uid(rngSerial()._generators[0]);
            yshift[tt]=uid(rngSerial()._generators[0]);
            zshift[tt]=uid(rngSerial()._generators[0]);
        }
    }
    else
    {
        for (int tt=0; tt<(int)nt; tt++)
        { xshift[tt]=0; yshift[tt]=0; zshift[tt]=0; }
    }
    CartesianCommunicator::BroadcastWorld(0,(void *)&xshift[0],sizeof(uint32_t)*xshift.size());
    CartesianCommunicator::BroadcastWorld(0,(void *)&yshift[0],sizeof(uint32_t)*yshift.size());
    CartesianCommunicator::BroadcastWorld(0,(void *)&zshift[0],sizeof(uint32_t)*zshift.size());

    std::vector<complex<double>> evalM(2*Nl_);

    int locx    = U.Grid()->_ldimensions[0];
    int locy    = U.Grid()->_ldimensions[1];
    int locz    = U.Grid()->_ldimensions[2];
    int loct    = U.Grid()->_ldimensions[3];
    int lstartx = U.Grid()->_lstart[0];
    int lstarty = U.Grid()->_lstart[1];
    int lstartz = U.Grid()->_lstart[2];
    int lstartt = U.Grid()->_lstart[3];
    LOG(Message) << "xshift" << xshift << std::endl;
    LOG(Message) << "yshift" << yshift << std::endl;
    LOG(Message) << "zshift" << zshift << std::endl;

    // For single-file mode: open once before the loop to avoid repeated
    // MPI_File_open/close cycles which exhaust GPFS/OMPIO resources (~974 opens).
    if (!par().multiFile)
    {
        binReader.open(binFile);
        EigenPackIo::readHeader(packRecord, binReader);
    }

    for (unsigned int il = 0; il < 2*Nl_; il++)
    {
        // read a new eigenvector from disk every other iteration
        if (il % 2 == 0)
        {
            int k = il / 2;
            startTimer("evec read");
            if (par().multiFile)
            {
                std::string fname = stem + "/v" + std::to_string(k) + ".bin";
                binReader.open(fname);
                EigenPackIo::readHeader(packRecord, binReader);
                EigenPackIo::readElement(tempEvec, currentEval, k, binReader);
                binReader.close();
            }
            else
            {
                EigenPackIo::readElement(tempEvec, currentEval, k, binReader);
            }
            stopTimer("evec read");
	    if (il == 0)
            {
                LOG(Message) << "tempEvec grid dimensions: " << tempEvec.Grid()->GlobalDimensions() << std::endl;
                LOG(Message) << "Full grid dimensions:     " << U.Grid()->GlobalDimensions() << std::endl;
                LOG(Message) << "RbGrid dimensions:        " << env().getRbGrid()->GlobalDimensions() << std::endl;
                LOG(Message) << "norm2(tempEvec)=          " << norm2(tempEvec) << std::endl;
	    	}	
        }

        std::complex<double> eval(mass, sqrt(currentEval));
        
        // double lambda;
        // if (currentEval < mass * mass)
        // {
        //     lambda = sqrt(currentEval);
        //     if (il == 0)
        //         LOG(Message) << "Eigenpack convention: massless DdagD (currentEval < m^2)" << std::endl;
        // }
        // else
        // {
        //     lambda = sqrt(currentEval - mass * mass);
        //     if (il == 0)
        //         LOG(Message) << "Eigenpack convention: massive (D+m)dag(D+m) (currentEval >= m^2)" << std::endl;
        // }
        // std::complex<double> eval(mass, lambda);

	startTimer("W low mode");
        LOG(Message) << "W vector i = " << il << " (low modes)" << std::endl;
        // don't divide by lambda — do it in contraction since it is complex
        a2a.makeLowModeW(temp, tempEvec, eval, il%2);
        stopTimer("W low mode");
        
        
	il%2 ? eval=conjugate(eval) : eval ;
        evalM[il]=eval;

        v[il]  = Zero();
        w0[il] = Zero();
        w1[il] = Zero();
        w2[il] = Zero();

        for (int mu=0; mu<3; mu++)
        {
            phases=1.0;
            if (mu==1)
            {
                phases = where( mod(x    ,2)==(Integer)0, phases,-phases);
            }
            else if (mu==2)
            {
                phases = where( mod(lin_z,2)==(Integer)0, phases,-phases);
            }
            LatticeColourMatrix Umu(U.Grid());
            Umu = PeekIndex<LorentzIndex>(U,mu);
            Umu *= phases;

            // v vec is shifted and * link for conserved current
            temp2 = Umu*Cshift(temp, mu, 1);

            thread_for(tt,loct,{
                int tglb=tt+lstartt;
                // same random shift for t, t+1 in same hypercube
                if (tt%2 == 1) continue;

                Coordinate site(Nd);
                Coordinate sparseSite(Nd);
                ColourVector vec;

                // loop over hypercubes on time slice and sparsen
                for (int z=0; z<(int)ns; z+=step) {
                    int zg=(zshift[tglb]+z)%ns;
                    for (int zl=0; zl<locz; zl++) {
                        int zgp=zl+lstartz;
                        if (zgp==zg || zgp==(zg+1)%ns) {
                            site[2]=zl;
                            if (par().inc==1) {
                                sparseSite[2]=site[2];
                            } else if (zshift[tglb]!=step-1) {
                                sparseSite[2]=2*int(site[2]/step) + (site[2]+zshift[tglb])%2;
                            } else {
                                sparseSite[2]=2*int(site[2]/step) + (site[2])%2;
                            }
                            for (int y=0; y<(int)ns; y+=step) {
                                int yg=(yshift[tglb]+y)%ns;
                                for (int yl=0; yl<locy; yl++) {
                                    int ygp=yl+lstarty;
                                    if (ygp==yg || ygp==(yg+1)%ns) {
                                        site[1]=yl;
                                        if (par().inc==1) {
                                            sparseSite[1]=site[1];
                                        } else if (yshift[tglb]!=step-1) {
                                            sparseSite[1]=2*int(site[1]/step) + (site[1]+yshift[tglb])%2;
                                        } else {
                                            sparseSite[1]=2*int(site[1]/step) + (site[1])%2;
                                        }
                                        for (int xx=0; xx<(int)ns; xx+=step) {
                                            int xg=(xshift[tglb]+xx)%ns;
                                            for (int xl=0; xl<locx; xl++) {
                                                int xgp=xl+lstartx;
                                                if (xgp==xg || xgp==(xg+1)%ns) {
                                                    site[0]=xl;
                                                    if (par().inc==1) {
                                                        sparseSite[0]=site[0];
                                                    } else if (xshift[tglb]!=step-1) {
                                                        sparseSite[0]=2*int(site[0]/step) + (site[0]+xshift[tglb])%2;
                                                    } else {
                                                        sparseSite[0]=2*int(site[0]/step) + (site[0])%2;
                                                    }
                                                    for (int that=0; that<2; that++) {
                                                        site[3]=tt+that;
                                                        sparseSite[3]=site[3];
                                                        if (mu==0) {
                                                            peekLocalSite(vec,temp,site);
                                                            pokeLocalSite(vec,v[il],sparseSite);
                                                            peekLocalSite(vec,temp2,site);
                                                            pokeLocalSite(vec,w0[il],sparseSite);
                                                        } else if (mu==1) {
                                                            peekLocalSite(vec,temp2,site);
                                                            pokeLocalSite(vec,w1[il],sparseSite);
                                                        } else if (mu==2) {
                                                            peekLocalSite(vec,temp2,site);
                                                            pokeLocalSite(vec,w2[il],sparseSite);
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            });
        }// end mu
    }// end evecs

    if (!par().multiFile)
        binReader.close();

    std::string dir = dirname(par().output);
    if (!par().output.empty())
    {
        int status = mkdir(dir);
        if (status)
        {
            HADRONS_ERROR(Io, "cannot create directory '" + dir
                          + "' ( " + std::strerror(errno) + ")");
        }
        startTimer("V I/O");
        A2AVectorsIo::write(par().output + "_v",  v,  par().multiFile, vm().getTrajectory());
        stopTimer("V I/O");
        startTimer("W I/O");
        A2AVectorsIo::write(par().output + "_w0", w0, par().multiFile, vm().getTrajectory());
        A2AVectorsIo::write(par().output + "_w1", w1, par().multiFile, vm().getTrajectory());
        A2AVectorsIo::write(par().output + "_w2", w2, par().multiFile, vm().getTrajectory());
        stopTimer("W I/O");
    }

    if (env().getGrid()->IsBoss())
    {
        std::string eval_filename;
        if (!par().output.empty())
            eval_filename = A2AVectorsIo::evalFilename(par().output, vm().getTrajectory());
        else
            eval_filename = A2AVectorsIo::evalFilename("evals", vm().getTrajectory());
        A2AVectorsIo::initEvalFile(eval_filename, evalM.size());
        A2AVectorsIo::saveEvalBlock(eval_filename, evalM.data(), 0, 2*Nl_);
    }
}


#if defined(USE_QLATTICE)
#include <qlat/qlat.h>


class StagSparseA2AVectorsEvecIoPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(StagSparseA2AVectorsEvecIoPar,
                                    std::string, action,
                                    std::string, gauge,
                                    std::string, evecPath,
                                    std::string, evecTag,
                                    std::string, output,
                                    int, numEvecs,
                                    int, inc,
                                    int, tinc,
                                    double, mass,
                                    bool, multiFile);
};

//
// Sparsened A2A vectors for staggered conserved current
//
template <typename FImpl>
class TStagSparseA2AVectorsEvecIo : public Module<StagSparseA2AVectorsEvecIoPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef A2AVectorsLowStaggered<FImpl> A2A;
    typedef typename Grid::NaiveStaggeredFermionD::FermionField SparseFermionField;
    
public:
    // constructor
    TStagSparseA2AVectorsEvecIo(const std::string name);
    // destructor
    virtual ~TStagSparseA2AVectorsEvecIo(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    unsigned int Nl_{0};
};

MODULE_REGISTER_TMP(StagSparseA2AVectorsEvecIo,
                    ARG(TStagSparseA2AVectorsEvecIo<STAGIMPL>),MSolver);

/******************************************************************************
 *                       TStagSparseA2AVectorsEvecIo implementation                           *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TStagSparseA2AVectorsEvecIo<FImpl>::TStagSparseA2AVectorsEvecIo(const std::string name)
: Module<StagSparseA2AVectorsEvecIoPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TStagSparseA2AVectorsEvecIo<FImpl>::getInput(void)
{
    std::vector<std::string> in;
    
    in.push_back(par().gauge);
    in.push_back(par().action);
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TStagSparseA2AVectorsEvecIo<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName() +"_v", getName() +"_w0", getName() +"_w1", getName() +"_w2"};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TStagSparseA2AVectorsEvecIo<FImpl>::setup(void)
{
    auto        &action     = envGet(FMat, par().action);
    
    Nl_ = par().numEvecs;
    envTmp(A2A, "a2a", 1, action);
    
    // Sparse Grid
    std::vector<int> blocksize(4);
    blocksize[0] = par().inc;
    blocksize[1] = par().inc;
    blocksize[2] = par().inc;
    blocksize[3] = par().tinc;
    envCreate(std::vector<SparseFermionField>, getName() + "_v", 1,
              2*Nl_, envGetCoarseGrid(SparseFermionField,blocksize));
    envCreate(std::vector<SparseFermionField>, getName() + "_w0", 1,
              2*Nl_, envGetCoarseGrid(SparseFermionField,blocksize));
    envCreate(std::vector<SparseFermionField>, getName() + "_w1", 1,
              2*Nl_, envGetCoarseGrid(SparseFermionField,blocksize));
    envCreate(std::vector<SparseFermionField>, getName() + "_w2", 1,
              2*Nl_, envGetCoarseGrid(SparseFermionField,blocksize));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TStagSparseA2AVectorsEvecIo<FImpl>::execute(void)
{
    auto    &action = envGet(FMat, par().action);
    auto    &U      = envGet(LatticeGaugeField, par().gauge);
    double  mass    = par().mass;
    uint64_t     nt      = env().getDim(Tp);
    uint64_t     ns      = env().getDim(Xp);
    uint64_t     glbsize = ns*ns*ns * nt;
    envGetTmp(A2A, a2a);
    
    int orthogdim=env().getNd()-1; // time dir
        
    auto &v = envGet(std::vector<SparseFermionField>, getName() + "_v");
    auto &w0 = envGet(std::vector<SparseFermionField>, getName() + "_w0");
    auto &w1 = envGet(std::vector<SparseFermionField>, getName() + "_w1");
    auto &w2 = envGet(std::vector<SparseFermionField>, getName() + "_w2");
    
    A2AVectorsIo::Record record;
    const int traj=vm().getTrajectory();
    
    LOG(Message) << "Computing sparse all-to-all vectors using evecs " << par().evecPath << " with " << 2*Nl_ << " low modes and " << par().inc << " sparse factor" << std::endl;
    LOG(Message) << " Full grid: " << std::endl;
    U.Grid()->show_decomposition();
    LOG(Message) << " Sparse grid: " << std::endl;
    v[0].Grid()->show_decomposition();
    
    // Staggered Phases. Do spatial gamma only
    Lattice<iScalar<vInteger> > x(U.Grid()); LatticeCoordinate(x,0);
    Lattice<iScalar<vInteger> > y(U.Grid()); LatticeCoordinate(y,1);
    Lattice<iScalar<vInteger> > lin_z(U.Grid()); lin_z=x+y;
        
    ComplexField phases(U.Grid());
    int shift;
    FermionField temp(U.Grid());
    FermionField temp2(U.Grid());
    FermionField tempRb(env().getRbGrid());
    tempRb.Checkerboard() = Odd;
    assert(tempRb.Checkerboard() == Odd);
    LOG(Message) << " Checkerboard grid: " << std::endl;
    tempRb.Grid()->show_decomposition();
    LOG(Message) << " Checkerboard dimension: "<< tempRb.Grid()->_checker_dim  << std::endl;

    // step size for hypercube loop
    int step=2*par().inc;
    //std::random_device rd;  // a seed source for the random number engine
    //std::mt19937 gen(rd()); // mersenne_twister_engine seeded with rd()
    // sparsen on time slice starting anywhere between 0 and ns-1, inclusive
    // seed with rngSerial (below)
    std::uniform_int_distribution<uint32_t> uid(0, step-1);
    std::vector<uint32_t> xshift(nt);
    std::vector<uint32_t> yshift(nt);
    std::vector<uint32_t> zshift(nt);
    if(par().inc != 1){
        for(int t=0;t<nt;t++){
            xshift[t]=uid(rngSerial()._generators[0]);
            yshift[t]=uid(rngSerial()._generators[0]);
            zshift[t]=uid(rngSerial()._generators[0]);
        }
    }else{//will miss site 0000 otherwise
        for(int t=0;t<nt;t++){
            xshift[t]=0;
            yshift[t]=0;
            zshift[t]=0;
        }
    }
    CartesianCommunicator::BroadcastWorld(0,(void *)&xshift[0],sizeof(uint32_t)*xshift.size());
    CartesianCommunicator::BroadcastWorld(0,(void *)&yshift[0],sizeof(uint32_t)*yshift.size());
    CartesianCommunicator::BroadcastWorld(0,(void *)&zshift[0],sizeof(uint32_t)*zshift.size());
    
    //save for later
    std::vector<complex<double>> evalM(2*Nl_);
    
    // global to local coord shifts
    int locx=U.Grid()->_ldimensions[0];
    int locy=U.Grid()->_ldimensions[1];
    int locz=U.Grid()->_ldimensions[2];
    int loct=U.Grid()->_ldimensions[3];
    int lstartx=U.Grid()->_lstart[0];
    int lstarty=U.Grid()->_lstart[1];
    int lstartz=U.Grid()->_lstart[2];
    int lstartt=U.Grid()->_lstart[3];
    
    {
        qlat::Coordinate qcoor;
        qlat::Coordinate qsizes;
        Coordinate  psizes = tempRb.Grid()->ProcessorGrid();
        Coordinate  pcoor  = tempRb.Grid()->ThisProcessorCoor();
        qsizes[0]=psizes[0];
        qsizes[1]=psizes[1];
        qsizes[2]=psizes[2];
        qsizes[3]=psizes[3];
        qcoor[0]=pcoor[0];
        qcoor[1]=pcoor[1];
        qcoor[2]=pcoor[2];
        qcoor[3]=pcoor[3];
        qlat::begin_once(qlat::index_from_coordinate(qcoor,qsizes),qsizes);
    }
    //printf("my node %d local tsites= %d\n", tempRb.Grid()->_processor_coor[3],loct);
    //printf("my node %d global t start= %d\n", U.Grid()->_processor_coor[3], lstartt);
    LOG(Message) << "xshift" << xshift << std::endl;
    LOG(Message) << "yshift" << yshift << std::endl;
    LOG(Message) << "zshift" << zshift << std::endl;
    qlat::Coordinate grid_layout;
    Coordinate  gLattice = tempRb.Grid()->GlobalDimensions();
    grid_layout[0]=gLattice[0];
    grid_layout[1]=gLattice[1];
    grid_layout[2]=gLattice[2];
    grid_layout[3]=gLattice[3];
    
    // for reading evals
    VecRecord _vecRecord;
    
    for (unsigned int il = 0; il < 2*Nl_; il++)
    {
        //startTimer("W low mode");
        LOG(Message) << "W vector i = " << il << " (low modes)" << std::endl;
        // don't divide by lambda. Do it in contraction since it is complex
        // read in evec...
        if(il%2==0){
            long evec_idx = il / 2;
            std::string tag = qlat::ssprintf("meta%ld.txt", evec_idx);
            //old tag
            std::string file = qlat::ssprintf("%s/meta%ld.txt",     par().evecTag.c_str(), evec_idx);
            qlat::Field<Complex> f;
            bool IsFileFloat = true;
            bool IsBigEndian = true;
            if (!IsFileFloat) {
                long total_bytes = qlat::read_field(f,par().evecPath,tag);
                if (total_bytes == 0) {
                    total_bytes = qlat::read_field(f,par().evecPath,file);
                }
                if (total_bytes == 0) {
                    LOG(Message) << "Failed to read i = " << il << " (low modes)" << std::endl;
                }
                if (IsBigEndian) {
                    qlat::to_from_big_endian(get_data(f));
                } else {
                    qlat::to_from_little_endian(get_data(f));
                }
            } else {
                qlat::Field<ComplexF> ff;
                long total_bytes = qlat::read_field(ff,par().evecPath,tag);
                LOG(Message) << "Total bytes read i = " << il << " " << total_bytes << std::endl;
                if (total_bytes == 0) {
                    total_bytes = qlat::read_field(ff,par().evecPath,file);
                }
                if (total_bytes == 0) {
                    LOG(Message) << "Failed to read single precision i = " << il << " (low modes)" << std::endl;
                }
                if (IsBigEndian) {
                    qlat::to_from_big_endian(get_data(ff));
                } else {
                    qlat::to_from_little_endian(get_data(ff));
                }
                const qlat::Geometry geo = ff.geo();
                f.init(geo, ff.multiplicity);
                qthread_for(index, geo.local_volume(), {
                    qlat::Vector<qlat::ComplexF> v1 = ff.get_elems(index);
                    qlat::Vector<qlat::ComplexD> v2 = f.get_elems(index);
                    for (int m = 0; m < f.multiplicity; ++m) {
                        v2[m] = v1[m];
                    }
                })
            }
            qassert(f.multiplicity == 3);
            {
                const qlat::Coordinate total_site = f.geo().total_site();
                qassert(total_site[0] == grid_layout[0]);
                qassert(total_site[1] == grid_layout[1]);
                qassert(total_site[2] == grid_layout[2]);
                qassert(total_site[3] == grid_layout[3]);
                qassert(U.Grid()->GlobalDimensions()[0] == total_site[0] * 2);
                qassert(U.Grid()->GlobalDimensions()[1] == total_site[1]);
                qassert(U.Grid()->GlobalDimensions()[2] == total_site[2]);
                qassert(U.Grid()->GlobalDimensions()[3] == total_site[3]);
                
            }
            {
                qthread_for(index, f.geo().local_volume(), {
                    qlat::Coordinate xl = f.geo().coordinate_from_index(index);
                    xl[0] *= 2;
                    Coordinate coor = qlat::grid_convert(xl);
                    qlat::Vector<ComplexD> v = f.get_elems(index);
                    qlat::array<ComplexD, 3> fs;
                    for (int m = 0; m < 3; ++m) {
                        fs[m] = v[m];
                    }
                    pokeLocalSite(fs, tempRb, coor);
                });
            }
            //LOG(Message) << "evec "<< il <<" norm^2 " << norm2(tempRb) << std::endl;
            // eval of unpreconditioned Dirac op from lime meta file
            GridLimeReader GLR;
            std::string limefile = qlat::ssprintf("%s/../meta%ld.txt",     par().evecPath.c_str(), evec_idx);
            GLR.open(limefile);
            GLR.readLimeObject(_vecRecord, _vecRecord.SerialisableClassName(),
                               std::string(SCIDAC_RECORD_XML));
            GLR.close();
            
        }
        
        std::complex<RealD> eval(mass,sqrt(_vecRecord.eval-mass*mass));
        a2a.makeLowModeW(temp, tempRb, eval, il%2);
        //stopTimer("W low mode");
        
        il%2 ? eval=conjugate(eval) : eval ;
        evalM[il]=eval;
        
        v[il] = Zero();
        w0[il] = Zero();
        w1[il] = Zero();
        w2[il] = Zero();
    
        for (int mu=0;mu<3;mu++){
            
            phases=1.0;
            if(mu==1){
                phases = where( mod(x    ,2)==(Integer)0, phases,-phases);
            } else if(mu==2){
                phases = where( mod(lin_z,2)==(Integer)0, phases,-phases);
            }
            LatticeColourMatrix Umu(U.Grid());
            Umu = PeekIndex<LorentzIndex>(U,mu);
            Umu *= phases;
        
            // v vec is shifted and * link for conserved current
            temp2 = Umu*Cshift(temp, mu, 1);
                        
            thread_for(t,loct,{
            //for(int t=0; t<loct;t++){//debug
                int tglb=t+lstartt;
                // same random shift for t, t+1 in same hypercube
                if(t%2 == 1) continue;
                
                Coordinate site(Nd);
                Coordinate sparseSite(Nd);
                ColourVector vec;
                
                
                // loop over hypercubes on time slice
                // and sparsen, keep all sites in hypercube
                for(int z=0;z<ns;z+=step){
                    int zg=(zshift[tglb]+z)%ns;
                    for(int zl=0;zl<locz;zl++){
                        int zgp=zl+lstartz;
                        if(zgp==zg || zgp==(zg+1)%ns){
                            site[2]=zl;
                            if(par().inc==1){
                                sparseSite[2]=site[2];
                            }else if(zshift[tglb]!=step-1){
                                sparseSite[2]=2*int(site[2]/step) + (site[2]+zshift[tglb])%2;
                            }else{
                                sparseSite[2]=2*int(site[2]/step) + (site[2])%2;
                            }
                            for(int y=0;y<ns;y+=step){
                                int yg=(yshift[tglb]+y)%ns;
                                for(int yl=0;yl<locy;yl++){
                                    int ygp=yl+lstarty;
                                    if(ygp==yg || ygp==(yg+1)%ns){
                                        site[1]=yl;
                                        if(par().inc==1){
                                            sparseSite[1]=site[1];
                                        }else if(yshift[tglb]!=step-1){
                                            sparseSite[1]=2*int(site[1]/step) + (site[1]+yshift[tglb])%2;
                                        }else{
                                            sparseSite[1]=2*int(site[1]/step) + (site[1])%2;
                                        }
                                        for(int x=0;x<ns;x+=step){
                                            int xg=(xshift[tglb]+x)%ns;
                                            for(int xl=0;xl<locx;xl++){
                                                int xgp=xl+lstartx;
                                                if(xgp==xg || xgp==(xg+1)%ns){
                                                    site[0]=xl;
                                                    if(par().inc==1){
                                                        sparseSite[0]=site[0];
                                                    }else if(xshift[tglb]!=step-1){
                                                        sparseSite[0]=2*int(site[0]/step) + (site[0]+xshift[tglb])%2;
                                                    }else{
                                                        sparseSite[0]=2*int(site[0]/step) + (site[0])%2;
                                                    }
                                                    for(int that=0;that<2;that++){
                                                        
                                                        site[3]=t+that;
                                                        sparseSite[3]=site[3];
                                                        
                                                        if(mu==0){// do v once
                                                            peekLocalSite(vec,temp,site);
                                                            pokeLocalSite(vec,v[il],sparseSite);
                                                            peekLocalSite(vec,temp2,site);
                                                            pokeLocalSite(vec,w0[il],sparseSite);
                                                        }else if(mu==1){
                                                            peekLocalSite(vec,temp2,site);
                                                            pokeLocalSite(vec,w1[il],sparseSite);
                                                        }else if(mu==2){
                                                            peekLocalSite(vec,temp2,site);
                                                            pokeLocalSite(vec,w2[il],sparseSite);
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            });
        }// end mu
    }// end evecs
    
    std::string dir = dirname(par().output);
    if (!par().output.empty())
    {
        int status = mkdir(dir);
        if (status)
        {
            HADRONS_ERROR(Io, "cannot create directory '" + dir
                          + "' ( " + std::strerror(errno) + ")");
        }
        startTimer("V I/O");
        A2AVectorsIo::write(par().output + "_v", v, par().multiFile, vm().getTrajectory());
        stopTimer("V I/O");
        startTimer("W I/O");
        A2AVectorsIo::write(par().output + "_w0", w0, par().multiFile, vm().getTrajectory());
        A2AVectorsIo::write(par().output + "_w1", w1, par().multiFile, vm().getTrajectory());
        A2AVectorsIo::write(par().output + "_w2", w2, par().multiFile, vm().getTrajectory());
        stopTimer("W I/O");
    }
    
    if ( env().getGrid()->IsBoss() ) {
        
        std::string eval_filename;
        if (!par().output.empty()){
            eval_filename=A2AVectorsIo::evalFilename(par().output,vm().getTrajectory());
        }else{
            eval_filename=A2AVectorsIo::evalFilename("evals",vm().getTrajectory());
        }
        A2AVectorsIo::initEvalFile(eval_filename,
                                   evalM.size());// total size
        A2AVectorsIo::saveEvalBlock(eval_filename,
                                    evalM.data(),
                                    0,// start of chunk
                                    2*Nl_);// size of chunk saved
    }
    qlat::clear_shuffled_fields_reader_cache();
}

#endif // USE_QLATTICCE

// Staggered A2A vectors

template <typename FImpl, typename Pack>
class TStagA2AVectors : public Module<A2AVectorsPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SOLVER_TYPE_ALIASES(FImpl,);
    typedef A2AVectorsSchurStaggered<FImpl> A2A;
public:
    // constructor
    TStagA2AVectors(const std::string name);
    // destructor
    virtual ~TStagA2AVectors(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    std::string  solverName_;
    unsigned int Nl_{0};
};

MODULE_REGISTER_TMP(StagA2AVectors,
                    ARG(TStagA2AVectors<STAGIMPL, BaseFermionEigenPack<STAGIMPL>>),
                    MSolver);

/******************************************************************************
 *                       TStagA2AVectors implementation                           *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
TStagA2AVectors<FImpl, Pack>::TStagA2AVectors(const std::string name)
: Module<A2AVectorsPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
std::vector<std::string> TStagA2AVectors<FImpl, Pack>::getInput(void)
{
    std::string              sub_string;
    std::vector<std::string> in;
    
    if (!par().eigenPack.empty())
    {
        in.push_back(par().eigenPack);
        sub_string = (!par().eigenPack.empty()) ? "_subtract" : "";
    }
    in.push_back(par().solver);
    in.push_back(par().noise);
    
    return in;
}

template <typename FImpl, typename Pack>
std::vector<std::string> TStagA2AVectors<FImpl, Pack>::getOutput(void)
{
    std::vector<std::string> out = {getName() + "_v", getName() + "_w"};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
void TStagA2AVectors<FImpl, Pack>::setup(void)
{
    bool        hasLowModes = (!par().eigenPack.empty());
    std::string sub_string  = (hasLowModes) ? "_subtract" : "";
    auto        &noise      = envGet(ColorDiagonalNoise<FImpl>, par().noise);
    auto        &action     = envGet(FMat, par().action);
    auto        &solver     = envGet(Solver, par().solver);
    int         Ls          = env().getObjectLs(par().action);
    
    
    if (hasLowModes)
    {
        auto &epack = envGet(Pack, par().eigenPack);
        Nl_ = epack.evec.size();
    }
    envCreate(std::vector<FermionField>, getName() + "_v", 1,
              2*Nl_ + noise.fermSize(), envGetGrid(FermionField));
    envCreate(std::vector<FermionField>, getName() + "_w", 1,
              2*Nl_ + noise.fermSize(), envGetGrid(FermionField));
    if (Ls > 1)
    {
        envTmpLat(FermionField, "f5", Ls);
    }
    envTmp(A2A, "a2a", 1, action, solver);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
void TStagA2AVectors<FImpl, Pack>::execute(void)
{
    auto        &action    = envGet(FMat, par().action);
    auto        &solver    = envGet(Solver, par().solver);
    auto        &noise     = envGet(ColorDiagonalNoise<FImpl>, par().noise);
    auto        &v         = envGet(std::vector<FermionField>, getName() + "_v");
    auto        &w         = envGet(std::vector<FermionField>, getName() + "_w");
    int         Ls         = env().getObjectLs(par().action);
    double      mass       = par().mass;
    
    envGetTmp(A2A, a2a);
    
    if (Nl_ > 0)
    {
        LOG(Message) << "Computing all-to-all vectors "
        << " using eigenpack '" << par().eigenPack << "' ("
        << 2*Nl_ << " low modes) and noise '"
        << par().noise << "' (" << noise.fermSize()
        << " noise vectors)" << std::endl;
    }
    else
    {
        LOG(Message) << "Computing all-to-all vectors "
        << " using noise '" << par().noise << "' (" << noise.size()
        << " noise vectors)" << std::endl;
    }
    // Low modes
    auto &epack  = envGet(Pack, par().eigenPack);
    for (unsigned int il = 0; il < Nl_; il++)
    {
        // eval of unpreconditioned Dirac op
        std::complex<double> eval(mass,sqrt(epack.eval[il]-mass*mass));
        
        startTimer("V low mode");
        LOG(Message) << "V vector i = " << 2*il << ", " << 2*il+1 << " (low modes)" << std::endl;
        if (Ls == 1)
        {
            a2a.makeLowModeV(v[2*il], epack.evec[il], eval);
            // construct -lambda evec
            a2a.makeLowModeV(v[2*il+1], epack.evec[il], eval, 1);
        }
        else
        {
            envGetTmp(FermionField, f5);
            a2a.makeLowModeV5D(v[2*il], f5, epack.evec[il], eval);
            // construct -lambda evec
            a2a.makeLowModeV5D(v[2*il+1], f5, epack.evec[il], eval, 1);
        }
        stopTimer("V low mode");
        startTimer("W low mode");
        LOG(Message) << "W vector i = " << 2*il << ", " << 2*il+1 << " (low modes)" << std::endl;
        if (Ls == 1)
        {
            a2a.makeLowModeW(w[2*il], epack.evec[il], eval);
            // construct -lambda evec
            a2a.makeLowModeW(w[2*il+1], epack.evec[il], eval, 1);
        }
        else
        {
            envGetTmp(FermionField, f5);
            a2a.makeLowModeW5D(w[2*il], f5, epack.evec[il], eval);
            // construct -lambda evec
            a2a.makeLowModeW5D(w[2*il+1], f5, epack.evec[il], eval, 1);
        }
        stopTimer("W low mode");
    }
    
    // High modes
#if 1
    
    FermionField sub(env().getGrid());
    FermionField noise_ih(env().getGrid());
    
    for (unsigned int ih = 0; ih < noise.fermSize(); ih++)
    {
        startTimer("V high mode");
        LOG(Message) << "V vector i = " << 2*Nl_ + ih
        << " (" << ((Nl_ > 0) ? "high " : "")
        << "stochastic mode)" << std::endl;
        if (Ls == 1)
        {
            a2a.makeHighModeV(v[2*Nl_ + ih], noise.getFerm(ih));
            // subtract the low modes
            sub = Zero();
            noise_ih = noise.getFerm(ih);
            for (int i=0;i<2*Nl_;i++) {
              const FermionField& tmp = v[i];
              // eval of unpreconditioned Dirac op
              std::complex<double> eval(mass,sqrt(epack.eval[i/2]-mass*mass));
              eval = i%2 ? eval : conjugate(eval);
              // need to subtract |l><l|/lambda_l |src>, so
              // mult, do not / by eval since v already has 1/eval
              axpy(sub,TensorRemove(innerProduct(tmp,noise_ih)) * eval,tmp,sub);
            }
            v[2*Nl_ + ih] -= sub;
        }
        else
        {
            envGetTmp(FermionField, f5);
            a2a.makeHighModeV5D(v[2*Nl_ + ih], f5, noise.getFerm(ih));
        }
        std::cout << "v high " << ih << std::endl;
        //std::cout << v[2*Nl_+ih] << std::endl;
        
        stopTimer("V high mode");
        startTimer("W high mode");
        LOG(Message) << "W vector i = " << 2*Nl_ + ih
        << " (" << ((Nl_ > 0) ? "high " : "")
        << "stochastic mode)" << std::endl;
        if (Ls == 1)
        {
            a2a.makeHighModeW(w[2*Nl_ + ih], noise.getFerm(ih));
        }
        else
        {
            envGetTmp(FermionField, f5);
            a2a.makeHighModeW5D(w[2*Nl_ + ih], f5, noise.getFerm(ih));
        }
        stopTimer("W high mode");
    }
#endif
    // I/O if necessary
    if (!par().output.empty())
    {
        startTimer("V I/O");
        A2AVectorsIo::write(par().output + "_v", v, par().multiFile, vm().getTrajectory());
        stopTimer("V I/O");
        startTimer("W I/O");
        A2AVectorsIo::write(par().output + "_w", w, par().multiFile, vm().getTrajectory());
        stopTimer("W I/O");
    }
    //printMem("End StagA2AVectors execute() ", env().getGrid()->ThisRank());
}
END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSolver_A2AVectors_hpp_
