/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/A2AMatrix.hpp

Copyright (C) 2015-2018

Author: Antonin Portelli <antonin.portelli@me.com>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

See the full license in the file "LICENSE" in the top level distribution directory
*************************************************************************************/
/*  END LEGAL */
#ifndef A2A_Matrix_Nucleon_hpp_
#define A2A_Matrix_Nucleon_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/TimerArray.hpp>
#include <Grid/Eigen/unsupported/CXX11/Tensor>
#ifdef USE_MKL
#include "mkl.h"
#include "mkl_cblas.h"
#endif

#ifndef HADRONS_A2AN_NAME 
#define HADRONS_A2AN_NAME "a2aMatrixNucleon"
#endif

// precision to store A2A nucleon field data as
#ifndef HADRONS_A2AN_IO_TYPE
#define HADRONS_A2AN_IO_TYPE ComplexF
#endif

// precision used for 
#ifndef HADRONS_A2AN_CALC_TYPE
#define HADRONS_A2AN_CALC_TYPE ComplexD
#endif

#define HADRONS_A2AN_PARALLEL_IO

BEGIN_HADRONS_NAMESPACE

// general A2A matrix set based on Eigen tensors and Grid-allocated memory
// Dimensions:
//   0 - ext - external field (momentum, EM field, ...)
//   1 - str - spin-color structure
//   2 - t   - timeslice
//   3 - i   - left  A2A mode index
//   4 - j   - right A2A mode index
//   5 - k   - q3 A2A mode index
template <typename T>
using A2AMatrixSetNuc = Eigen::TensorMap<Eigen::Tensor<T, 6, Eigen::RowMajor>>;

// MCA -- testing - using rank 4 tensor here for mu(spin),i,j,k notation, using Tensor to be able to resize
template <typename T>
using A2AMatrixNuc = Eigen::Tensor<T,4,Eigen::RowMajor>;

template <typename T>
using A2AMatrixNucTr = Eigen::Tensor<T,4,Eigen::ColMajor>;

// MCA - Keeping the old calls just in case I need them
template <typename T>
using A2AMatrix = Eigen::Matrix<T, -1, -1, Eigen::RowMajor>;

template <typename T>
using A2AMatrixTr = Eigen::Matrix<T, -1, -1, Eigen::ColMajor>;

/******************************************************************************
 *                      Abstract class for A2A kernels                        *
 ******************************************************************************/
template <typename T, typename Field>
class A2AKernelNucleon
{
public:
    A2AKernelNucleon(void) = default;
    virtual ~A2AKernelNucleon(void) = default;
    virtual void operator()(A2AMatrixSetNuc<T> &m, const Field *left, const Field *right, const Field *q3,
                          const unsigned int orthogDim, double &time) = 0;
    virtual double flops(const unsigned int blockSizei, const unsigned int blockSizej, const unsigned int blockSizek) = 0;
    virtual double bytes(const unsigned int blockSizei, const unsigned int blockSizej, const unsigned int blockSizek) = 0;
};

// MCA - I will need to update most of this part
/******************************************************************************
 *                  Class to handle A2A matrix block HDF5 I/O                 *
 ******************************************************************************/
template <typename T>
class A2AMatrixNucIo
{
public:
    // constructors
    A2AMatrixNucIo(void) = default;
    A2AMatrixNucIo(std::string filename, std::string dataname, 
                const unsigned int nt, const unsigned int ni = 0,
                const unsigned int nj = 0, const unsigned int nk = 0);
    // destructor
    ~A2AMatrixNucIo(void) = default;
    // access
    unsigned int getNi(void) const;
    unsigned int getNj(void) const;
    unsigned int getNk(void) const;
    unsigned int getNt(void) const;
    size_t       getSize(void) const;
    // file allocation
    template <typename MetadataType>
    void initFile(const MetadataType &d, const unsigned int chunkSize);
    // block I/O
    void saveBlock(const T *data, const unsigned int str, const unsigned int i, const unsigned int j, const unsigned int k,
                   const unsigned int blockSizei, const unsigned int blockSizej, const unsigned int blockSizek);
    void saveBlock(const A2AMatrixSetNuc<T> &m, const unsigned int ext, const unsigned int str,
                   const unsigned int i, const unsigned int j, const unsigned int k);
    template <template <class> class Vec, typename VecT>
    void load(Vec<VecT> &v, double *tRead = nullptr);
    template <template <class> class Vec, typename VecT>
    void load(Vec<VecT> &v, GridCartesian *grid, double *tRead = nullptr);
private:
    std::string  filename_{""}, dataname_{""};
    unsigned int nt_{0}, ni_{0}, nj_{0}, nk_{0};
};

/******************************************************************************
 *                  Wrapper for A2A matrix block computation                  *
 ******************************************************************************/
template <typename T, typename Field, typename MetadataType, typename TIo = T>
class A2AMatrixNucleonBlockComputation
{
private:
    struct IoHelper
    {
        A2AMatrixNucIo<TIo> io;
        MetadataType     md;
        unsigned int     e, s, i, j, k;
    };
    typedef std::function<std::string(const unsigned int)>  FilenameFn;
    typedef std::function<MetadataType(const unsigned int)> MetadataFn;
public:
    // constructor
    A2AMatrixNucleonBlockComputation(GridBase *grid,
                              const unsigned int orthogDim,
                              const unsigned int next,
                              const unsigned int nstr,
                              const unsigned int blockSize,
                              const unsigned int cacheBlockSize,
                              TimerArray *tArray = nullptr);
    // execution
    void execute(const std::vector<Field> &left, 
                 const std::vector<Field> &right,
                 const std::vector<Field> &q3,
                 A2AKernelNucleon<T, Field> &kernel,
                 const FilenameFn &ionameFn,
                 const FilenameFn &filenameFn,
                 const MetadataFn &metadataFn);
private:
    // I/O handler
    void saveBlock(const A2AMatrixSetNuc<TIo> &m, IoHelper &h);
private:
    TimerArray            *tArray_;
    GridBase              *grid_;
    unsigned int          orthogDim_, nt_, next_, nstr_, blockSize_, cacheBlockSize_;
    Vector<T>             mCache_;
    Vector<TIo>           mBuf_;
    std::vector<IoHelper> nodeIo_;
};

/******************************************************************************
 *                       A2A matrix contraction kernels                       *
 ******************************************************************************/
class A2AContractionNucleon
{
public:
	// function to contract two rank 4 A2A nucleon fields into a spin matrix
	// spinMat_(mu,nu) = a_(ijk,mu)*b(ijk,nu)
	template <typename Mat, typename TenLeft, typename TenRight>
	static inline void contNucTen(Mat &spinMat, const TenLeft &a, const TenRight &b)
	{
		//std::cout << "DEBUG: Nucleon contraction 1 started" << std::endl;

		Eigen::Index nRow_mu, nRow_nu, nRow_i, nRow_j, nRow_k;
		
		// Make sure spinMat is a square matrix
		if ((a.dimension(0) != b.dimension(0)) or (spinMat.rows() != spinMat.cols()) or (a.dimension(0) != spinMat.rows()))
		{
			// error here -- dimensions mismatch
		}
		else
		{
			nRow_mu = a.dimension(0);
			nRow_nu = b.dimension(0);
		}
		
		// Make sure a and b dimensions match
		if ((a.dimension(1) != b.dimension(1)) or (a.dimension(2) != b.dimension(2)) or (a.dimension(3) != b.dimension(3)))
		{
			// error here -- index dimensions mismatch
		}
		else
		{
			nRow_i = a.dimension(1);
			nRow_j = a.dimension(2);
			nRow_k = a.dimension(3);
		}
		
	       // MCA - very basic parallelization for now
               thread_for(mu, Ns, {
                   for(int nu=0; nu < Ns; nu++)
                      for (int i = 0; i < nRow_i; i++)
                        for (int j = 0; j < nRow_j; j++)
                          for (int k = 0; k < nRow_k; k++)
                             {
                                //ComplexD tmp;
                                //tmp = a(mu,i,j,k) * (std::conj(b(nu,i,j,k)) - std::conj(b(nu,k,j,i)));
                                //ComplexD tmp2(b(nu,i,j,k).real()-b(nu,k,j,i).real(),b(nu,k,j,i).imag()-b(nu,i,j,k).imag());
                                //tmp = a(mu,i,j,k) * tmp2;
                                //spinMat(mu,nu) += tmp;
                                spinMat(mu,nu) += a(mu,i,j,k) * (std::conj(b(nu,i,j,k)) - std::conj(b(nu,k,j,i)));
                             }
               });
	}

	static void contNucTenTest()
	{
		std::cout << "DEBUG: Starting test of nucleon contraction 1" << std::endl;

		static const ComplexD one(1., 0.), zero(0., 0.);
		static const std::complex<double> I(0., 1.);

		int imax = 50;
		int jmax = 50;
		int kmax = 50;
		
		std::cout << "Constructing rank 2 dummy tensor" << std::endl;
		A2AMatrix<ComplexD> spinMat(Ns,Ns);
		
		std::cout << "Constructing rank 4 dummy tensors" << std::endl;
		A2AMatrixNuc<ComplexD> A(Ns,imax,jmax,kmax);
		A2AMatrixNuc<ComplexD> B(Ns,imax,jmax,kmax);
		
		std::cout << "Initializing dummy tensors" << std::endl;
		// initialize
		for (int mu = 0; mu < Ns; mu++)
		for (int i = 0; i < imax; i++)
		for (int j = 0; j < jmax; j++)
		for (int k = 0; k < kmax; k++)
		{
			A(mu, i, j, k) = 1. + (1.0*I);
			B(mu, i, j, k) = 1. + (1.0*I);
		}
		std::cout << "Initialized successfully" << std::endl;
		
		std::cout << "Initializing result matrix" << std::endl;
		for (int mu = 0; mu < Ns; mu++)
		for (int nu = 0; nu < Ns; nu++)
		{
			spinMat(mu,nu) = zero;
		}
		std::cout << "Initialized successfully" << std::endl;
		
		std::cout << A(0,1,2,3) << std::endl;
		
		std::cout << "Calling contraction 1 function" << std::endl;
		A2AContractionNucleon::contNucTen(spinMat, A, B);
		std::cout << "Successfully called" << std::endl;
		
		for (int mu = 0; mu < Ns; mu++)
		{
			for (int nu = 0; nu < Ns; nu++)
			{
				std::cout << spinMat(mu,nu) << " ";
			}
			std::cout << std::endl;
		}
		
		std::cout << "DEBUG: Finished test of nucleon contraction 1" << std::endl;
	}

	// function to contract two rank 4 A2A nucleon fields and a rank 2 meson field into a spin matrix
	// spinMat_(mu,nu) = a_(ijk,mu)*c(jm)*b(imk,nu)
	template <typename Mat, typename TenLeft, typename MatMid, typename TenRight>
	static inline void contNuc3ptUp(Mat &spinMat, const TenLeft &a, const MatMid &c, const TenRight &b)
	{

		Eigen::Index nRow_mu, nRow_nu, nRow_i, nRow_j, nRow_k, nRow_m;
		
		// Make sure spinMat is a square matrix
		if ((a.dimension(0) != b.dimension(0)) or (spinMat.rows() != spinMat.cols()) or (a.dimension(0) != spinMat.rows()))
		{
			// error here -- dimensions mismatch
		}
		else
		{
			nRow_mu = a.dimension(0);
			nRow_nu = b.dimension(0);
		}
		
		// Make sure a and b dimensions match
		if ((a.dimension(1) != b.dimension(1)) or (a.dimension(2) != b.dimension(2)) or (a.dimension(3) != b.dimension(3)))
		{
			// error here -- index dimensions mismatch
		}
		else
		{
			nRow_i = a.dimension(1);
			nRow_j = a.dimension(2);
			nRow_k = a.dimension(3);
		}
		
		// Make sure c is a square matrix
		if ((c.rows() != c.cols()) or (c.rows() != nRow_i))
		{
			// error here -- dimensions mismatch
		}
		else
		{
			nRow_m = c.rows();
		}
		
		// MCA - very basic parallelization for now
                thread_for(mu,Ns,{
                       thread_for(nu,Ns,{
                                for (int i = 0; i < nRow_i; i++)
                                for (int j = 0; j < nRow_j; j++)
                                for (int k = 0; k < nRow_k; k++)
                                for (int m = 0; m < nRow_m; m++)
                                {
                                        ComplexD tmp;
                                        tmp = a(mu,i,j,k) * c(k,m) * std::conj(b(nu,i,j,m))
                                                        - a(mu,i,j,k) * c(k,m) * std::conj(b(nu,m,j,i))
                                                        + a(mu,k,j,i) * c(k,m) * std::conj(b(nu,m,j,i))
                                                        - a(mu,k,j,i) * c(k,m) * std::conj(b(nu,i,j,m));
                                        spinMat(mu,nu) += tmp;

                                }
                        });
                });
	}

	static void contNuc3ptUpTest()
	{
		std::cout << "DEBUG: Starting test of nuc3pt up contraction 1" << std::endl;

		static const ComplexD one(1., 0.), zero(0., 0.);
		static const std::complex<double> I(0., 1.);

		int imax = 50;
		int jmax = 50;
		int kmax = 50;
		int mmax = 50;
		
		std::cout << "Constructing rank 2 dummy tensors" << std::endl;
		A2AMatrix<ComplexD> spinMat(Ns,Ns);
		A2AMatrix<ComplexD> C(jmax, mmax);
		
		std::cout << "Constructing rank 4 dummy tensors" << std::endl;
		A2AMatrixNuc<ComplexD> A(Ns,imax,jmax,kmax);
		A2AMatrixNuc<ComplexD> B(Ns,imax,jmax,kmax);
		
		std::cout << "Initializing dummy tensors" << std::endl;
		// initialize
		for (int mu = 0; mu < Ns; mu++)
		for (int i = 0; i < imax; i++)
		for (int j = 0; j < jmax; j++)
		for (int k = 0; k < kmax; k++)
		{
			A(mu, i, j, k) = 1. + (1. * I);
			B(mu, i, j, k) = 1. + (1. * I);
		}
		
		for (int j = 0; j < jmax; j++)
		for (int m = 0; m < mmax; m++)
		{
			C(j, m) = 1. + (1. * I);
		}
		std::cout << "Initialized successfully" << std::endl;
		
		std::cout << "Initializing result matrix" << std::endl;
		for (int mu = 0; mu < Ns; mu++)
		for (int nu = 0; nu < Ns; nu++)
		{
			spinMat(mu,nu) = zero;
		}
		std::cout << "Initialized successfully" << std::endl;
		
		std::cout << A(0,1,2,3) << std::endl;
		
		std::cout << "Calling contraction 1 function" << std::endl;
		A2AContractionNucleon::contNuc3ptUp(spinMat, A, C, B);
		std::cout << "Successfully called" << std::endl;
		
		for (int mu = 0; mu < Ns; mu++)
		{
			for (int nu = 0; nu < Ns; nu++)
			{
				std::cout << spinMat(mu,nu) << " ";
			}
			std::cout << std::endl;
		}
		
		std::cout << "DEBUG: Finished test of nuc3pt contraction 1" << std::endl;
	}

	// MCA - TO DO? (may remove)
	// TrProjMat(acc, p, m): acc = tr(p*m)
	// This is essentially the same as accTrMul but designed for use with projection
	// of spin matrices.
	template <typename C, typename MatLeft, typename MatRight>
	static inline void TrProjMat(C &acc, const MatLeft &p, const MatRight &m)
	{
		std::cout << "Testing TrProjMat function" << std::endl;
	}
	
	// MCA - TO DO
	// ProjTPlus: res = (1/2)(1+Gt)*a
	template <typename Mat>
	static inline void ProjTPlus(Mat &res, Mat &a)
	{

		// do projector manually here (1+Gt)
		A2AMatrix<ComplexD> projMat(Ns,Ns);
		projMat = A2AMatrix<ComplexD>::Identity(Ns,Ns);
		projMat.bottomLeftCorner(Ns/2,Ns/2) = A2AMatrix<ComplexD>::Identity(Ns/2,Ns/2);
		projMat.topRightCorner(Ns/2,Ns/2) = A2AMatrix<ComplexD>::Identity(Ns/2,Ns/2);

		res = 0.5 * (projMat*a);

	}

	// ProjTPlusPxy: res = (1/2)(1+Gt)(1-GxGy)*a
	template <typename Mat>
	static inline void ProjTPlusPxy(Mat &res, Mat &a)
	{

		// do projector manually here (1+Gt)(1-GxGy)
		A2AMatrix<ComplexD> projMat(Ns,Ns);
		A2AMatrix<ComplexD> subMat(Ns/2, Ns/2);
		
		subMat = A2AMatrix<ComplexD>::Zero(Ns/2,Ns/2);
		subMat(1,1) = 2;
		
		projMat = A2AMatrix<ComplexD>::Zero(Ns,Ns);
		projMat.topLeftCorner(Ns/2,Ns/2) = subMat;
		projMat.bottomLeftCorner(Ns/2,Ns/2) = subMat;
		projMat.topRightCorner(Ns/2,Ns/2) = subMat;
		projMat.bottomRightCorner(Ns/2,Ns/2) = subMat;

		res = projMat*a;

	}

	// TrProjTPlusPxy
	template <typename C, typename Mat>
	static inline void TrProjTPlusPxy(C &acc, Mat &a)
	{
		A2AMatrix<ComplexD> tmpMat(Ns,Ns);
		C tmp;
		
		ProjTPlusPxy(tmpMat,a);
		tmp = tmpMat.trace();
		
		//std::cout << tmp << std::endl;
		
		acc += tmp;
		
	}

	// testProjTPlusPxy

	template <typename C, typename Mat>
	static inline void TrProjTPlus(C &acc, Mat &a)
	{
		A2AMatrix<ComplexD> tmpMat(Ns,Ns);
		C tmp;
		
		ProjTPlus(tmpMat,a);
		tmp = tmpMat.trace();
		
		//std::cout << tmp << std::endl;
		
		acc += tmp;
		
	}

	// MCA - Test Function
	static inline void testProjTPlus()
	{
		std::cout << "DEBUG: Starting test of ProjTPlus" << std::endl;
		
		static const std::complex<double> I(0., 1.);

		//A2AMatrix<ComplexD> res(Ns,Ns);
		//res = A2AMatrix<ComplexD>::Zero(Ns,Ns);
		ComplexD total = 0.;
		
		A2AMatrix<ComplexD> dummyMat(Ns,Ns);
		for (int mu = 0; mu < Ns; mu++)
		for (int nu = 0; nu < Ns; nu++)
		{
			dummyMat(mu,nu) = 1. + (1. * I);
		}
				
		TrProjTPlus(total, dummyMat);
		
		std::cout << total << std::endl;
		
		std::cout << "DEBUG: Ending test of ProjTPlus" << std::endl;
	}

	static inline void testProjTPlusPxy()
	{
		std::cout << "DEBUG: Starting test of ProjTPlusPxy" << std::endl;
		
		static const std::complex<double> I(0., 1.);

		//A2AMatrix<ComplexD> res(Ns,Ns);
		//res = A2AMatrix<ComplexD>::Zero(Ns,Ns);
		ComplexD total = 0.;
		
		A2AMatrix<ComplexD> dummyMat(Ns,Ns);
		for (int mu = 0; mu < Ns; mu++)
		for (int nu = 0; nu < Ns; nu++)
		{
			dummyMat(mu,nu) = 1. + (1. * I);
		}
				
		TrProjTPlusPxy(total, dummyMat);
		
		std::cout << total << std::endl;
		
		std::cout << "DEBUG: Ending test of ProjTPlus" << std::endl;
	}

	template <typename C, typename TenLeft, typename TenRight>
	static inline void ContractNucleonTPlus(C &acc, const TenLeft &a, const TenRight &b)
	{
		// Temp matrix
		A2AMatrix<ComplexD> spinMat(Ns,Ns);
		
		// Initialize		
		for (int mu = 0; mu < Ns; mu++)
		for (int nu = 0; nu < Ns; nu++)
		{
			spinMat(mu,nu) = 0.;
		}

		// Do contraction (1)
		contNucTen(spinMat, a, b);
		
		// Take the PPAR projector and trace
		TrProjTPlus(acc, spinMat);
	}

	template <typename C, typename TenLeft, typename MatMid, typename TenRight>
	static inline void ContractNuc3ptTPlusPxy(C &acc, const TenLeft &a, const MatMid &c, const TenRight &b)
	{
		// Temp matrix
		A2AMatrix<ComplexD> spinMat(Ns,Ns);
		
		// Initialize		
		for (int mu = 0; mu < Ns; mu++)
		for (int nu = 0; nu < Ns; nu++)
		{
			spinMat(mu,nu) = 0.;
		}

		// Do contraction (1)
		contNuc3ptUp(spinMat, a, c, b);
		
		// Take the PPAR projector and trace
		TrProjTPlusPxy(acc, spinMat);
	}

    // accTrMul(acc, a, b): acc += tr(a*b)
    template <typename C, typename MatLeft, typename MatRight>
    static inline void accTrMul(C &acc, const MatLeft &a, const MatRight &b)
    {
        if ((MatLeft::Options == Eigen::RowMajor) and
            (MatRight::Options == Eigen::ColMajor))
        {
            thread_for(r,a.rows(),{
                C tmp;
#ifdef USE_MKL
                dotuRow(tmp, r, a, b);
#else
                tmp = a.row(r).conjugate().dot(b.col(r));
#endif
                thread_critical
                {
                    acc += tmp;
                }
            });
        }
        else
        {
            thread_for(c,a.cols(),{
                C tmp;
#ifdef USE_MKL 
                dotuCol(tmp, c, a, b);
#else
                tmp = a.col(c).conjugate().dot(b.row(c));
#endif
                thread_critical
                {
                    acc += tmp;
                }
            });
        }
    }

	// MCA - TO DO
	template <typename TenLeft, typename TenRight>
	static inline double contNucTenFlops(const TenLeft &a, const TenRight &b)
	{
		// Add flops data here
		
		return 0.;
	}

    template <typename MatLeft, typename MatRight>
    static inline double accTrMulFlops(const MatLeft &a, const MatRight &b)
    {
        double n = a.rows()*a.cols();

        return 8.*n;
    }

    // mul(res, a, b): res = a*b
#ifdef USE_MKL
    template <template <class, int...> class Mat, int... Opts>
    static inline void mul(Mat<ComplexD, Opts...> &res, 
                           const Mat<ComplexD, Opts...> &a, 
                           const Mat<ComplexD, Opts...> &b)
    {
        static const ComplexD one(1., 0.), zero(0., 0.);

        if ((res.rows() != a.rows()) or (res.cols() != b.cols()))
        {
            res.resize(a.rows(), b.cols());
        }
        if (Mat<ComplexD, Opts...>::Options == Eigen::RowMajor)
        {
            cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, a.rows(), b.cols(),
                        a.cols(), &one, a.data(), a.cols(), b.data(), b.cols(), &zero,
                        res.data(), res.cols());
        }
        else if (Mat<ComplexD, Opts...>::Options == Eigen::ColMajor)
        {
            cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, a.rows(), b.cols(),
                        a.cols(), &one, a.data(), a.rows(), b.data(), b.rows(), &zero,
                        res.data(), res.rows());
        }
    }

    template <template <class, int...> class Mat, int... Opts>
    static inline void mul(Mat<ComplexF, Opts...> &res, 
                           const Mat<ComplexF, Opts...> &a, 
                           const Mat<ComplexF, Opts...> &b)
    {
        static const ComplexF one(1., 0.), zero(0., 0.);

        if ((res.rows() != a.rows()) or (res.cols() != b.cols()))
        {
            res.resize(a.rows(), b.cols());
        }
        if (Mat<ComplexF, Opts...>::Options == Eigen::RowMajor)
        {
            cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, a.rows(), b.cols(),
                        a.cols(), &one, a.data(), a.cols(), b.data(), b.cols(), &zero,
                        res.data(), res.cols());
        }
        else if (Mat<ComplexF, Opts...>::Options == Eigen::ColMajor)
        {
            cblas_cgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, a.rows(), b.cols(),
                        a.cols(), &one, a.data(), a.rows(), b.data(), b.rows(), &zero,
                        res.data(), res.rows());
        }
    }
#else
    template <typename Mat>
    static inline void mul(Mat &res, const Mat &a, const Mat &b)
    {
        res = a*b;
    }
#endif
    template <typename Mat>
    static inline double mulFlops(const Mat &a, const Mat &b)
    {
        double nr = a.rows(), nc = a.cols();

        return nr*nr*(6.*nc + 2.*(nc - 1.));
    }
private:
    template <typename C, typename MatLeft, typename MatRight>
    static inline void makeDotRowPt(C * &aPt, unsigned int &aInc, C * &bPt, 
                                    unsigned int &bInc, const unsigned int aRow, 
                                    const MatLeft &a, const MatRight &b)
    {
        if (MatLeft::Options == Eigen::RowMajor)
        {
            aPt  = a.data() + aRow*a.cols();
            aInc = 1;
        }
        else if (MatLeft::Options == Eigen::ColMajor)
        {
            aPt  = a.data() + aRow;
            aInc = a.rows();
        }
        if (MatRight::Options == Eigen::RowMajor)
        {
            bPt  = b.data() + aRow;
            bInc = b.cols();
        }
        else if (MatRight::Options == Eigen::ColMajor)
        {
            bPt  = b.data() + aRow*b.rows();
            bInc = 1;
        }
    }

#ifdef USE_MKL
    template <typename C, typename MatLeft, typename MatRight>
    static inline void makeDotColPt(C * &aPt, unsigned int &aInc, C * &bPt, 
                                    unsigned int &bInc, const unsigned int aCol, 
                                    const MatLeft &a, const MatRight &b)
    {
        if (MatLeft::Options == Eigen::RowMajor)
        {
            aPt  = a.data() + aCol;
            aInc = a.cols();
        }
        else if (MatLeft::Options == Eigen::ColMajor)
        {
            aPt  = a.data() + aCol*a.rows();
            aInc = 1;
        }
        if (MatRight::Options == Eigen::RowMajor)
        {
            bPt  = b.data() + aCol*b.cols();
            bInc = 1;
        }
        else if (MatRight::Options == Eigen::ColMajor)
        {
            bPt  = b.data() + aCol;
            bInc = b.rows();
        }
    }

    template <typename MatLeft, typename MatRight>
    static inline void dotuRow(ComplexF &res, const unsigned int aRow,
                               const MatLeft &a, const MatRight &b)
    {
        const ComplexF *aPt, *bPt;
        unsigned int   aInc, bInc;

        makeDotRowPt(aPt, aInc, bPt, bInc, aRow, a, b);
        cblas_cdotu_sub(a.cols(), aPt, aInc, bPt, bInc, &res);
    }

    template <typename MatLeft, typename MatRight>
    static inline void dotuCol(ComplexF &res, const unsigned int aCol,
                               const MatLeft &a, const MatRight &b)
    {
        const ComplexF *aPt, *bPt;
        unsigned int   aInc, bInc;

        makeDotColPt(aPt, aInc, bPt, bInc, aCol, a, b);
        cblas_cdotu_sub(a.rows(), aPt, aInc, bPt, bInc, &res);
    }

    template <typename MatLeft, typename MatRight>
    static inline void dotuRow(ComplexD &res, const unsigned int aRow,
                               const MatLeft &a, const MatRight &b)
    {
        const ComplexD *aPt, *bPt;
        unsigned int   aInc, bInc;

        makeDotRowPt(aPt, aInc, bPt, bInc, aRow, a, b);
        cblas_zdotu_sub(a.cols(), aPt, aInc, bPt, bInc, &res);
    }

    template <typename MatLeft, typename MatRight>
    static inline void dotuCol(ComplexD &res, const unsigned int aCol,
                               const MatLeft &a, const MatRight &b)
    {
        const ComplexD *aPt, *bPt;
        unsigned int   aInc, bInc;

        makeDotColPt(aPt, aInc, bPt, bInc, aCol, a, b);
        cblas_zdotu_sub(a.rows(), aPt, aInc, bPt, bInc, &res);
    }
#endif
};

/******************************************************************************
 *                     A2AMatrixIo template implementation                    *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename T>
A2AMatrixNucIo<T>::A2AMatrixNucIo(std::string filename, std::string dataname, 
                            const unsigned int nt, const unsigned int ni,
                            const unsigned int nj, const unsigned int nk)
: filename_(filename), dataname_(dataname)
, nt_(nt), ni_(ni), nj_(nj), nk_(nk)
{}

// access //////////////////////////////////////////////////////////////////////
template <typename T>
unsigned int A2AMatrixNucIo<T>::getNt(void) const
{
    return nt_;
}

template <typename T>
unsigned int A2AMatrixNucIo<T>::getNi(void) const
{
    return ni_;
}

template <typename T>
unsigned int A2AMatrixNucIo<T>::getNj(void) const
{
    return nj_;
}

template <typename T>
unsigned int A2AMatrixNucIo<T>::getNk(void) const
{
    return nk_;
}

template <typename T>
size_t A2AMatrixNucIo<T>::getSize(void) const
{
    return static_cast<size_t>(Ns)*
           static_cast<size_t>(nt_)*
           static_cast<size_t>(ni_)*
           static_cast<size_t>(nj_)*
           static_cast<size_t>(nk_)*
           static_cast<size_t>(sizeof(T));
}

// file allocation /////////////////////////////////////////////////////////////
template <typename T>
template <typename MetadataType>
void A2AMatrixNucIo<T>::initFile(const MetadataType &d, const unsigned int chunkSize)
{
#ifdef HAVE_HDF5
    std::vector<hsize_t>    dim = {static_cast<hsize_t>(nt_),
			                       static_cast<hsize_t>(Ns), 
                                   static_cast<hsize_t>(ni_), 
                                   static_cast<hsize_t>(nj_),
                                   static_cast<hsize_t>(nk_)},
                            chunk = {static_cast<hsize_t>(1),
								     static_cast<hsize_t>(Ns), 
                                     static_cast<hsize_t>(chunkSize), 
                                     static_cast<hsize_t>(chunkSize),
                                     static_cast<hsize_t>(chunkSize)};
    H5NS::DataSpace         dataspace(dim.size(), dim.data());
    H5NS::DataSet           dataset;
    H5NS::DSetCreatPropList plist;
    
    // create empty file just with metadata
    {
        Hdf5Writer writer(filename_);
        write(writer, dataname_, d);
    }

    // create the dataset
    Hdf5Reader reader(filename_, false);

    push(reader, dataname_);
    auto &group = reader.getGroup();
    plist.setChunk(chunk.size(), chunk.data());
    plist.setFletcher32();
    dataset = group.createDataSet(HADRONS_A2AN_NAME, Hdf5Type<T>::type(), dataspace, plist);
#else
    HADRONS_ERROR(Implementation, "all-to-all matrix I/O needs HDF5 library");
#endif
}

// block I/O ///////////////////////////////////////////////////////////////////
template <typename T>
void A2AMatrixNucIo<T>::saveBlock(const T *data, 
							   const unsigned int str,
                               const unsigned int i, 
                               const unsigned int j,
                               const unsigned int k,
                               const unsigned int blockSizei,
                               const unsigned int blockSizej,
                               const unsigned int blockSizek)
{
#ifdef HAVE_HDF5
    Hdf5Reader           reader(filename_, false);
    std::vector<hsize_t> count = {nt_, Ns, blockSizei, blockSizej, blockSizek},
                         offset = {0, static_cast<hsize_t>(str), static_cast<hsize_t>(i),
                                   static_cast<hsize_t>(j),
                                   static_cast<hsize_t>(k)},
                         stride = {1, 1, 1, 1, 1},
                         block  = {1, 1, 1, 1, 1}; 
    H5NS::DataSpace      memspace(count.size(), count.data()), dataspace;
    H5NS::DataSet        dataset;
    size_t               shift;

    push(reader, dataname_);
    auto &group = reader.getGroup();
    dataset     = group.openDataSet(HADRONS_A2AN_NAME);
    dataspace   = dataset.getSpace();
    dataspace.selectHyperslab(H5S_SELECT_SET, count.data(), offset.data(),
                              stride.data(), block.data());
    dataset.write(data, Hdf5Type<T>::type(), memspace, dataspace);
#else
    HADRONS_ERROR(Implementation, "all-to-all matrix I/O needs HDF5 library");
#endif
}

template <typename T>
void A2AMatrixNucIo<T>::saveBlock(const A2AMatrixSetNuc<T> &m,
                               const unsigned int ext, const unsigned int str,
                               const unsigned int i, const unsigned int j, const unsigned int k)
{
    unsigned int blockSizei = m.dimension(3);
    unsigned int blockSizej = m.dimension(4);
    unsigned int blockSizek = m.dimension(5);
    unsigned int nstr       = m.dimension(1);
    size_t       offset     = (ext*nstr + str)*nt_*blockSizei*blockSizej*blockSizek;

	std::cout << "Still working inside NucIo::saveBlock with external " << std::to_string(ext) << " and spin " << std::to_string(str)
				<< " and offset " << std::to_string(offset) << std::endl;
    saveBlock(m.data() + offset, str, i, j, k, blockSizei, blockSizej, blockSizek);
}

template <typename T>
template <template <class> class Vec, typename VecT>
void A2AMatrixNucIo<T>::load(Vec<VecT> &v, double *tRead)
{
#ifdef HAVE_HDF5
    Hdf5Reader           reader(filename_);
    std::vector<hsize_t> hdim;
    H5NS::DataSet        dataset;
    H5NS::DataSpace      dataspace;
    H5NS::CompType       datatype;
    H5NS::DSetCreatPropList plist;
    
    push(reader, dataname_);
    auto &group = reader.getGroup();
    dataset     = group.openDataSet(HADRONS_A2AN_NAME);
    plist       = dataset.getCreatePlist();
    datatype    = dataset.getCompType();
    dataspace   = dataset.getSpace();
    hdim.resize(dataspace.getSimpleExtentNdims());
    dataspace.getSimpleExtentDims(hdim.data());
        
    // chunking info
    bool isChunked = false;
    hsize_t chunk_dim[5];
    int rank_chunk;
    unsigned int chunkCount = 1;
    
    // adjust the chunk cache to try to improve performance
    hid_t dapl_id = H5Dget_access_plist(dataset.getId());
    hsize_t max_objects_in_chunk_cache = 12799;
    hsize_t max_bytes_in_cache = 128*1024*1024;
    
    H5Pset_chunk_cache(dapl_id, max_objects_in_chunk_cache, max_bytes_in_cache, H5D_CHUNK_CACHE_W0_DEFAULT);
    
    if (plist.getLayout() == H5D_CHUNKED)
    {
        isChunked = true;
        rank_chunk = plist.getChunk(5, chunk_dim);
        std::cout << "Data is chunked with rank " << std::to_string(rank_chunk) << std::endl;
        std::cout << "and chunk dimensions " << std::to_string(chunk_dim[0]) << " " <<
                                                std::to_string(chunk_dim[1]) << " " <<
                                                std::to_string(chunk_dim[2]) << " " <<
                                                std::to_string(chunk_dim[3]) << " " <<
                                                std::to_string(chunk_dim[4]) << std::endl;
    }
    
    std::cout << "DEBUG: Ns nt_ ni_ nj_ nk_ are " << std::to_string(Ns) << " " << std::to_string(nt_) << " " << std::to_string(ni_) << " "
			  << std::to_string(nj_) << " " << std::to_string(nk_) << std::endl;
	std::cout << "hdims are " <<std::to_string(hdim[0]) << " " << std::to_string(hdim[1]) << " " << std::to_string(hdim[2]) << " "
			  << std::to_string(hdim[3]) << " " << std::to_string(hdim[4]) << std::endl;
    if ((Ns*nt_*ni_*nj_*nk_ != 0) and
        ((hdim[0] != Ns) or (hdim[1] != nt_) or (hdim[2] != ni_) or (hdim[3] != nj_) or (hdim[4] != nk_)))
    {
        HADRONS_ERROR(Size, "all-to-all matrix size mismatch (got "
            + std::to_string(hdim[0]) + "x" + std::to_string(hdim[1]) + "x"
            + std::to_string(hdim[2]) + "x" + std::to_string(hdim[3]) + "x"
            + std::to_string(hdim[4]) + ", expected "
            + std::to_string(Ns) + "x"
            + std::to_string(nt_) + "x" + std::to_string(ni_) + "x"
            + std::to_string(nj_) + "x" + std::to_string(nk_));
    }
    else if (ni_*nj_*nk_ == 0)
    {
        if ((hdim[0] != nt_) or (hdim[1] != Ns))
        {
            HADRONS_ERROR(Size, "all-to-all time size mismatch (got "
                + std::to_string(hdim[0]) + "x"
                + std::to_string(hdim[1]) + ", expected "
                + std::to_string(nt_) + "x"
                + std::to_string(Ns) + ")");
        }
        ni_ = hdim[2];
        nj_ = hdim[3];
        nk_ = hdim[4];
    }
    
    // check to make sure that mode number is divisible by chunk dimensions
    if (isChunked == true)
    {
        if ((hdim[2] % chunk_dim[2] != 0) or (hdim[3] % chunk_dim[3] != 0) or (hdim[4] % chunk_dim[4] != 0))
        {
            HADRONS_ERROR(Size, "all-to-all chunk size mismatch (not a divisor of overall dimensions): got chunk dims"
            + std::to_string(chunk_dim[0]) + "x" + std::to_string(chunk_dim[1]) + "x"
            + std::to_string(chunk_dim[2]) + "x" + std::to_string(chunk_dim[3]) + "x"
            + std::to_string(chunk_dim[4]) + ", for overall dims "
            + std::to_string(hdim[0]) + "x"
            + std::to_string(hdim[1]) + "x" + std::to_string(hdim[2]) + "x"
            + std::to_string(hdim[3]) + "x" + std::to_string(hdim[4]));
        }
        else
        {
            chunkCount = hdim[2] / chunk_dim[2];
        }
    }

    A2AMatrixNuc<T>         buf(Ns, ni_, nj_, nk_);
 
    // Default values
    std::vector<hsize_t> count    = {1, static_cast<hsize_t>(Ns),
                                     static_cast<hsize_t>(ni_),
                                     static_cast<hsize_t>(nj_),
                                     static_cast<hsize_t>(nk_)},
                         stride   = {1, 1, 1, 1, 1},
                         block    = {1, 1, 1, 1, 1},
                         memCount = {static_cast<hsize_t>(Ns),
                                     static_cast<hsize_t>(ni_),
                                     static_cast<hsize_t>(nj_),
                                     static_cast<hsize_t>(nk_)};
    
    // Adjust hyperslab selection values if data is chunked
    if (isChunked == true)
    {
        count    = {1, 1,
                    static_cast<hsize_t>(chunkCount),
                    static_cast<hsize_t>(chunkCount),
                    static_cast<hsize_t>(chunkCount)},
        stride   = {1, 1, chunk_dim[2], chunk_dim[3], chunk_dim[4]},
        block    = {1, static_cast<hsize_t>(Ns), chunk_dim[2], chunk_dim[3], chunk_dim[4]},
        memCount = {static_cast<hsize_t>(Ns),
                    static_cast<hsize_t>(ni_),
                    static_cast<hsize_t>(nj_),
                    static_cast<hsize_t>(nk_)};
    }

    H5NS::DataSpace      memspace(memCount.size(), memCount.data());

    LOG(Message) << "Starting to load data." << std::endl;
    std::cout.flush();
    *tRead = 0.;
    for (unsigned int tp1 = nt_; tp1 > 0; --tp1)
    {
        unsigned int         t      = tp1 - 1;
        std::vector<hsize_t> offset = {static_cast<hsize_t>(t), 0, 0, 0, 0};
        
        // changed from t % 10 to t % 1 for debugging
        if (t % 1 == 0)
        {
            LOG(Message) << "Loading timeslice " << t << std::endl;
        }
        dataspace.selectHyperslab(H5S_SELECT_SET, count.data(), offset.data(),
                                  stride.data(), block.data());
        if (tRead) *tRead -= usecond();    
        dataset.read(buf.data(), datatype, memspace, dataspace);
        if (tRead) *tRead += usecond();
        v[t] = buf.template cast<VecT>();
    }
    std::cout << std::endl;
#else
    HADRONS_ERROR(Implementation, "all-to-all matrix I/O needs HDF5 library");
#endif
}

template <typename T>
template <template <class> class Vec, typename VecT>
void A2AMatrixNucIo<T>::load(Vec<VecT> &v, GridCartesian *grid, double *tRead)
{
#ifdef HAVE_HDF5
    Hdf5Reader           reader(filename_);
    std::vector<hsize_t> hdim;
    H5NS::DataSet        dataset;
    H5NS::DataSpace      dataspace;
    H5NS::CompType       datatype;
    H5NS::DSetCreatPropList plist;
    
    push(reader, dataname_);
    auto &group = reader.getGroup();
    dataset     = group.openDataSet(HADRONS_A2AN_NAME);
    plist       = dataset.getCreatePlist();
    datatype    = dataset.getCompType();
    dataspace   = dataset.getSpace();
    hdim.resize(dataspace.getSimpleExtentNdims());
    dataspace.getSimpleExtentDims(hdim.data());
        
    // chunking info
    bool isChunked = false;
    hsize_t chunk_dim[5];
    int rank_chunk;
    unsigned int chunkCount = 1;
    
    // adjust the chunk cache to try to improve performance
    hid_t dapl_id = H5Dget_access_plist(dataset.getId());
    hsize_t max_objects_in_chunk_cache = 12799;
    hsize_t max_bytes_in_cache = 128*1024*1024;
    
    H5Pset_chunk_cache(dapl_id, max_objects_in_chunk_cache, max_bytes_in_cache, H5D_CHUNK_CACHE_W0_DEFAULT);
    
    if (plist.getLayout() == H5D_CHUNKED)
    {
        isChunked = true;
        rank_chunk = plist.getChunk(5, chunk_dim);
        std::cout << "Data is chunked with rank " << std::to_string(rank_chunk) << std::endl;
        std::cout << "and chunk dimensions " << std::to_string(chunk_dim[0]) << " " <<
                                                std::to_string(chunk_dim[1]) << " " <<
                                                std::to_string(chunk_dim[2]) << " " <<
                                                std::to_string(chunk_dim[3]) << " " <<
                                                std::to_string(chunk_dim[4]) << std::endl;
    }
    
    std::cout << "DEBUG: Ns nt_ ni_ nj_ nk_ are " << std::to_string(Ns) << " " << std::to_string(nt_) << " " << std::to_string(ni_) << " "
              << std::to_string(nj_) << " " << std::to_string(nk_) << std::endl;
    std::cout << "hdims are " <<std::to_string(hdim[0]) << " " << std::to_string(hdim[1]) << " " << std::to_string(hdim[2]) << " "
              << std::to_string(hdim[3]) << " " << std::to_string(hdim[4]) << std::endl;
    if ((Ns*nt_*ni_*nj_*nk_ != 0) and
        ((hdim[0] != Ns) or (hdim[1] != nt_*grid->_processors[3]) or (hdim[2] != ni_) or (hdim[3] != nj_) or (hdim[4] != nk_)))
    {
        HADRONS_ERROR(Size, "all-to-all matrix size mismatch (got "
            + std::to_string(hdim[0]) + "x" + std::to_string(hdim[1]) + "x"
            + std::to_string(hdim[2]) + "x" + std::to_string(hdim[3]) + "x"
            + std::to_string(hdim[4]) + ", expected "
            + std::to_string(Ns) + "x"
            + std::to_string(nt_*grid->_processors[3]) + "x" + std::to_string(ni_) + "x"
            + std::to_string(nt_) + "x" + std::to_string(ni_) + "x"
            + std::to_string(nj_) + "x" + std::to_string(nk_));
    }
    else if (ni_*nj_*nk_ == 0)
    {
        if ((hdim[0] != nt_*grid->_processors[3]) or (hdim[1] != Ns))
        {
            HADRONS_ERROR(Size, "all-to-all time size mismatch (got "
                + std::to_string(hdim[0]) + "x"
                + std::to_string(hdim[1]) + ", expected "
                + std::to_string(nt_*grid->_processors[3]) + "x"
                + std::to_string(Ns) + ")");
        }
        ni_ = hdim[2];
        nj_ = hdim[3];
        nk_ = hdim[4];
    }
    
    // check to make sure that mode number is divisible by chunk dimensions
    if (isChunked == true)
    {
        if ((hdim[2] % chunk_dim[2] != 0) or (hdim[3] % chunk_dim[3] != 0) or (hdim[4] % chunk_dim[4] != 0))
        {
            HADRONS_ERROR(Size, "all-to-all chunk size mismatch (not a divisor of overall dimensions): got chunk dims"
            + std::to_string(chunk_dim[0]) + "x" + std::to_string(chunk_dim[1]) + "x"
            + std::to_string(chunk_dim[2]) + "x" + std::to_string(chunk_dim[3]) + "x"
            + std::to_string(chunk_dim[4]) + ", for overall dims "
            + std::to_string(hdim[0]) + "x"
            + std::to_string(hdim[1]) + "x" + std::to_string(hdim[2]) + "x"
            + std::to_string(hdim[3]) + "x" + std::to_string(hdim[4]));
        }
        else
        {
            chunkCount = hdim[2] / chunk_dim[2];
        }
    }

    A2AMatrixNuc<T>         buf(Ns, ni_, nj_, nk_);
 
    // Default values
    std::vector<hsize_t> count    = {1, static_cast<hsize_t>(Ns),
                                     static_cast<hsize_t>(ni_),
                                     static_cast<hsize_t>(nj_),
                                     static_cast<hsize_t>(nk_)},
                         stride   = {1, 1, 1, 1, 1},
                         block    = {1, 1, 1, 1, 1},
                         memCount = {static_cast<hsize_t>(Ns),
                                     static_cast<hsize_t>(ni_),
                                     static_cast<hsize_t>(nj_),
                                     static_cast<hsize_t>(nk_)};
    
    // Adjust hyperslab selection values if data is chunked
    if (isChunked == true)
    {
        count    = {1, 1,
                    static_cast<hsize_t>(chunkCount),
                    static_cast<hsize_t>(chunkCount),
                    static_cast<hsize_t>(chunkCount)},
        stride   = {1, 1, chunk_dim[2], chunk_dim[3], chunk_dim[4]},
        block    = {1, static_cast<hsize_t>(Ns), chunk_dim[2], chunk_dim[3], chunk_dim[4]},
        memCount = {static_cast<hsize_t>(Ns),
                    static_cast<hsize_t>(ni_),
                    static_cast<hsize_t>(nj_),
                    static_cast<hsize_t>(nk_)};
    }

    H5NS::DataSpace      memspace(memCount.size(), memCount.data());

    std::cout << "Loading timeslice";
    std::cout.flush();
    std::cout.flush();
    *tRead = 0.;
    int localNt = grid->LocalDimensions()[3];
    int globalNt = grid->GlobalDimensions()[3];
    int tshift = grid->_processor_coor[3] * localNt;
    for (unsigned int tp1 = localNt; tp1 > 0; --tp1)
    {
        unsigned int         t      = tp1 - 1;
        std::vector<hsize_t> offset = {static_cast<hsize_t>(t+tshift)%globalNt, 0, 0, 0, 0};
        
        if (t % 1 == 0)
        {
            std::cout << " " << t;
            std::cout.flush();
        }
        dataspace.selectHyperslab(H5S_SELECT_SET, count.data(), offset.data(),
                                  stride.data(), block.data());
        if (tRead) *tRead -= usecond();
        dataset.read(buf.data(), datatype, memspace, dataspace);
        if (tRead) *tRead += usecond();
        v[t] = buf.template cast<VecT>();
    }
    std::cout << std::endl;
#else
    HADRONS_ERROR(Implementation, "all-to-all matrix I/O needs HDF5 library");
#endif
}



/******************************************************************************
 *               A2AMatrixBlockComputation template implementation            *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename T, typename Field, typename MetadataType, typename TIo>
A2AMatrixNucleonBlockComputation<T, Field, MetadataType, TIo>
::A2AMatrixNucleonBlockComputation(GridBase *grid,
                            const unsigned int orthogDim,
                            const unsigned int next, 
                            const unsigned int nstr,
                            const unsigned int blockSize, 
                            const unsigned int cacheBlockSize,
                            TimerArray *tArray)
: grid_(grid), nt_(grid->GlobalDimensions()[orthogDim]), orthogDim_(orthogDim)
, next_(next), nstr_(nstr), blockSize_(blockSize), cacheBlockSize_(cacheBlockSize)
, tArray_(tArray)
{
    mCache_.resize(nt_*next_*nstr_*cacheBlockSize_*cacheBlockSize_*cacheBlockSize_);
    mBuf_.resize(nt_*next_*nstr_*blockSize_*blockSize_*blockSize_);
}

#define START_TIMER(name) if (tArray_) tArray_->startTimer(name)
#define STOP_TIMER(name)  if (tArray_) tArray_->stopTimer(name)
#define GET_TIMER(name)   ((tArray_ != nullptr) ? tArray_->getDTimer(name) : 0.)

// execution ///////////////////////////////////////////////////////////////////
template <typename T, typename Field, typename MetadataType, typename TIo>
void A2AMatrixNucleonBlockComputation<T, Field, MetadataType, TIo>
::execute(const std::vector<Field> &left, const std::vector<Field> &right, const std::vector<Field> &q3,
          A2AKernelNucleon<T, Field> &kernel, const FilenameFn &ionameFn,
          const FilenameFn &filenameFn, const MetadataFn &metadataFn)
{
    //////////////////////////////////////////////////////////////////////////
    // i,j,k   is first  loop over blockSize_ factors
    // ii,jj,kk is second loop over cacheBlockSize_ factors for high perf contractions
    // iii,jjj,kkk are loops within cacheBlock
    // Total index is sum of these  i+ii+iii etc...
    //////////////////////////////////////////////////////////////////////////
    int    N_i = left.size();
    int    N_j = right.size();
    int    N_k = q3.size();
    double flops, bytes, t_kernel;
    double nodes = grid_->NodeCount();
    
    int NBlock_i = N_i/blockSize_ + (((N_i % blockSize_) != 0) ? 1 : 0);
    int NBlock_j = N_j/blockSize_ + (((N_j % blockSize_) != 0) ? 1 : 0);
    int NBlock_k = N_k/blockSize_ + (((N_k % blockSize_) != 0) ? 1 : 0);

    for(int i=0;i<N_i;i+=blockSize_)
    for(int j=0;j<N_j;j+=blockSize_)
    for(int k=0;k<N_k;k+=blockSize_)
    {
        // Get the W and V vectors for this block^3 set of terms
        int N_ii = MIN(N_i-i,blockSize_);
        int N_jj = MIN(N_j-j,blockSize_);
        int N_kk = MIN(N_k-k,blockSize_);
        // MCA - need to change this function to tensor from matrix???
        A2AMatrixSetNuc<TIo> mBlock(mBuf_.data(), next_, nt_, nstr_, N_ii, N_jj, N_kk);

        LOG(Message) << "All-to-all matrix block " 
                     << k/blockSize_ + NBlock_k*j/blockSize_ + NBlock_k*NBlock_j*i/blockSize_ + 1 
                     << "/" << NBlock_i*NBlock_j*NBlock_k << " [" << i <<" .. " 
                     << i+N_ii-1 << ", " << j <<" .. " << j+N_jj-1 << ", "
                     << k << " .. " << k+N_kk-1 <<"]" 
                     << std::endl;
        // Series of cache blocked chunks of the contractions within this block
        flops    = 0.0;
        bytes    = 0.0;
        t_kernel = 0.0;
        for(int ii=0;ii<N_ii;ii+=cacheBlockSize_)
        for(int jj=0;jj<N_jj;jj+=cacheBlockSize_)
        for(int kk=0;kk<N_kk;kk+=cacheBlockSize_)
        {
            double t;
            int N_iii = MIN(N_ii-ii,cacheBlockSize_);
            int N_jjj = MIN(N_jj-jj,cacheBlockSize_);
            int N_kkk = MIN(N_kk-kk,cacheBlockSize_);
            // MCA - need to change this function too
            A2AMatrixSetNuc<T> mCacheBlock(mCache_.data(), next_, nstr_, nt_, N_iii, N_jjj, N_kkk);

            START_TIMER("kernel");
            kernel(mCacheBlock, &left[i+ii], &right[j+jj], &q3[k+kk], orthogDim_, t);
            STOP_TIMER("kernel");
            t_kernel += t;
            flops    += kernel.flops(N_iii, N_jjj, N_kkk);
            bytes    += kernel.bytes(N_iii, N_jjj, N_kkk);

            START_TIMER("cache copy");
            thread_for_collapse(6, e, next_, {
                for(int s=0;s< nstr_;s++)
                    for(int t=0;t< nt_;t++)
                        for(int iii=0;iii< N_iii;iii++)
                            for(int jjj=0;jjj< N_jjj;jjj++)
                                for(int kkk=0;kkk< N_kkk;kkk++){
                                    mBlock(e,t,s,ii+iii,jj+jjj,kk+kkk) = mCacheBlock(e,s,t,iii,jjj,kkk);
                                }
            });
            STOP_TIMER("cache copy");
        }

	LOG(Message) << "t | i | j | k | mu Value" << std::endl;
        // perf
        LOG(Message) << "Kernel perf " << flops/t_kernel/1.0e3/nodes 
                     << " Gflop/s/node " << std::endl;
        LOG(Message) << "Kernel perf " << bytes/t_kernel*1.0e6/1024/1024/1024/nodes 
                     << " GB/s/node "  << std::endl;

        // IO
        double       blockSize, ioTime;
        unsigned int myRank = grid_->ThisRank(), nRank  = grid_->RankCount();
    
        LOG(Message) << "Writing block to disk" << std::endl;
        ioTime = -GET_TIMER("IO: write block");
        START_TIMER("IO: total");
        makeFileDir(filenameFn(0), grid_);
#ifdef HADRONS_A2AN_PARALLEL_IO
        grid_->Barrier();
        // make task list for current node
        nodeIo_.clear();
        for(int f = myRank; f < next_; f += nRank)
        {
            IoHelper h;

            h.i  = i;
            h.j  = j;
            h.k  = k;
            h.e  = f/nstr_;
            h.s  = 0;
            h.io = A2AMatrixNucIo<TIo>(filenameFn(h.e), 
                                    ionameFn(h.e), nt_, N_i, N_j, N_k);
            h.md = metadataFn(h.e);
            nodeIo_.push_back(h);
        }
        std::cout << "Still working" << std::endl;
        // parallel IO
        for (auto &h: nodeIo_)
        {
            saveBlock(mBlock, h);
        }
        grid_->Barrier();
        std::cout << "Still working after saveBlock" << std::endl;
#else
        // serial IO, for testing purposes only
        for(int e = 0; e < next_; e++)
        //for(int s = 0; s < nstr_; s++)
        {
            IoHelper h;

            h.i  = i;
            h.j  = j;
            h.k  = k;
            h.e  = e;
            h.s  = 0;
            h.io = A2AMatrixNucIo<TIo>(filenameFn(h.e), 
                                    ionameFn(h.e), nt_, N_i, N_j, N_k);
            h.md = metadataFn(h.e);
            saveBlock(mBlock, h);
        }
#endif
        STOP_TIMER("IO: total");
        blockSize  = static_cast<double>(next_*nt_*nstr_*N_ii*N_jj*N_kk*sizeof(TIo));
        ioTime    += GET_TIMER("IO: write block");
        LOG(Message) << "HDF5 IO done " << sizeString(blockSize) << " in "
                     << ioTime  << " us (" 
                     << blockSize/ioTime*1.0e6/1024/1024
                     << " MB/s)" << std::endl;
    }
}

// I/O handler /////////////////////////////////////////////////////////////////
template <typename T, typename Field, typename MetadataType, typename TIo>
void A2AMatrixNucleonBlockComputation<T, Field, MetadataType, TIo>
::saveBlock(const A2AMatrixSetNuc<TIo> &m, IoHelper &h)
{
    if ((h.i == 0) and (h.j == 0) and (h.k == 0))
    {
        START_TIMER("IO: file creation");
        h.io.initFile(h.md, blockSize_);
        STOP_TIMER("IO: file creation");
    }
    START_TIMER("IO: write block");
    h.io.saveBlock(m, h.e, h.s, h.i, h.j, h.k);
    STOP_TIMER("IO: write block");
}

#undef START_TIMER
#undef STOP_TIMER
#undef GET_TIMER

END_HADRONS_NAMESPACE

#endif // A2A_Matrix_hpp_
