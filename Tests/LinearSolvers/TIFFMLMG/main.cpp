
#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>

#include <prob_par.H>

using namespace amrex;

void init_prob_parms ();
void init_prob (const Vector<Geometry>& geom, Vector<MultiFab>& alpha, Vector<MultiFab>& beta,
                Vector<MultiFab>& rhs, Vector<MultiFab>& exact);
void solve_with_mlmg (const Vector<Geometry>& geom, int rr,
                      Vector<MultiFab>& soln,
                      const Vector<MultiFab>& alpha, const Vector<MultiFab>& beta,
                      Vector<MultiFab>& rhs, const Vector<MultiFab>& exact);
void write_plotfile (const Vector<Geometry>& geom, int rr,
                     const Vector<MultiFab>& soln, const Vector<MultiFab>& exact,
                     const Vector<MultiFab>& alpha, const Vector<MultiFab>& beta,
                     const Vector<MultiFab>& rhs);

namespace {
    void build_geometry_and_grids (Vector<Geometry>& geom, Vector<BoxArray>& grids);
    BoxArray readBoxList (const std::string& file, Box& domain);
    void init_from_tiff2D (const std::string& file, 
        const Vector<Geometry>& vecgeom, Vector<MultiFab>& vecmf);
}

namespace {
    int max_level     = 1;
    int nlevels       = 2;
    int n_cell        = 64;
    int max_grid_size = 32;
    int ref_ratio     = 2;
    std::string boxes_file;
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    
    {
        BL_PROFILE("main()");

        init_prob_parms();

        Vector<Geometry> geom;
        Vector<BoxArray> grids;
        build_geometry_and_grids(geom, grids);
        
        Vector<MultiFab> soln(nlevels);
        Vector<MultiFab> exact(nlevels);
        Vector<MultiFab> alpha(nlevels);
        Vector<MultiFab> beta(nlevels);
        Vector<MultiFab> rhs(nlevels);
        
        for (int ilev = 0; ilev < nlevels; ++ilev)
        {
            DistributionMapping dm{grids[ilev]};
            // 2 ghost cell to store boundary conditions for HO operator
            soln [ilev].define(grids[ilev], dm, 1, 2);  
            // 1 ghost cell for 4th-order averaging
            exact[ilev].define(grids[ilev], dm, 1, 1);
            // 1 ghost cell for 4th-order averaging
            alpha[ilev].define(grids[ilev], dm, 1, 1);
            // 1 ghost cell for averaging to faces
            beta [ilev].define(grids[ilev], dm, 1, 1);  
            // 1 ghost cell for 4th-order averaging
            rhs  [ilev].define(grids[ilev], dm, 1, 1);
        }

        init_prob (geom, alpha, beta, rhs, exact);

        for (auto& mf : soln) {
            mf.setVal(0.0); // initial guess
        }
        init_from_tiff2D("Particle_based0000.tif", geom, soln);

        // solve_with_mlmg (geom, ref_ratio, soln, alpha, beta, rhs, exact);

        write_plotfile (geom, ref_ratio, soln, exact, alpha, beta, rhs);
    }

    amrex::Finalize();
}

namespace {

void build_geometry_and_grids (Vector<Geometry>& geom, Vector<BoxArray>& grids)
{
    ParmParse pp;
    pp.query("n_cell", n_cell);
    pp.query("max_level", max_level);
    pp.query("max_grid_size", max_grid_size);
    pp.query("ref_ratio", ref_ratio);
    pp.query("boxes", boxes_file);

    if (!boxes_file.empty())
    {
        Box dmn;
        const BoxArray& ba = readBoxList(boxes_file, dmn);

        const BoxArray& uncovered = amrex::complementIn(dmn, ba);

        if (uncovered.empty() && dmn.isSquare() && dmn.smallEnd() == IntVect::TheZeroVector())
        {
            max_level = 0;
            nlevels = max_level + 1;
            n_cell = dmn.longside();

            geom.resize(nlevels);
            grids.resize(nlevels);

            grids[0] = ba;
        }
        else
        {
            amrex::Print() << "Add a coarse level\n";
            max_level = 1;
            nlevels = max_level + 1;

            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(dmn.coarsenable(ref_ratio), "Domain must be coarsenable");
            dmn.coarsen(ref_ratio);
            dmn.setSmall(IntVect::TheZeroVector());
            n_cell = dmn.longside();

            geom.resize(nlevels);
            grids.resize(nlevels);

            IntVect dom0_lo {IntVect::TheZeroVector()};
            IntVect dom0_hi {AMREX_D_DECL(n_cell-1, n_cell-1, n_cell-1)};
            Box dom0 {dom0_lo, dom0_hi};
            grids[0] = BoxArray{dom0};
            grids[0].maxSize(max_grid_size);

            grids[1] = ba;
        }
    }
    else
    {
        nlevels = max_level + 1;
        geom.resize(nlevels);
        grids.resize(nlevels);
        
        IntVect dom0_lo {IntVect::TheZeroVector()};
        IntVect dom0_hi {AMREX_D_DECL(n_cell-1, n_cell-1, n_cell-1)};
        
        Box dom0 {dom0_lo, dom0_hi};
        BoxArray ba0{dom0};
        
        grids[0] = ba0;
        grids[0].maxSize(max_grid_size);
        
        for (int ilev=1, n=grids.size(); ilev < n; ++ilev)
        {
            ba0.grow(-n_cell/4);
            ba0.refine(ref_ratio);
            grids[ilev] = ba0;
            grids[ilev].maxSize(max_grid_size);
        }
    }
    
    std::array<Real,AMREX_SPACEDIM> prob_lo{AMREX_D_DECL(0.,0.,0.)};
    std::array<Real,AMREX_SPACEDIM> prob_hi{AMREX_D_DECL(1.,1.,1.)};
    RealBox real_box{prob_lo, prob_hi};
    
    const int coord = 0;  // Cartesian coordinates
    std::array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0,0,0)};
    if (prob::bc_type == MLLinOp::BCType::Periodic)
    {
        std::fill(is_periodic.begin(), is_periodic.end(), 1);
    }

    IntVect dom0_lo {IntVect::TheZeroVector()};
    IntVect dom0_hi {AMREX_D_DECL(n_cell-1, n_cell-1, n_cell-1)};
    Box dom0 {dom0_lo, dom0_hi};
    
    geom[0].define(dom0, &real_box, coord, is_periodic.data());
    for (int ilev=1, n=grids.size(); ilev < n; ++ilev)
    {
        dom0.refine(ref_ratio);
        geom[ilev].define(dom0, &real_box, coord, is_periodic.data());
    }
}

BoxArray
readBoxList (const std::string& file, Box& domain)
{
    BoxList retval;

    std::ifstream boxspec;

    boxspec.open(file.c_str(), std::ios::in);

    if( !boxspec )
    {
        std::string msg = "readBoxList: unable to open ";
        msg += file;
        amrex::Error(msg.c_str());
    }
    boxspec >> domain;
    
    int numbox = 0;
    boxspec >> numbox;

    for ( int i=0; i<numbox; i++ )
    {
        Box tmpbox;
        boxspec >> tmpbox;
        if( !domain.contains(tmpbox) )
	{
            std::cerr << "readBoxList: bogus box " << tmpbox << '\n';
            exit(1);
        }
        retval.push_back(tmpbox);
    }

    return BoxArray(retval);
}

#include <tiffio.h>

void
init_from_tiff2D (const std::string& file, const Vector<Geometry>& vecgeom, 
    Vector<MultiFab>& vecmf)
{
  std::cout << "TIFF file " << file << '\n';
  TIFF *tif=TIFFOpen(file.c_str(), "r");
  uint32 width, height;
  TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width);
  TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &height);
  std::cout << "Size: " << width << "w x " << height << "h\n";
  uint32 npixels = width*height;
  uint32* raster=(uint32 *) _TIFFmalloc(npixels *sizeof(uint32));
  int ok = TIFFReadRGBAImage(tif, width, height, raster, 0);
  std::cout << "Ok read?: " << ok << "\n";
  
  int refRatio = vecgeom[1].Domain().size()[0]/vecgeom[0].Domain().size()[0];
  std::cout << "Ref ratio: " << refRatio << "\n";
  const int levels = vecgeom.size();
  for (int ilev = levels-1; ilev >= 0; --ilev)
  {
    std::cout << "Level " << ilev << "\n";
    const Geometry& geom = vecgeom[ilev];
    MultiFab& mf = vecmf[ilev];
    Box domainBox = geom.Domain();
    IntVect domainSize = domainBox.size();
    IntVect origin = domainBox.smallEnd();
    std::cout << "Domain size " << domainSize[0] << "w x " << domainSize[1] << "h\n";
    std::cout << "MultiFab size " << mf.size() << "\n";

    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.validbox();
      IntVect rasterOffset = bx.smallEnd() - origin;
      rasterOffset += IntVect(20,50); // FIXME - just shifting to middle
      const Box fabBx = mf[mfi].box();
      IntVect offset = (fabBx.size() - bx.size())/2; // ghost cells
      int jstride = fabBx.size()[0];
      IntVect size = bx.size();
      for (int j=0; j < size[1]; j++)
      {
        for (int i=0; i < size[0]; i++)
        {
          int ix = i + rasterOffset[0] + (j + rasterOffset[1])* width;
          // ix = floor(ix/(refRatio);
          double* dataPtr = mf[mfi].dataPtr();

          int ixdata = i + offset[0] + (j + offset[1])* jstride;
          // if (i >= 20) break;
          int val = raster[ix];
          if (val == -1)
            dataPtr[ixdata] = 1;
          else
            dataPtr[ixdata] = .1;

          // std::cout << "(" << i << "," << j << ") = " << val << "\n";
        }
        // if (j >= 20) break;
      }
    }
  }

  _TIFFfree(raster);
}

}

