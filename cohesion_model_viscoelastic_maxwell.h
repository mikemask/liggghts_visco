/* ----------------------------------------------------------------------
    This is the

    ██╗     ██╗ ██████╗  ██████╗  ██████╗ ██╗  ██╗████████╗███████╗
    ██║     ██║██╔════╝ ██╔════╝ ██╔════╝ ██║  ██║╚══██╔══╝██╔════╝
    ██║     ██║██║  ███╗██║  ███╗██║  ███╗███████║   ██║   ███████╗
    ██║     ██║██║   ██║██║   ██║██║   ██║██╔══██║   ██║   ╚════██║
    ███████╗██║╚██████╔╝╚██████╔╝╚██████╔╝██║  ██║   ██║   ███████║
    ╚══════╝╚═╝ ╚═════╝  ╚═════╝  ╚═════╝ ╚═╝  ╚═╝   ╚═╝   ╚══════╝®

    DEM simulation engine, released by
    DCS Computing Gmbh, Linz, Austria
    http://www.dcs-computing.com, office@dcs-computing.com

    LIGGGHTS® is part of CFDEM®project:
    http://www.liggghts.com | http://www.cfdem.com

    Core developer and main author:
    Christoph Kloss, christoph.kloss@dcs-computing.com

    LIGGGHTS® is open-source, distributed under the terms of the GNU Public
    License, version 2 or later. It is distributed in the hope that it will
    be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. You should have
    received a copy of the GNU General Public License along with LIGGGHTS®.
    If not, see http://www.gnu.org/licenses . See also top-level README
    and LICENSE files.

    LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
    the producer of the LIGGGHTS® software and the CFDEM®coupling software
    See http://www.cfdem.com/terms-trademark-policy for details.

-------------------------------------------------------------------------
    Contributing author and copyright for this file:

    Michael Mascara (DCS Computing GmbH, TU Graz)
------------------------------------------------------------------------- */

#ifdef COHESION_MODEL
COHESION_MODEL(VISCOELASTIC_MAXWELL,maxwell,5)
#else
#ifndef VISCOELASTIC_MAXWELL_H_
#define VISCOELASTIC_MAXWELL_H_
#include "contact_models.h"
#include "cohesion_model_base.h"
#include <cmath>
#include "math_extra_liggghts.h"
#include "update.h"
#include "global_properties.h"
#include "fix_property_atom.h"
#include <iomanip>
#include <iostream>

//namespace MODEL_PARAMS
//{
    //inline static ScalarProperty* createTimestepCreateBond(PropertyRegistry & registry, const char * caller, bool sanity_checks)
    //{
    //  ScalarProperty* timestepCreateBondScalar = MODEL_PARAMS::createScalarProperty(registry, "timestepCreateBond", caller);
    //  return timestepCreateBondScalar;
    //}

    //inline static ScalarProperty* createRadiusMultiplier(PropertyRegistry & registry, const char * caller, bool sanity_checks)
    //{
    //  ScalarProperty* radiusMultiplierScalar = MODEL_PARAMS::createScalarProperty(registry, "radiusMultiplier", caller);
    //  return radiusMultiplierScalar;
    //}

    //inline static ScalarProperty* createMaxDistance(PropertyRegistry & registry, const char * caller, bool sanity_checks)
    //{
    //  ScalarProperty* maxDistanceScalar = MODEL_PARAMS::createScalarProperty(registry, "maxDistance", caller);
    //  return maxDistanceScalar;
    //}

    //inline static ScalarProperty* createYMax1(PropertyRegistry & registry, const char * caller, bool sanity_checks)
    //{
    //  ScalarProperty* yMax1Scalar = MODEL_PARAMS::createScalarProperty(registry, "yMax1", caller);
    //  return yMax1Scalar;
    //}

    //inline static ScalarProperty* createYMax2(PropertyRegistry & registry, const char * caller, bool sanity_checks)
    //{
    //  ScalarProperty* yMax2Scalar = MODEL_PARAMS::createScalarProperty(registry, "yMax2", caller);
    //  return yMax2Scalar;
    //}

    //inline static ScalarProperty* createYMax3(PropertyRegistry & registry, const char * caller, bool sanity_checks)
    //{
    //  ScalarProperty* yMax3Scalar = MODEL_PARAMS::createScalarProperty(registry, "yMax3", caller);
    //  return yMax3Scalar;
    //}

    //inline static ScalarProperty* createYMax4(PropertyRegistry & registry, const char * caller, bool sanity_checks)
    //{
    //  ScalarProperty* yMax4Scalar = MODEL_PARAMS::createScalarProperty(registry, "yMax4", caller);
    //  return yMax4Scalar;
    //}

    //inline static ScalarProperty* createYSpring(PropertyRegistry & registry, const char * caller, bool sanity_checks)
    //{
    //  ScalarProperty* ySpringScalar = MODEL_PARAMS::createScalarProperty(registry, "ySpring1", caller);
    //  return ySpringScalar;
    //}

    //inline static ScalarProperty* createMuMax1(PropertyRegistry & registry, const char * caller, bool sanity_checks)
    //{
    //  ScalarProperty* muMax1Scalar = MODEL_PARAMS::createScalarProperty(registry, "muMax1", caller);
    //  return muMax1Scalar;
    //}

    //inline static ScalarProperty* createMuMax2(PropertyRegistry & registry, const char * caller, bool sanity_checks)
    //{
    //  ScalarProperty* muMax2Scalar = MODEL_PARAMS::createScalarProperty(registry, "muMax2", caller);
    //  return muMax2Scalar;
    //}

    //inline static ScalarProperty* createMuMax3(PropertyRegistry & registry, const char * caller, bool sanity_checks)
    //{
    //  ScalarProperty* muMax3Scalar = MODEL_PARAMS::createScalarProperty(registry, "muMax3", caller);
    //  return muMax3Scalar;
    //}

    //inline static ScalarProperty* createMuMax4(PropertyRegistry & registry, const char * caller, bool sanity_checks)
    //{
    //  ScalarProperty* muMax4Scalar = MODEL_PARAMS::createScalarProperty(registry, "muMax4", caller);
    //  return muMax4Scalar;
    //}

//}

namespace LIGGGHTS {
namespace ContactModels
{
  template<>
  class CohesionModel<VISCOELASTIC_MAXWELL> : public CohesionModelBase
  {
    public:
      CohesionModel(LAMMPS * lmp, IContactHistorySetup * hsetup,class ContactModelBase *cmb) :
        CohesionModelBase(lmp, hsetup, cmb)
        //timestepCreateBond(-1),
        //radiusMultiplier(1.),
        //maxDistance(0.),
        //yMax1(0.),
        //yMax2(0.),
        //yMax3(0.),
        //yMax4(0.),
        //ySpring(0.),
        //muMax1(0.),
        //muMax2(0.),
        //muMax3(0.),
        //muMax4(0.)
    {
      history = hsetup->add_history_value("bond_exists", "0");
      hsetup->add_history_value("r_0x", "1");
      hsetup->add_history_value("r_0y", "1");
      hsetup->add_history_value("r_0z", "1");
      hsetup->add_history_value("r_oldx", "1");
      hsetup->add_history_value("r_oldy", "1");
      hsetup->add_history_value("r_oldz", "1");
      hsetup->add_history_value("f_n1x", "1");
      hsetup->add_history_value("f_n1y", "1");
      hsetup->add_history_value("f_n1z", "1");
      hsetup->add_history_value("f_n2x", "1");
      hsetup->add_history_value("f_n2y", "1");
      hsetup->add_history_value("f_n2z", "1");
      hsetup->add_history_value("f_n3x", "1");
      hsetup->add_history_value("f_n3y", "1");
      hsetup->add_history_value("f_n3z", "1");
      hsetup->add_history_value("f_n4x", "1");
      hsetup->add_history_value("f_n4y", "1");
      hsetup->add_history_value("f_n4z", "1");
      hsetup->add_history_value("f_s1x", "1");
      hsetup->add_history_value("f_s1y", "1");
      hsetup->add_history_value("f_s1z", "1");
      hsetup->add_history_value("f_s2x", "1");
      hsetup->add_history_value("f_s2y", "1");
      hsetup->add_history_value("f_s2z", "1");
      hsetup->add_history_value("f_s3x", "1");
      hsetup->add_history_value("f_s3y", "1");
      hsetup->add_history_value("f_s3z", "1");
      hsetup->add_history_value("f_s4x", "1");
      hsetup->add_history_value("f_s4y", "1");
      hsetup->add_history_value("f_s4z", "1");
    }

    void registerSettings(Settings & settings)
    {
      settings.registerOnOff("hertz_active", hertz_on_off, false);
    }

    inline void postSettings(IContactHistorySetup * hsetup, ContactModelBase *cmb) {}

    void connectToProperties(PropertyRegistry & registry)
    {

      //registry.registerProperty("timestepCreateBond", &MODEL_PARAMS::createTimestepCreateBond);
      //registry.registerProperty("radiusMultiplier", &MODEL_PARAMS::createRadiusMultiplier);
      //registry.registerProperty("maxDistance", &MODEL_PARAMS::createMaxDistance);
      //registry.registerProperty("yMax1", &MODEL_PARAMS::createYMax1);
      //registry.registerProperty("yMax2", &MODEL_PARAMS::createYMax2);
      //registry.registerProperty("yMax3", &MODEL_PARAMS::createYMax3);
      //registry.registerProperty("yMax4", &MODEL_PARAMS::createYMax4);
      //registry.registerProperty("ySpring", &MODEL_PARAMS::createYSpring);
      //registry.registerProperty("muMax1", &MODEL_PARAMS::createMuMax1);
      //registry.registerProperty("muMax2", &MODEL_PARAMS::createMuMax2);
      //registry.registerProperty("muMax3", &MODEL_PARAMS::createMuMax3);
      //registry.registerProperty("muMax4", &MODEL_PARAMS::createMuMax4);

      //registry.connect("timestepCreateBond", timestepCreateBond, "cohesion_model maxwell");
      //registry.connect("radiusMultiplier", radiusMultiplier, "cohesion_model maxwell");
      //registry.connect("maxDistance", maxDistance, "cohesion_model maxwell");
      //registry.connect("yMax1", yMax1, "cohesion_model maxwell");
      //registry.connect("yMax2", yMax2, "cohesion_model maxwell");
      //registry.connect("yMax3", yMax3, "cohesion_model maxwell");
      //registry.connect("yMax4", yMax4, "cohesion_model maxwell");
      //registry.connect("ySpring", ySpring, "cohesion_model maxwell");
      //registry.connect("muMax1", muMax1, "cohesion_model maxwell");
      //registry.connect("muMax2", muMax2, "cohesion_model maxwell");
      //registry.connect("muMax3", muMax3, "cohesion_model maxwell");
      //registry.connect("muMax4", muMax4, "cohesion_model maxwell");

      // error checks on coarsegraining
      if(force->cg_active())
        error->cg(FLERR,"cohesion model maxwell");

      class FixPropertyGlobal* max_dist=static_cast<FixPropertyGlobal*>(modify->find_fix_property("maxDistance","property/global","scalar",0,0,force->pair_style));

      const double max_rad = registry.max_radius();
      const double max_dist_ratio = 0.5*max_dist->compute_scalar()/max_rad;

      neighbor->register_contact_dist_factor(1.1*max_dist_ratio);
    }

    inline void endSurfacesIntersect(SurfacesIntersectData &sidata, ForceData&, ForceData&) {}
    void beginPass(SurfacesIntersectData&, ForceData&, ForceData&){}
    void endPass(SurfacesIntersectData&, ForceData&, ForceData&){}

    void surfacesIntersect(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces)
    {
      if (hertz_on_off==true && sidata.is_wall==false)
      {

        fprintf(screen, "hertz!\n");

        int max_type = atom->get_properties()->max_type();

        const int itype = sidata.itype;
        const int jtype = sidata.jtype;

        const double ri = sidata.radi;
        const double rj = sidata.radj;
        const double reff = ri*rj/(ri+rj);
        const double meff=sidata.meff;

        const double *Y=static_cast<FixPropertyGlobal*>(modify->find_fix_property("youngsModulus","property/global","peratomtype",max_type,0,force->pair_style))->get_values();
        const double *v=static_cast<FixPropertyGlobal*>(modify->find_fix_property("poissonsRatio","property/global","peratomtype",max_type,0,force->pair_style))->get_values();
        //class FixPropertyGlobal* cr=static_cast<FixPropertyGlobal*>(modify->find_fix_property("coefficientRestitution","property/global","peratomtypepair",0,0,force->pair_style));

        double Yeff;

        //pre-calculate parameters for possible contact material combinations
        for(int i=1; i<max_type+1; i++)
        {
          for(int j=1; j<max_type+1; j++)
          {
            Yeff = 1./((1.-pow(v[i-1],2.))/Y[i-1]+(1.-pow(v[j-1],2.))/Y[j-1]);
            //double betaEff = log(cri)/sqrt(log(cri)*log(cri)+M_PI*M_PI);
          }
        }


        const double r = sidata.r;
        const double deltan = ri+rj-r;

        //fprintf(screen, "dn = %f \n", deltan);
        //fprintf(screen, "dn = %f \n", r);

        const double sqrtval = sqrt(reff*deltan);

        const double Sn = 2.*Yeff*sqrtval;
        double kn = 4./3.*Yeff*sqrtval;

        const double sqrt_five_six = 0.91287092917527685576161630466800355658790782499663875;

        //const double gamman = -2.*sqrt_five_six*betaEff*sqrt(Sn*meff);

        kn /= force->nktv2p;

        //const double Fn_damping = -gamman*sidata.vn;
        const double Fn_contact = kn*deltan;
        double Fn = Fn_contact;
        //double Fn = Fn_damping + Fn_contact;

        if (Fn<0.0) Fn = 0.0;

        sidata.Fn = Fn;

        // apply normal force
        i_forces.delta_F[0] += sidata.Fn * sidata.en[0];
        i_forces.delta_F[1] += sidata.Fn * sidata.en[1];
        i_forces.delta_F[2] += sidata.Fn * sidata.en[2];

        j_forces.delta_F[0] += -i_forces.delta_F[0];
        j_forces.delta_F[1] += -i_forces.delta_F[1];
        j_forces.delta_F[2] += -i_forces.delta_F[2];
      }

      const double dist = sidata.r;

      class FixPropertyGlobal* max_dist=static_cast<FixPropertyGlobal*>(modify->find_fix_property("maxDistance","property/global","scalar",0,0,force->pair_style));

      if (dist >= max_dist->compute_scalar() || sidata.is_wall)
          return;

      double * hist_arr = &sidata.contact_history[history];

      class FixPropertyGlobal* tsCreate=static_cast<FixPropertyGlobal*>(modify->find_fix_property("timestepCreateBond","property/global","scalar",0,0,force->pair_style));
      class FixPropertyGlobal* radMult=static_cast<FixPropertyGlobal*>(modify->find_fix_property("radiusMultiplier","property/global","scalar",0,0,force->pair_style));

      const int ts_create = tsCreate->compute_scalar();
      const double lambda = radMult->compute_scalar();

      const double dt = update->dt;

      const int itype = sidata.itype;
      const int jtype = sidata.jtype;

      const double ri = sidata.radi;
      const double rj = sidata.radj;

      const double r[3] = {sidata.delta[0], sidata.delta[1], sidata.delta[2]};
      const double r_0[3] = {hist_arr[1],hist_arr[2],hist_arr[3]};
      double r_old[3] = {hist_arr[4],hist_arr[5],hist_arr[6]};

      const double d_0 = sqrt(r_0[0]*r_0[0]+r_0[1]*r_0[1]+r_0[2]*r_0[2]);
      const double d_old = sqrt(r_old[0]*r_old[0]+r_old[1]*r_old[1]+r_old[2]*r_old[2]);

      const double n[3] = {r[0]/dist,r[1]/dist,r[2]/dist};

      if (hist_arr[0] < 0.5)
      {
          const bool create_now = (ts_create < 0 || update->ntimestep == ts_create);
          const double dist_create = lambda*(ri+rj);
          const bool is_close = dist_create >= dist;

          if (create_now && is_close)
          {
              createBond(sidata, r);
              return;
          }
          else
              return;
      }

      double f_n1[3] = {hist_arr[7],hist_arr[8],hist_arr[9]};
      double f_n2[3] = {hist_arr[10],hist_arr[11],hist_arr[12]};
      double f_n3[3] = {hist_arr[13],hist_arr[14],hist_arr[15]};
      double f_n4[3] = {hist_arr[16],hist_arr[17],hist_arr[18]};

      double f_s1[3] = {hist_arr[19],hist_arr[20],hist_arr[21]};
      double f_s2[3] = {hist_arr[22],hist_arr[23],hist_arr[24]};
      double f_s3[3] = {hist_arr[25],hist_arr[26],hist_arr[27]};
      double f_s4[3] = {hist_arr[28],hist_arr[29],hist_arr[30]};

      const double rb = std::min(ri,rj);
      const double Ab = M_PI*rb*rb;

      class FixPropertyGlobal* K_1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("yMax1","property/global","scalar",0,0,force->pair_style));
      class FixPropertyGlobal* K_2=static_cast<FixPropertyGlobal*>(modify->find_fix_property("yMax2","property/global","scalar",0,0,force->pair_style));
      class FixPropertyGlobal* K_3=static_cast<FixPropertyGlobal*>(modify->find_fix_property("yMax3","property/global","scalar",0,0,force->pair_style));
      class FixPropertyGlobal* K_4=static_cast<FixPropertyGlobal*>(modify->find_fix_property("yMax4","property/global","scalar",0,0,force->pair_style));
      class FixPropertyGlobal* K_k=static_cast<FixPropertyGlobal*>(modify->find_fix_property("ySpring","property/global","scalar",0,0,force->pair_style));
      class FixPropertyGlobal* C_1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("muMax1","property/global","scalar",0,0,force->pair_style));
      class FixPropertyGlobal* C_2=static_cast<FixPropertyGlobal*>(modify->find_fix_property("muMax2","property/global","scalar",0,0,force->pair_style));
      class FixPropertyGlobal* C_3=static_cast<FixPropertyGlobal*>(modify->find_fix_property("muMax3","property/global","scalar",0,0,force->pair_style));
      class FixPropertyGlobal* C_4=static_cast<FixPropertyGlobal*>(modify->find_fix_property("muMax4","property/global","scalar",0,0,force->pair_style));

      const double k1_n = Ab/d_0*K_1->compute_scalar();
      const double k2_n = Ab/d_0*K_2->compute_scalar();
      const double k3_n = Ab/d_0*K_3->compute_scalar();
      const double k4_n = Ab/d_0*K_4->compute_scalar();
      const double kk_n = Ab/d_0*K_k->compute_scalar();
      const double c1_n = Ab/d_0*C_1->compute_scalar();
      const double c2_n = Ab/d_0*C_2->compute_scalar();
      const double c3_n = Ab/d_0*C_3->compute_scalar();
      const double c4_n = Ab/d_0*C_4->compute_scalar();

      const double k1_s = k1_n;
      const double k2_s = k2_n;
      const double k3_s = k3_n;
      const double k4_s = k4_n;
      const double kk_s = kk_n;
      const double c1_s = c1_n;
      const double c2_s = c2_n;
      const double c3_s = c3_n;
      const double c4_s = c4_n;

      const double fold_n1 = vectorDot3D(f_n1,n);
      const double fold_n2 = vectorDot3D(f_n2,n);
      const double fold_n3 = vectorDot3D(f_n3,n);
      const double fold_n4 = vectorDot3D(f_n4,n);

      const double n_0[3] = {r_0[0]/d_0,r_0[1]/d_0,r_0[2]/d_0};

      double delta_u[3] = {0., 0., 0.};
      vectorSubtract3D(r, r_old, delta_u);

      double delta_u0[3] = {0., 0., 0.};
      vectorSubtract3D(r, r_0, delta_u0);

      const double dn = delta_u[0]*n[0]+delta_u[1]*n[1]+delta_u[2]*n[2];

      const double fnew_n1 = (dn/dt + fold_n1*(1./(k1_n*dt)-0.5/c1_n))/(1./(k1_n*dt)+0.5/c1_n);
      const double fnew_n2 = (dn/dt + fold_n2*(1./(k2_n*dt)-0.5/c2_n))/(1./(k2_n*dt)+0.5/c2_n);
      const double fnew_n3 = (dn/dt + fold_n3*(1./(k3_n*dt)-0.5/c3_n))/(1./(k3_n*dt)+0.5/c3_n);
      const double fnew_n4 = (dn/dt + fold_n4*(1./(k4_n*dt)-0.5/c4_n))/(1./(k4_n*dt)+0.5/c4_n);

      double sh[3] = {0., 0., 0.};
      vectorCross3D(n, n_0, sh);

      double shear[3] = {0., 0., 0.};
      vectorCross3D(n, sh, shear);

      const double s_abs = vectorMag3D(shear);

      double s[3] = {0., 0., 0.};

      if (s_abs > 1.e-20*d_0)
      {
          s[0] = shear[0]/s_abs;
          s[1] = shear[1]/s_abs;
          s[2] = shear[2]/s_abs;
      }

      const double ds = delta_u[0]*s[0]+delta_u[1]*s[1]+delta_u[2]*s[2];
      const double ds_0 = delta_u0[0]*s[0]+delta_u0[1]*s[1]+delta_u0[2]*s[2];

      const double fold_s1 = vectorDot3D(f_s1,s);
      const double fold_s2 = vectorDot3D(f_s2,s);
      const double fold_s3 = vectorDot3D(f_s3,s);
      const double fold_s4 = vectorDot3D(f_s4,s);

      const double fnew_s1 = (ds/dt + fold_s1*(1./(k1_s*dt)-0.5/c1_s))/(1./(k1_s*dt)+0.5/c1_s);
      const double fnew_s2 = (ds/dt + fold_s2*(1./(k2_s*dt)-0.5/c2_s))/(1./(k2_s*dt)+0.5/c2_s);
      const double fnew_s3 = (ds/dt + fold_s3*(1./(k3_s*dt)-0.5/c3_s))/(1./(k3_s*dt)+0.5/c3_s);
      const double fnew_s4 = (ds/dt + fold_s4*(1./(k4_s*dt)-0.5/c4_s))/(1./(k4_s*dt)+0.5/c4_s);

      hist_arr[4] = r[0];
      hist_arr[5] = r[1];
      hist_arr[6] = r[2];

      hist_arr[7] = fnew_n1*n[0];
      hist_arr[8] = fnew_n1*n[1];
      hist_arr[9] = fnew_n1*n[2];
      hist_arr[10] = fnew_n2*n[0];
      hist_arr[11] = fnew_n2*n[1];
      hist_arr[12] = fnew_n2*n[2];
      hist_arr[13] = fnew_n3*n[0];
      hist_arr[14] = fnew_n3*n[1];
      hist_arr[15] = fnew_n3*n[2];
      hist_arr[16] = fnew_n4*n[0];
      hist_arr[17] = fnew_n4*n[1];
      hist_arr[18] = fnew_n4*n[2];

      hist_arr[19] = fnew_s1*s[0];
      hist_arr[20] = fnew_s1*s[1];
      hist_arr[21] = fnew_s1*s[2];
      hist_arr[22] = fnew_s2*s[0];
      hist_arr[23] = fnew_s2*s[1];
      hist_arr[24] = fnew_s2*s[2];
      hist_arr[25] = fnew_s3*s[0];
      hist_arr[26] = fnew_s3*s[1];
      hist_arr[27] = fnew_s3*s[2];
      hist_arr[28] = fnew_s4*s[0];
      hist_arr[29] = fnew_s4*s[1];
      hist_arr[30] = fnew_s4*s[2];

      double f_nk[3] = {0., 0., 0.};
      f_nk[0] = kk_n*(dist-d_0)*n[0];
      f_nk[1] = kk_n*(dist-d_0)*n[1];
      f_nk[2] = kk_n*(dist-d_0)*n[2];

      double f_sk[3] = {0., 0., 0.};
      f_sk[0] = kk_s*ds_0*s[0];
      f_sk[1] = kk_s*ds_0*s[1];
      f_sk[2] = kk_s*ds_0*s[2];

      const double fn_tot = fnew_n1+fnew_n2+fnew_n3+fnew_n4;
      const double fs_tot = fnew_s1+fnew_s2+fnew_s3+fnew_s4;

      //i_forces.delta_F[0] += 0.;
      //i_forces.delta_F[1] += 0.;
      //i_forces.delta_F[2] += 0.;

      i_forces.delta_F[0] += -(fn_tot*n[0]+fs_tot*s[0]+f_nk[0]+f_sk[0]);
      i_forces.delta_F[1] += -(fn_tot*n[1]+fs_tot*s[1]+f_nk[1]+f_sk[1]);
      i_forces.delta_F[2] += -(fn_tot*n[2]+fs_tot*s[2]+f_nk[2]+f_sk[2]);

      //fprintf(screen, "debug i_forces: %f \n", i_forces.delta_F[2]);
      //fprintf(screen, "debug n: %f \n", n[2]);

      j_forces.delta_F[0] += -i_forces.delta_F[0];
      j_forces.delta_F[1] += -i_forces.delta_F[1];
      j_forces.delta_F[2] += -i_forces.delta_F[2];
    }

    void surfacesClose(SurfacesCloseData & scdata, ForceData & i_forces, ForceData & j_forces)
    {
        if(scdata.contact_flags) *scdata.contact_flags &= ~CONTACT_COHESION_MODEL;

        const double dist = sqrt(scdata.rsq);

        class FixPropertyGlobal* max_dist=static_cast<FixPropertyGlobal*>(modify->find_fix_property("maxDistance","property/global","scalar",0,0,force->pair_style));

        if (dist >= max_dist->compute_scalar() || scdata.is_wall)
            return;

        double * hist_arr = &scdata.contact_history[history];

        class FixPropertyGlobal* tsCreate=static_cast<FixPropertyGlobal*>(modify->find_fix_property("timestepCreateBond","property/global","scalar",0,0,force->pair_style));
        class FixPropertyGlobal* radMult=static_cast<FixPropertyGlobal*>(modify->find_fix_property("radiusMultiplier","property/global","scalar",0,0,force->pair_style));

        const int ts_create = tsCreate->compute_scalar();
        const double lambda = radMult->compute_scalar();

        const double dt = update->dt;

        const int itype = scdata.itype;
        const int jtype = scdata.jtype;

        const double ri = scdata.radi;
        const double rj = scdata.radj;

        const double r[3] = {scdata.delta[0], scdata.delta[1], scdata.delta[2]};
        const double r_0[3] = {hist_arr[1],hist_arr[2],hist_arr[3]};
        double r_old[3] = {hist_arr[4],hist_arr[5],hist_arr[6]};

        const double d_0 = sqrt(r_0[0]*r_0[0]+r_0[1]*r_0[1]+r_0[2]*r_0[2]);
        const double d_old = sqrt(r_old[0]*r_old[0]+r_old[1]*r_old[1]+r_old[2]*r_old[2]);

        const double n[3] = {r[0]/dist,r[1]/dist,r[2]/dist};

        if (hist_arr[0] < 0.5)
        {
            const bool create_now = (ts_create < 0 || update->ntimestep == ts_create);
            const double dist_create = lambda*(ri+rj);
            const bool is_close = dist_create >= dist;

            if (create_now && is_close)
            {
                createBond(scdata, r);
                return;
            }
            else
                return;
        }

        double f_n1[3] = {hist_arr[7],hist_arr[8],hist_arr[9]};
        double f_n2[3] = {hist_arr[10],hist_arr[11],hist_arr[12]};
        double f_n3[3] = {hist_arr[13],hist_arr[14],hist_arr[15]};
        double f_n4[3] = {hist_arr[16],hist_arr[17],hist_arr[18]};

        double f_s1[3] = {hist_arr[19],hist_arr[20],hist_arr[21]};
        double f_s2[3] = {hist_arr[22],hist_arr[23],hist_arr[24]};
        double f_s3[3] = {hist_arr[25],hist_arr[26],hist_arr[27]};
        double f_s4[3] = {hist_arr[28],hist_arr[29],hist_arr[30]};

        const double rb = std::min(ri,rj);
        const double Ab = M_PI*rb*rb;

        class FixPropertyGlobal* K_1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("yMax1","property/global","scalar",0,0,force->pair_style));
        class FixPropertyGlobal* K_2=static_cast<FixPropertyGlobal*>(modify->find_fix_property("yMax2","property/global","scalar",0,0,force->pair_style));
        class FixPropertyGlobal* K_3=static_cast<FixPropertyGlobal*>(modify->find_fix_property("yMax3","property/global","scalar",0,0,force->pair_style));
        class FixPropertyGlobal* K_4=static_cast<FixPropertyGlobal*>(modify->find_fix_property("yMax4","property/global","scalar",0,0,force->pair_style));
        class FixPropertyGlobal* K_k=static_cast<FixPropertyGlobal*>(modify->find_fix_property("ySpring","property/global","scalar",0,0,force->pair_style));
        class FixPropertyGlobal* C_1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("muMax1","property/global","scalar",0,0,force->pair_style));
        class FixPropertyGlobal* C_2=static_cast<FixPropertyGlobal*>(modify->find_fix_property("muMax2","property/global","scalar",0,0,force->pair_style));
        class FixPropertyGlobal* C_3=static_cast<FixPropertyGlobal*>(modify->find_fix_property("muMax3","property/global","scalar",0,0,force->pair_style));
        class FixPropertyGlobal* C_4=static_cast<FixPropertyGlobal*>(modify->find_fix_property("muMax4","property/global","scalar",0,0,force->pair_style));

        const double k1_n = Ab/d_0*K_1->compute_scalar();
        const double k2_n = Ab/d_0*K_2->compute_scalar();
        const double k3_n = Ab/d_0*K_3->compute_scalar();
        const double k4_n = Ab/d_0*K_4->compute_scalar();
        const double kk_n = Ab/d_0*K_k->compute_scalar();
        const double c1_n = Ab/d_0*C_1->compute_scalar();
        const double c2_n = Ab/d_0*C_2->compute_scalar();
        const double c3_n = Ab/d_0*C_3->compute_scalar();
        const double c4_n = Ab/d_0*C_4->compute_scalar();

        const double k1_s = k1_n;
        const double k2_s = k2_n;
        const double k3_s = k3_n;
        const double k4_s = k4_n;
        const double kk_s = kk_n;
        const double c1_s = c1_n;
        const double c2_s = c2_n;
        const double c3_s = c3_n;
        const double c4_s = c4_n;

        const double fold_n1 = vectorDot3D(f_n1,n);
        const double fold_n2 = vectorDot3D(f_n2,n);
        const double fold_n3 = vectorDot3D(f_n3,n);
        const double fold_n4 = vectorDot3D(f_n4,n);

        const double n_0[3] = {r_0[0]/d_0,r_0[1]/d_0,r_0[2]/d_0};

        double delta_u[3] = {0., 0., 0.};
        vectorSubtract3D(r, r_old, delta_u);

        double delta_u0[3] = {0., 0., 0.};
        vectorSubtract3D(r, r_0, delta_u0);

        const double dn = delta_u[0]*n[0]+delta_u[1]*n[1]+delta_u[2]*n[2];

        const double fnew_n1 = (dn/dt + fold_n1*(1./(k1_n*dt)-0.5/c1_n))/(1./(k1_n*dt)+0.5/c1_n);
        const double fnew_n2 = (dn/dt + fold_n2*(1./(k2_n*dt)-0.5/c2_n))/(1./(k2_n*dt)+0.5/c2_n);
        const double fnew_n3 = (dn/dt + fold_n3*(1./(k3_n*dt)-0.5/c3_n))/(1./(k3_n*dt)+0.5/c3_n);
        const double fnew_n4 = (dn/dt + fold_n4*(1./(k4_n*dt)-0.5/c4_n))/(1./(k4_n*dt)+0.5/c4_n);

        double sh[3] = {0., 0., 0.};
        vectorCross3D(n, n_0, sh);

        double shear[3] = {0., 0., 0.};
        vectorCross3D(n, sh, shear);

        const double s_abs = vectorMag3D(shear);

        double s[3] = {0., 0., 0.};

        if (s_abs > 1.e-20*d_0)
        {
            s[0] = shear[0]/s_abs;
            s[1] = shear[1]/s_abs;
            s[2] = shear[2]/s_abs;
        }

        const double ds = delta_u[0]*s[0]+delta_u[1]*s[1]+delta_u[2]*s[2];
        const double ds_0 = delta_u0[0]*s[0]+delta_u0[1]*s[1]+delta_u0[2]*s[2];

        const double fold_s1 = vectorDot3D(f_s1,s);
        const double fold_s2 = vectorDot3D(f_s2,s);
        const double fold_s3 = vectorDot3D(f_s3,s);
        const double fold_s4 = vectorDot3D(f_s4,s);

        const double fnew_s1 = (ds/dt + fold_s1*(1./(k1_s*dt)-0.5/c1_s))/(1./(k1_s*dt)+0.5/c1_s);
        const double fnew_s2 = (ds/dt + fold_s2*(1./(k2_s*dt)-0.5/c2_s))/(1./(k2_s*dt)+0.5/c2_s);
        const double fnew_s3 = (ds/dt + fold_s3*(1./(k3_s*dt)-0.5/c3_s))/(1./(k3_s*dt)+0.5/c3_s);
        const double fnew_s4 = (ds/dt + fold_s4*(1./(k4_s*dt)-0.5/c4_s))/(1./(k4_s*dt)+0.5/c4_s);

        hist_arr[4] = r[0];
        hist_arr[5] = r[1];
        hist_arr[6] = r[2];

        hist_arr[7] = fnew_n1*n[0];
        hist_arr[8] = fnew_n1*n[1];
        hist_arr[9] = fnew_n1*n[2];
        hist_arr[10] = fnew_n2*n[0];
        hist_arr[11] = fnew_n2*n[1];
        hist_arr[12] = fnew_n2*n[2];
        hist_arr[13] = fnew_n3*n[0];
        hist_arr[14] = fnew_n3*n[1];
        hist_arr[15] = fnew_n3*n[2];
        hist_arr[16] = fnew_n4*n[0];
        hist_arr[17] = fnew_n4*n[1];
        hist_arr[18] = fnew_n4*n[2];

        hist_arr[19] = fnew_s1*s[0];
        hist_arr[20] = fnew_s1*s[1];
        hist_arr[21] = fnew_s1*s[2];
        hist_arr[22] = fnew_s2*s[0];
        hist_arr[23] = fnew_s2*s[1];
        hist_arr[24] = fnew_s2*s[2];
        hist_arr[25] = fnew_s3*s[0];
        hist_arr[26] = fnew_s3*s[1];
        hist_arr[27] = fnew_s3*s[2];
        hist_arr[28] = fnew_s4*s[0];
        hist_arr[29] = fnew_s4*s[1];
        hist_arr[30] = fnew_s4*s[2];

        double f_nk[3] = {0., 0., 0.};
        f_nk[0] = kk_n*(dist-d_0)*n[0];
        f_nk[1] = kk_n*(dist-d_0)*n[1];
        f_nk[2] = kk_n*(dist-d_0)*n[2];

        double f_sk[3] = {0., 0., 0.};
        f_sk[0] = kk_s*ds_0*s[0];
        f_sk[1] = kk_s*ds_0*s[1];
        f_sk[2] = kk_s*ds_0*s[2];

        const double fn_tot = fnew_n1+fnew_n2+fnew_n3+fnew_n4;
        const double fs_tot = fnew_s1+fnew_s2+fnew_s3+fnew_s4;

        scdata.has_force_update = true;

        //i_forces.delta_F[0] += 0.;
        //i_forces.delta_F[1] += 0.;
        //i_forces.delta_F[2] += 0.;

        int nts = 0.;

        if (dist<=(ri+rj))
        {
            nts = update->ntimestep;
            fprintf(screen, "cont_ts = %d\n", nts);
        }

        i_forces.delta_F[0] += -(fn_tot*n[0]+fs_tot*s[0]+f_nk[0]+f_sk[0]);
        i_forces.delta_F[1] += -(fn_tot*n[1]+fs_tot*s[1]+f_nk[1]+f_sk[1]);
        i_forces.delta_F[2] += -(fn_tot*n[2]+fs_tot*s[2]+f_nk[2]+f_sk[2]);

        if (update->ntimestep == nts || update->ntimestep == nts+1)
            fprintf(screen, "fz = %f \n", i_forces.delta_F[2]);

        //fprintf(screen, "debug i_forces: %f \n", i_forces.delta_F[2]);
        //fprintf(screen, "debug n: %f \n", n[2]);

        j_forces.delta_F[0] += -i_forces.delta_F[0];
        j_forces.delta_F[1] += -i_forces.delta_F[1];
        j_forces.delta_F[2] += -i_forces.delta_F[2];
    }

    private:
        inline void createBond(SurfacesCloseData & scdata, const double *delta)
        {
            double * hist_arr = &scdata.contact_history[history];

            hist_arr[0] = 1.;
            hist_arr[1] = delta[0];
            hist_arr[2] = delta[1];
            hist_arr[3] = delta[2];
            hist_arr[4] = delta[0];
            hist_arr[5] = delta[1];
            hist_arr[6] = delta[2];

            fprintf(screen, "bond created\n");
        }

  protected:

    int history;

    bool hertz_on_off;
  };
}
}
#endif // NORMAL_MODEL_LUDING_H_
#endif
