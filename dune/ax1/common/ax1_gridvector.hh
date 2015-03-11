#ifndef GRIDVECTOR_HH
#define GRIDVECTOR_HH

#include<math.h>
#include<vector>

/**
 * Class GridVector, modified version from Dominic Kempf
 */
template<typename ctype>
class GridVector : public std::vector<ctype>
{
  public:

  void equidistant_n_h(int n,ctype h)
  {
    if (this->empty()) 
      this->push_back(0.0);
    for (int i=0; i<n; i++)
      this->push_back(this->back()+h);
  }

  void equidistant_n_xend(int n, ctype xend)
  {
    if (this->empty())
      this->push_back(0.0);
    ctype h = (xend-this->back())/n;
    for (int i=0; i<n-1; i++)
      this->push_back(this->back()+h);
    this->push_back(xend);
  }
  void equidistant_h_xend(ctype h, ctype xend)
  {
    if (this->empty()) 
      this->push_back(0.0);
    while (this->back() < xend)
      this->push_back(this->back()+h);
  }
  void linear_n_h0_dh(int n, ctype h0, ctype dh)
  {
    if (this->empty())
      this->push_back(0.0);
    for (int i=0; i<n; i++)
      this->push_back(this->back()+h0+i*dh);
  }
  void linear_h0_dh_xend(ctype h0, ctype dh, ctype xend)
  {
    if (this->empty())
      this->push_back(0.0);
    ctype hi=h0;
    while (this->back() < xend)
    {
      this->push_back(this->back()+hi);
      hi += dh;
    }
  }
  void linear_n_h0_xend(int n, ctype h0, ctype xend)
  {
    if (this->empty())
      this->push_back(0.0);
    ctype dh = 2*(xend-this->back()-n*h0)/(n*(n-1));
    for (int i=0; i<n-1; i++)
      this->push_back(this->back()+h0+i*dh);
    this->push_back(xend);
  }
  void linear_h0_hend_xend(ctype h0, ctype hend, ctype xend)
  {
    if (this->empty())
      this->push_back(0.0);
    int n = (int)(2*(xend-this->back())/(h0+hend));
    this->linear_n_h0_xend(n,h0,xend);
  }
  void geometric_n_h0_q(int n, ctype h0, ctype q)
  {
    if (this->empty())
      this->push_back(0.0);
    ctype h = h0;
    for (int i=0; i<n; i++)
    {
      this->push_back(this->back()+h);
      h *= q;
    }
  }
  void geometric_h0_q_xend(ctype h0, ctype q, ctype xend)
  {
    if (this->empty())
      this->push_back(0.0);
    ctype h=h0;
    while (this->back() < xend)
    {
      this->push_back(this->back()+h);
      h *= q;
    }
  }
  void geometric_n_q_xend(int n, ctype q, ctype xend)
  {
    if (this->empty())
      this->push_back(0.0);
    ctype h = (xend-this->back())*(1-q)/(1-pow(q,n));
    for (int i=0; i<n-1; i++)
    {
      this->push_back(this->back()+h);
      h *= q;
    }
    this->push_back(xend);
  }

  ctype newton(int n,ctype x_s,ctype x_e,ctype h)
  {
    ctype m = (x_e - x_s)/h;
    ctype xold = 0.0;
    ctype xnew = x_e-x_s;
    while (fabs(xnew-xold) > 1E-8)
    {
      xold = xnew;
      xnew = xold - (-pow(xold,n)+m*xold-m+1)/(-n*pow(xold,n-1)+m);
    }
    if (fabs(xnew-1) < 1E-6)
    {
      xold = x_e-x_s;
      xnew = 0.0;
      while (fabs(xnew-xold) > 1E-8)
      {
        xold = xnew;
        xnew = xold - (-pow(xold,n)+m*xold-m+1)/(-n*pow(xold,n-1)+m);
      }
    }
    return xnew;
  }

  void geometric_n_h0_xend(int n, ctype h0, ctype xend)
  {
    if (this->empty())
      this->push_back(0.0);
    ctype q = newton(n,this->back(),xend,h0);
    ctype h = h0;
    for (int i=0; i<n-1; i++)
    {
      this->push_back(this->back()+h);
      h *= q;
    }
    this->push_back(xend);
  }

  void geometric_n_hend_xend(int n, ctype hend, ctype xend)
  {
    if (this->empty())
      this->push_back(0.0);
    ctype p = newton(n,this->back(),xend,hend);
    ctype h = hend*pow(p,n-1);
    for (int i=0; i<n-1; i++)
    {
      this->push_back(this->back()+h);
      h /= p;
    }
    this->push_back(xend);
  }

  // Optional parameter guarantee (0=h0, 1=hend) to make sure either h0 or hend is guaranteed
  void geometric_h0_hend_xend(ctype h0, ctype hend, ctype xend, int guarantee = 0)
  {
    if (this->empty())
      this->push_back(0.0);
    ctype q = (xend-this->back()-h0)/(xend-this->back()-hend);
    int n = (int)(log(hend/h0)/log(q))+1;
    if(guarantee == 0)
      this->geometric_n_h0_xend(n,h0,xend);
    else
      this->geometric_n_hend_xend(n,hend,xend);
  }
  void start(ctype x)
  {
    if (this->empty())
      this->push_back(x);
  }
  void setEntityChange()
  {
    entityChange.push_back(this->size()-1);
  }

  std::vector<int> entityChange;
};

/**
 * A specialization for a special form of grid where
 * one coordinate axis represents the radial axis of
 * a cylinder. A 2D tensor grid is transformed to a
 * 3D cylinder by assuming symmetry in the angular
 * direction.
 * This class represents the radial direction and
 * provides methods for grid generation for this
 * very special case.
 */
template<typename ctype>
class CylinderGridVector : public GridVector<ctype>
{
  public:
    typedef GridVector<ctype> BaseT;

    CylinderGridVector(ctype H_)
    : H(H_)
    {}

    /**
     * Fille the range until rend with intervals such that the
     * volumes of all cells will be equal
     *
     * @param h This is the start/end interval, depending on whether guarantee is 0 or 1
     * @param q This is the factor by which the volume is to be increased in each interval
     * @param rend The end point of the interval
     * @param allowExceedRange As only one range border from the interval [r0 rend] can
     * be guaranteed, this flag tells wether the range may be exceeded at the other end
     * or if it should lie inside the range in every case, even if it wasn't fully exhausted
     * @param guarantee 0=guarantee lower border r0; 1=guarantee upper border rend
     */
    void geometric_volume_h_q_rend(ctype h, ctype q, ctype rend, bool allowExceedRange = true,
        int guarantee = 0, bool verbose = true)
    {
      if (this->empty())
        this->push_back(0.0);

      ctype r = this->back();

      // Swap limits in case the upper limit is to be guaranteed
      if(guarantee == 1)
      {
        r = -rend;
        rend = -this->back();
        this->push_back(-r);
      }

      // Lower limit for h
      ctype h_limit = 1;

      // Start with equidistant h
      while(rend-r > 0)
      {
        // Do not add the last interval if it would exceed the range [r0 rend]
        if(! allowExceedRange && (r+h > rend-1e-6))
        {
          if(verbose)
          {
            debug_jochen << "Next step: r =" << r << ", dr=" << h
              << " -> R = " << (r+h)
              << " would exceed interval border, omitting this point!" << std::endl;
          }
          break;
        }

        double volume2 = con_pi*(r+h)*(r+h)*H - con_pi*r*r*H;
        double volume = (2*con_pi*H*r + con_pi*H*h)*h;

        if(verbose)
        {
          debug_jochen << "Current step: r =" << r << ", dr=" << h
              << " -> R = " << (r+h) << std::endl;
          debug_jochen << "   V1 = "  << volume << std::endl;
          debug_jochen << "   V2 = "  << volume2 << std::endl;
        }

        this->push_back(std::abs(r + h));

        // Solve quadratic equation with p-q-formula
        ctype radicand = r*r + (1+q)*2*r*h + (1+q)*h*h;
        if(radicand < 0)
        {
          debug_warn << "Could not calculated square root of negative radicand!" << std::endl;
          break;
        }
        ctype h1 = -(r+h)+std::sqrt(radicand);
        ctype h2 = -(r+h)-std::sqrt(radicand);
        if(verbose)
        {
          debug_jochen << "  h1 = " << h1 << std::endl;
          debug_jochen << "  h2 = " << h2 << std::endl;
        }

        ctype hnew;

        if(guarantee == 0)
          hnew = h1;
        else
          hnew = h2;

        // Update r before updating h
        r = r + h;
        h = std::max(hnew, h_limit);
      }

      std::sort(this->begin(), this->end());
    }

    /**
    * Fill the range until rend with intervals such that the
    * volumes of all cells will be equal
    */
   void equivolume_h_rend(ctype h, ctype rend, bool allowExceedRange = true,
       int guarantee = 0, bool verbose = true)
   {
     geometric_volume_h_q_rend(h, 1.0, rend, allowExceedRange, guarantee, verbose);
   }


   // The 'simulate' flag just calculates h0 and returns it, without actually doing
   // anything with the vector
    void equivolume_n_rend(int nPoints, ctype rend, bool allowExceedRange = true,
        int guarantee = 0, bool verbose = true)
    {
      // Calculate h first
      ctype r0 = this->back();

      ctype remaining_vol = remaining_volume(rend);

      ctype volume_per_element = remaining_vol / nPoints;

      if(verbose)
      {
        debug_jochen << "remaining volume: " << remaining_vol << std::endl;
        debug_jochen << "volume per element: " << volume_per_element << std::endl;
      }

      ctype h;
      ctype volume;
      if(guarantee == 0)
      {
        h = get_h_for_volume(volume_per_element);

        volume = con_pi*(r0+h)*(r0+h)*H
                  - con_pi*r0*r0*H;
      } else {
        h = get_h_for_volume(volume_per_element, false, rend);
        volume = con_pi*rend*rend*H
               - con_pi*(rend-h)*(rend-h)*H;
      }

      if(verbose)
      {
        debug_jochen << "calculated h0=" << h << std::endl;
        debug_jochen << " --> volume=" << volume << std::endl;
      }

      equivolume_h_rend(h, rend, allowExceedRange, guarantee, verbose);
    }



    ctype remaining_volume(ctype rend)
    {
      return volume_r0_rend(this->back(), rend);
    }

    ctype next_volume(ctype h)
    {
      assert(! this->empty());

      return volume_r0_h(this->back(), h);
    }

    ctype volume_r0_rend(ctype r0, ctype rend)
    {
      return (con_pi*rend*rend*H-con_pi*r0*r0*H);
    }

    ctype volume_r0_h(ctype r0, ctype h)
    {
      return (2*con_pi*H*r0 + con_pi*H*h)*h;
    }

    ctype get_h_for_volume(ctype V, bool up = true, ctype r = std::numeric_limits<ctype>::max())
    {
      if(r == std::numeric_limits<ctype>::max())
        r = this->back();

      if(! up)
      {
        r *= -1;
        V *= -1;
      }

      // Solve quadratic equation with p-q-formula
      ctype radicand = r*r + V/(con_pi*H);
      if(radicand < 0)
        DUNE_THROW(Dune::Exception, "Could not calculate square root of negative radicand!");

      ctype h1 = -r+std::sqrt(radicand);
      ctype h2 = -r-std::sqrt(radicand);

      debug_jochen << "h1 = " << h1 << std::endl;
      debug_jochen << "h2 = " << h2 << std::endl;

      if(up)
        return h1;
      else
        return h2;
    }

  private:
    ctype H; // 'height' of the cylinder
};

#endif
