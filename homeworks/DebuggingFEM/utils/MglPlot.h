#ifndef MGL_PLOT_H
#define MGL_PLOT_H

#include <mgl2/mgl.h>

namespace mgl {

class MglPlot {
public:

  MglPlot(const std::string & style)
    : style_{style}
    , legend_{""}
  {}
  virtual void plot(mglGraph* gr) = 0;
  virtual bool is_3d() = 0;

  MglPlot& label(const std::string &l) {
    legend_ = l;
    return *this;
  }

  MglPlot& style(const std::string &s) {
    style_ = s;
    return *this;
  }

  MglPlot& width(int w) {
    w = w < 0 ? 0 : (w > 9 ? 9 : w);
    if(style_.size() == 3)
      style_[3] = char( '0' + w );
    else
      style_ += {char( '0' + w )};
    return *this;
  }

protected:
  std::string style_;
  std::string legend_;
};

class MglPlot2d : public MglPlot {
public:

  MglPlot2d(const mglData& xd, const mglData& yd, const std::string &style)
    : MglPlot(style)
    , xd_(xd)
    , yd_(yd)
  {}

  void plot(mglGraph * gr) {
    gr->Plot(xd_, yd_, style_.c_str());
    gr->AddLegend(legend_.c_str(), style_.c_str());
  }

  bool is_3d() {
    return false;
  }

private:
  mglData xd_;
  mglData yd_;
};

class MglPlot3d : public MglPlot {
public:

  MglPlot3d(const mglData& xd, const mglData& yd, const mglData& zd, const std::string& style)
    : MglPlot(style)
    , xd_(xd)
    , yd_(yd)
    , zd_(zd)
  {}

  void plot(mglGraph * gr) {
    gr->Plot(xd_, yd_, zd_, style_.c_str());
    gr->AddLegend(legend_.c_str(), style_.c_str());
  }

  bool is_3d() {
    return true;
  }

private:
  mglData xd_;
  mglData yd_;
  mglData zd_;
};

class MglFPlot : public MglPlot {
public:

  MglFPlot(const std::string & fplot_str, const std::string &style)
    : MglPlot(style)
    , fplot_str_(fplot_str)
  {}

  void plot(mglGraph * gr) {
    gr->FPlot(fplot_str_.c_str(), style_.c_str());
    gr->AddLegend(legend_.c_str(), style_.c_str());
  }

  bool is_3d() {
    return false;
  }

private:
  std::string fplot_str_;
};


} // end namespace
#endif
