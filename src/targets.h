#include <Eigen/Core>

typedef const Eigen::Ref<const Eigen::MatrixXd> constRefMat;
typedef const Eigen::Ref<const Eigen::VectorXd> constRefVec;

class Target
{
public:
    virtual ~Target(){};

    virtual int64_t NumDensityEvaluations() const = 0;
    virtual int64_t NumGradientEvaluations() const = 0;

    virtual double LogDensity(constRefVec &theta) const = 0;
    virtual Eigen::VectorXd GradLogDensity(constRefVec &theta) const = 0;
    virtual Eigen::MatrixXd HessLogDensity(constRefVec &theta) const = 0;
};


// Logistic regression target with spherical Gaussian prior N(0, lambda^2)
class LogisticRegression : public Target
{
public:
  LogisticRegression(constRefMat &yX, const double &lambda)
        : yX_(yX), lambda_(lambda), num_density_evaluations_(0), num_gradient_evaluations_(0)
    {
    }
    // Density evaluation
    virtual double LogDensity(constRefVec &theta) const 
    {
        num_density_evaluations_ += 1;
        return (yX_ * theta).unaryExpr(&logcdf).sum() - 0.5 * theta.squaredNorm() / (lambda_ * lambda_);
    }
    virtual int64_t NumDensityEvaluations() const
    {
        return num_density_evaluations_;
    }

    // Gradient evaluation
    virtual Eigen::VectorXd GradLogDensity(constRefVec &theta) const 
    {
        num_gradient_evaluations_ += 1;
        Eigen::VectorXd weight = (yX_ * theta).unaryExpr(&logcdfdash);
        Eigen::VectorXd gradloglike = (weight.asDiagonal() * yX_).colwise().sum();
        return gradloglike - theta / (lambda_ * lambda_);
    }
    virtual int64_t NumGradientEvaluations() const
    {
        return num_gradient_evaluations_;
    }

    // Hessian evaluation
    virtual Eigen::MatrixXd HessLogLikelihood(const int &n, constRefVec &theta) const
    {
        double weight = logcdfdoubledash(yX_.row(n).dot(theta));
        return (weight * yX_.row(n).transpose()) * yX_.row(n);
    }
    virtual Eigen::MatrixXd HessLogLikelihood(constRefVec &theta) const
    {
        Eigen::MatrixXd hess = Eigen::MatrixXd::Zero(yX_.cols(), yX_.cols());
        for (int n = 0; n < yX_.rows(); n++) hess += HessLogLikelihood(n, theta);
        return hess;
    }
    virtual Eigen::MatrixXd HessLogDensity(constRefVec &theta) const 
    {
        return HessLogLikelihood(theta) - (1 / lambda_) * Eigen::MatrixXd::Identity(yX_.cols(),yX_.cols());
    }

private:
    constRefMat yX_;
    const double lambda_;
    mutable int64_t num_density_evaluations_;
    mutable int64_t num_gradient_evaluations_;

  //log-sigmoid
  static double logcdf(const double &t)
  {
    if (t < -33.3)
    {
      return t;
    }
    else if (t <= -18)
    {
      return t - exp(t);
    }
    else if (t <= 37)
    {
      return -log1p(exp(-t));
    }
    else
    {
      return -exp(-t);
    }
  }
  
  //log-sigmoid dashed
  static double logcdfdash(const double &t) 
  {
    if (t < 0)
    {
      return 1. / (1. + exp(t));
    }
    else
    {
      return exp(-t) / (1. + exp(-t));
    }
  }
  
  // log-sigmoid double-dashed
  static double logcdfdoubledash(const double &t) 
  {
    if (t < 0)
    {
      return -exp(t) / pow(1. + exp(t), 2);
    }
    else
    {
      return -exp(-t) / pow(1. + exp(-t), 2);
    }
  }

};



