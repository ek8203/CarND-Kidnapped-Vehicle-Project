/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles.
  num_particles = 100;

  // Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1.
	// Add random Gaussian noise to each particle.
  // NOTE: Consult particle_filter.h for more information about this method (and others in this file).
  default_random_engine gen;

  // std[] - GPS measurement uncertainty [x [m], y [m], theta [rad]]
  // Creates a normal (Gaussian) distribution for x, y and theta.
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);

  // init particles
  for (int i = 0; i < num_particles; i++) {
    // init a particle
    Particle  particle;
    particle.id       = i;
    particle.x        = dist_x(gen);  // "gen" is the random engine initialized earlier
    particle.y        = dist_y(gen);;
    particle.theta    = dist_theta(gen);;
    particle.weight   = 1.0;

    // add the new particle
    particles.push_back(particle);
  }

  // init wieights to 1
  weights.resize(num_particles, 1.0);

  // done with init
  is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

  // generate random Gaussian noise
  default_random_engine gen;
  normal_distribution<double> dist_x(0, std_pos[0]);
  normal_distribution<double> dist_y(0, std_pos[1]);
  normal_distribution<double> dist_yaw(0, std_pos[2]);

  for (auto particle = particles.begin() ; particle != particles.end(); ++particle) {
    // initial position and heading
    double x = particle->x;
    double y = particle->y;
    double yaw = particle->theta;

    // predict for yaw_rate != zero
    if(fabs(yaw_rate) > 0.0001) {
      // new position and heading mean
      x += (velocity/yaw_rate) * (sin(yaw + yaw_rate*delta_t) - sin(yaw));
      y += (velocity/yaw_rate) * (cos(yaw) - cos(yaw + yaw_rate*delta_t));
      yaw += yaw_rate*delta_t;
    }
    else  {
      x += velocity * delta_t * cos(yaw);
      y += velocity * delta_t * sin(yaw);
    }

    // a new prediction with added noise
    particle->x        = x + dist_x(gen);  // "gen" is the random engine initialized earlier
    particle->y        = y + dist_y(gen);;
    particle->theta    = yaw + dist_yaw(gen);;

    //cout << "Predict" << "\t" << n << "\t" << particle->weight << "\t" << particle->x << "\t" << particle->y << endl;
  }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

  // Transform observations from VEHICL's coordinate system to the particles MAP's coordinate system
  for (int n = 0; n < num_particles; n++) {

    Particle particle = particles[n];

    vector<LandmarkObs> trans_observations; // a vector of transformed observations
    vector<LandmarkObs> lm_in_range;        // a vector of landmarks in the sensor range
    vector<int> associations;
    vector<double> sense_x;
    vector<double> sense_y;

    // find landmarks in the sensor range assuming the sensor is in the position of the particle
    for (auto lm = map_landmarks.landmark_list.begin(); lm != map_landmarks.landmark_list.end(); ++lm)  {

      double range = dist(particle.x, particle.y, lm->x_f, lm->y_f); // always positive range
      if(range <= sensor_range) {
        LandmarkObs lm_obs;
        lm_obs.id = lm->id_i;
        lm_obs.x  = lm->x_f;
        lm_obs.y  = lm->y_f;
        lm_in_range.push_back(lm_obs);
      }
    }

    double particle_weight = 1.0; // initial particle weight

    // transform observations and associate with the closest landmark in range
    for(auto obs = observations.begin(); obs != observations.end(); ++obs)  {
      LandmarkObs trans_obs;
      double  x = obs->x;
      double  y = obs->y;
      double  theta = particle.theta;

      // transformed coordinates
      double  trans_x = particle.x + cos(theta) * x - sin(theta) * y;
      double  trans_y = particle.y + sin(theta) * x + cos(theta) * y;

      // find the nearest neighbor landmark association
      auto lm = lm_in_range.begin();
      auto nearest_lm = lm;;
      double min_range = dist(trans_x, trans_y, lm->x, lm->y);
      // find the closest
      while(lm != lm_in_range.end()) {
        double range = dist(trans_x, trans_y, lm->x, lm->y);
        if(range < min_range) {
          min_range = range;
          nearest_lm = lm;
        }
        lm++;
      }

      trans_obs.x = trans_x;
      trans_obs.y = trans_y;
      trans_observations.push_back(trans_obs);

      sense_x.push_back(trans_x);
      sense_y.push_back(trans_y);
      associations.push_back(nearest_lm->id);

      // Calculate a Multivariate-Gaussian probability density value for each observation
      // Adopted from the forum post:
      // https://discussions.udacity.com/t/transformations-and-associations-calculating-the-particles-final-weight/308602/2
      double std_x = std_landmark[0];
      double std_y = std_landmark[1];
      double multiplier = 1.0/(2*M_PI*std_x*std_y);
      double cov_x = pow(std_x, 2.0);
      double cov_y = pow(std_y, 2.0);

      double observation_prob = multiplier*exp(-pow(trans_x - nearest_lm->x, 2.0)/(2.0*cov_x) -
                                                pow(trans_y - nearest_lm->y, 2.0)/(2.0*cov_y));

      // The final weight is a product of all the calculated measurement probabilities
      particle_weight *= observation_prob;
    }

    // Update weight of the particle
    particle.weight = particle_weight;
    weights[n]      = particle_weight;

    // update associations
    particles[n] = SetAssociations(particle, associations, sense_x, sense_y);

    //cout  << "Update" << "\t" << n << "\t" << particles[n].weight << "\t" << particles[n].x << "\t" << particles[n].y << endl;
  }
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

  // The following resample method is adopted from the forum post example:
  // https://discussions.udacity.com/t/resampling-algorithm-using-resampling-wheel/241313

  // discrete_distribution produces random integers on the interval [0, n), where the probability of each
  // individual integer i is defined as w[i]/Sum, that is the weight of the i-th integer divided by the sum
  // of all n weights.
  discrete_distribution<> dist(weights.begin(), weights.end());
  default_random_engine gen;
  vector<Particle> new_particles; // new vector of the selected particles

  for (int i = 0; i < num_particles; ++i)
  {
    int index = dist(gen);
    new_particles.push_back(particles[index]);
  }

  particles = new_particles;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations = associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;

    return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
