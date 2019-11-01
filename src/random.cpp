#include "random.h"


RandomNumbers::RandomNumbers(unsigned long int s)
	: seed(s)
{
	if (seed == 0) {
        std::random_device rd;
        seed = rd();
    }
    rng = std::mt19937(seed);
}

void RandomNumbers::uniform_double(std::vector<double>& tab, double lower, double upper){
	//on pourrait utiliser la surcharge uniform_double, mais ca créerai un nouveau uniform_real_... à chaque iteration sur la boucle
	std::uniform_real_distribution<> dis(lower, upper);
		for(size_t i(0); i < tab.size(); ++i){
			tab[i] = dis(rng);
		}
}
	
double RandomNumbers::uniform_double(double lower, double upper){
	std::uniform_real_distribution<> dis(lower, upper);
	return dis(rng);
	}
	
void RandomNumbers::normal(std::vector<double>& tab, double mean, double sd){
	std::normal_distribution<> dis{mean, sd};
		for(size_t i(0); i < tab.size(); ++i){
			tab[i] = dis(rng);
		}
	}
	
double RandomNumbers::normal(double mean, double sd){
	std::normal_distribution<> dis{mean, sd};
	return dis(rng);
	}
	
void RandomNumbers::poisson(std::vector<int>& tab, double mean){
	std::poisson_distribution<int> distribution(mean);
	for(size_t i(0); i < tab.size(); ++i){
			tab[i] = distribution(rng);
		}
	}
	
int RandomNumbers::poisson(double mean){
	std::poisson_distribution<int> distribution(mean);
	return distribution(rng);
	}
